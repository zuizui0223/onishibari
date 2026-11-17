#!/usr/bin/env python3
"""Audit iNaturalist observation/photo availability for Daphne pseudomezereum.

Photo and plant-phenology annotations are volunteer metadata, not sex diagnoses.
No precise coordinates, usernames, or downloaded image bytes are exported.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.inaturalist.org/v1"
TAXON_ID = 447647
PLACE_JP = 6737
HEADERS = {"User-Agent": "onishibari-evidence-audit/0.2 (research citation)"}
BASE = {
    "taxon_id": TAXON_ID,
    "place_id": PLACE_JP,
    "d1": "2016-01-01",
    "d2": "2026-11-14",
    "per_page": 1,
}

def get(path, filters):
    url = API + path + "?" + urlencode(filters, doseq=True)
    last = None
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=40) as f:
                return json.loads(f.read().decode("utf-8")), url
        except (HTTPError, URLError, TimeoutError) as err:
            last = str(err)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"iNaturalist API request failed: {url}, {last}")

def count(filters):
    res, url = get("/observations", dict(BASE, **filters))
    return res.get("total_results", 0), url

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/derived/inat_onishibari_opportunity_snapshot.json")
    parser.add_argument("--sample-size", type=int, default=100)
    args = parser.parse_args()
    if not 0 <= args.sample_size <= 200:
        parser.error("--sample-size must be 0..200")
    tx, taxon_url = get("/taxa/" + str(TAXON_ID), {})
    rows = tx.get("results", [])
    if not rows or "pseudomezereum" not in rows[0].get("name", "").lower():
        raise ValueError(f"Unexpected iNaturalist taxon: {rows}")
    counts = {}
    urls = {"taxon": taxon_url}
    for label, filters in (
        ("all_verifiable_or_not", {}),
        ("photos_any_grade", {"photos": "true"}),
        ("research_photos", {"photos": "true", "quality_grade": "research"}),
        ("flower_annotated", {"photos": "true", "quality_grade": "research",
                              "term_id": 12, "term_value_id": 13}),
        ("fruit_annotated", {"photos": "true", "quality_grade": "research",
                             "term_id": 12, "term_value_id": 14}),
        ("female_annotated", {"photos": "true", "quality_grade": "research",
                              "term_id": 9, "term_value_id": 10}),
        ("male_annotated", {"photos": "true", "quality_grade": "research",
                            "term_id": 9, "term_value_id": 11}),
    ):
        counts[label], urls[label] = count(filters)
    months = {}
    month_urls = {}
    for month in range(1,13):
        months[str(month)], month_urls[str(month)] = count({
            "photos": "true", "quality_grade": "research", "month": month})
        time.sleep(0.25)
    preview, preview_url = get("/observations", dict(
        BASE, photos="true", quality_grade="research",
        per_page=args.sample_size if args.sample_size else 1,
        order_by="observed_on", order="desc"))
    sample = []
    for item in preview.get("results", [])[:args.sample_size]:
        photos = item.get("photos") or []
        sample.append({
            "observation_id": item.get("id"),
            "observed_on": item.get("observed_on"),
            "quality_grade": item.get("quality_grade"),
            "taxon_id": (item.get("taxon") or {}).get("id"),
            "identification": (item.get("taxon") or {}).get("name"),
            "annotated_terms": [{
                "term_id": ann.get("controlled_attribute_id") or (ann.get("controlled_attribute") or {}).get("id"),
                "value_id": ann.get("controlled_value_id") or (ann.get("controlled_value") or {}).get("id"),
            } for ann in item.get("annotations") or []],
            "photo_count": len(photos),
            "photo_refs": [{
                "url": p.get("url"),
                "license_code": p.get("license_code"),
                "attribution": p.get("attribution"),
            } for p in photos[:3]],
            "observation_url": "https://www.inaturalist.org/observations/" + str(item.get("id")),
        })
    male_obs, male_url = get("/observations", dict(BASE, photos="true", quality_grade="research",
                               term_id=9, term_value_id=11, per_page=1))
    doc = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source": "iNaturalist v1 API",
        "taxon_name": rows[0].get("name"),
        "taxon_id": TAXON_ID,
        "japan_place_id": PLACE_JP,
        "base_filters": BASE,
        "counts": counts,
        "research_photo_observations_by_month": months,
        "monthly_counts_sum": sum(months.values()),
        "month_sum_matches_photos": sum(months.values()) == counts["research_photos"],
        "sample": sample,
        "male_annotated_example_ids": [r.get("id") for r in male_obs.get("results", [])],
        "urls": dict(urls, photo_sample=preview_url, male_example=male_url, month_queries=month_urls),
        "limitations": (
            "Each count is observation events, NOT unique plants or flowering. "
            "Only research-grade photo observations are used for month counts. "
            "Flowering annotations can be missing; term 9 male/female may not be offered for plants. "
            "Month is observation effort plus biology. The sample is recent, not random."
        ),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("INAT_AUDIT=" + json.dumps({
        "taxon_id": TAXON_ID, "taxon_name": rows[0].get("name"),
        "counts": counts, "monthly": months,
        "month_sum_valid": doc["month_sum_matches_photos"],
        "sample": len(sample),
    }, sort_keys=True), flush=True)
    print("INAT_AUDIT_OUTPUT=" + str(output), flush=True)

if __name__ == "__main__":
    main()
