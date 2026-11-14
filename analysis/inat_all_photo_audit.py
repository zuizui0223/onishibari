#!/usr/bin/env python3
"""Capture all public Japanese iNaturalist D. pseudomezereum photo opportunities.

Preserves original observation/photo IDs, dates and Creative Commons attribution.
Never infers sex morph, pollination, unique shrubs or flowering phenology.
No coordinates, user handles, or precise place descriptions are exported.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.inaturalist.org/v1"
TAXON = 447647
PLACE = 6737
BASE = {
    "taxon_id": TAXON, "place_id": PLACE,
    "d1": "2016-01-01", "d2": "2026-11-14", "photos": "true"
}
HEADERS = {"User-Agent": "onishibari-evidence-audit/0.3 (academic research)"}
MAX_PER_PAGE = 200

def get(query, path="/observations"):
    url = API + path + "?" + urlencode(query)
    for attempt in range(4):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=40) as stream:
                return json.loads(stream.read().decode("utf-8")), url
        except (HTTPError, URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)

def fetch():
    result = []
    urls = []
    first, url = get(dict(BASE, page=1, per_page=MAX_PER_PAGE, order="asc", order_by="id"))
    urls.append(url)
    total = int(first["total_results"])
    result.extend(first.get("results", []))
    page = 1
    while len(result) < total:
        page += 1
        payload, url = get(dict(BASE, page=page, per_page=MAX_PER_PAGE,
                                order="asc", order_by="id"))
        result.extend(payload.get("results", []))
        urls.append(url)
        if not payload.get("results"):
            raise RuntimeError("Unexpected empty iNaturalist page while records remain")
    if len(result) != total:
        raise ValueError("API results have inconsistent pagination")
    return result, urls, total

def transform(items):
    obs = []
    photo_ids = set()
    for row in items:
        if not row.get("id"):
            raise ValueError("Missing observation ID")
        photos = []
        for p in row.get("photos") or []:
            pid = p.get("id")
            if not pid:
                raise ValueError("Missing photo ID")
            if pid in photo_ids:
                raise ValueError(f"Duplicate image photo_id: {pid}")
            photo_ids.add(pid)
            photos.append({
                "photo_id": pid,
                "url": p.get("url"),
                "license_code": p.get("license_code"),
                "attribution": p.get("attribution"),
                "original_photo_count": len(row.get("photos") or []),
            })
        if not photos:
            raise ValueError("Non-photo record returned from photos=true filter")
        obs.append({
            "observation_id": row["id"],
            "observed_on": row.get("observed_on"),
            "quality_grade": row.get("quality_grade"),
            "taxon_id_at_observation": (row.get("taxon") or {}).get("id"),
            "taxon_name_at_observation": (row.get("taxon") or {}).get("name"),
            "annotations": [{
                "term_id": a.get("controlled_attribute_id") or
                           (a.get("controlled_attribute") or {}).get("id"),
                "value_id": a.get("controlled_value_id") or
                            (a.get("controlled_value") or {}).get("id"),
            } for a in row.get("annotations") or []],
            "photos": photos,
            "url": "https://www.inaturalist.org/observations/" + str(row["id"]),
        })
    return obs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/derived/inat_all_photo_snapshot.json")
    args = parser.parse_args()
    items, urls, total = fetch()
    observations = transform(items)
    photo_count = sum(len(o["photos"]) for o in observations)
    counts = Counter(o["quality_grade"] for o in observations)
    dates = Counter((o["observed_on"] or "unknown")[:7] for o in observations)
    doc = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source": "iNaturalist API /v1/observations",
        "filters": BASE,
        "api_urls": urls,
        "total_observations": total,
        "total_photo_refs": photo_count,
        "quality_grades": dict(counts),
        "observations": observations,
        "limitations": (
            "All public image-bearing observation events, not unique plant IDs. "
            "No independent sex determination, flowering phenology or pollinator visits. "
            "Observation grade can change; source URLs and retrieval timestamp retained. "
            "Photograph reuse must follow its license and attribution."
        ),
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("INAT_ALL_PHOTOS=" + json.dumps({
        "observations": total, "photo_refs": photo_count,
        "grades": dict(counts), "year_month_bins": len(dates),
        "flower_or_sex_verified": False, "output": str(out)
    }, ensure_ascii=False, sort_keys=True), flush=True)

if __name__ == "__main__":
    main()
