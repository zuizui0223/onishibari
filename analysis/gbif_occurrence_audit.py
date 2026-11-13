#!/usr/bin/env python3
"""Audit public GBIF records for *photo opportunities*, not plant sex or pollination.

This script requests live occurrence metadata and writes a dated snapshot.
It cannot infer flower sex, scent, pollinator identity, phenophase or selfing
from the month of a GBIF record.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = "https://api.gbif.org/v1"
SPECIES = [
    "Daphne pseudomezereum",
    "Daphne jezoensis",
    "Daphne kiusiana",
]
HEADERS = {"User-Agent": "onishibari-evidence-audit/0.1 (GitHub research)"}
FILTERS = {
    "country": "JP",
    "occurrenceStatus": "PRESENT",
    "year": "2016,2026",
    "hasCoordinate": "true",
    "mediaType": "StillImage",
}
ALL_BASE = {"country": "JP", "occurrenceStatus": "PRESENT"}

def api_get(path, query):
    url = ROOT + path + "?" + urlencode(query, doseq=True)
    last = None
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers=HEADERS), timeout=40) as response:
                return json.loads(response.read().decode("utf-8")), url
        except (HTTPError, URLError, TimeoutError) as error:
            last = str(error)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"GBIF API request failed after retries: {url}; {last}")

def occ_count(query):
    result, url = api_get("/occurrence/search", dict(query, limit=0))
    return result.get("count", 0), url

def facet_month(query):
    response, url = api_get("/occurrence/search", dict(
        query, limit=0, facet="month", facetLimit=12, facetMincount=1
    ))
    hits = {}
    for facet in response.get("facets", []):
        if str(facet.get("field", "")).lower() == "month":
            for cell in facet.get("counts", []):
                raw = str(cell.get("name", ""))
                if raw.isdigit():
                    hits[int(raw)] = int(cell.get("count", 0))
    # Some API variants omit month facets; explicit queries retain reproducibility.
    if not hits and response.get("count", 0):
        for month in range(1, 13):
            hits[month] = occ_count(dict(query, month=month))[0]
    hits = {str(m): hits.get(m, 0) for m in range(1, 13)}
    return hits, url

def taxon(name):
    match, url = api_get("/species/match", {"name": name, "strict": "true"})
    key = match.get("usageKey")
    epithet = name.split()[-1].lower()
    canonical = str(match.get("canonicalName", "")).lower()
    scientific = str(match.get("scientificName", "")).lower()
    if not key or epithet not in (canonical + " " + scientific):
        raise ValueError(f"Taxon match did not confirm {name}: {match}")
    return {"name_query": name, "gbif_taxon_key": key,
            "matched_name": match.get("scientificName"),
            "canonical_name": match.get("canonicalName"),
            "rank": match.get("rank"),
            "status": match.get("status"), "match_type": match.get("matchType"),
            "lookup_url": url}

def one_species(name, sample_size):
    identity = taxon(name)
    key = identity["gbif_taxon_key"]
    broad = dict(ALL_BASE, taxonKey=key)
    strict = dict(FILTERS, taxonKey=key)
    total_jp, total_url = occ_count(broad)
    georeferenced_jp, coord_url = occ_count(dict(broad, hasCoordinate="true"))
    photo_jp_2016_26, photo_url = occ_count(strict)
    photo_obs_jp_2016_26, human_url = occ_count(
        dict(strict, basisOfRecord="HUMAN_OBSERVATION"))
    months, facet_url = facet_month(strict)
    sample_resp, sample_url = api_get(
        "/occurrence/search", dict(strict, limit=sample_size, offset=0))
    sample = []
    for r in sample_resp.get("results", []):
        media = r.get("media") or []
        url = next((m.get("identifier") for m in media
                    if str(m.get("type", "")).lower() == "stillimage"
                    and m.get("identifier")), None)
        sample.append({
            "gbif_id": r.get("key"),
            "event_date": r.get("eventDate"),
            "month": r.get("month"),
            "basis_of_record": r.get("basisOfRecord"),
            "license": r.get("license"),
            "media_url": url,
            "record_url": "https://www.gbif.org/occurrence/" + str(r.get("key")),
            "dataset_key": r.get("datasetKey"),
            "recorded_scientific_name": r.get("scientificName"),
        })
    monthly_total = sum(months.values())
    return {
        "taxon": identity,
        "counts": {
            "all_jp_present": total_jp,
            "georeferenced_jp_present": georeferenced_jp,
            "jp_georeferenced_photo_2016_2026": photo_jp_2016_26,
            "jp_georeferenced_photo_human_2016_2026": photo_obs_jp_2016_26,
        },
        "monthly_counts_in_photo_subset": months,
        "monthly_sum_matches_count": monthly_total == photo_jp_2016_26,
        "month_01_to_04_fraction": (sum(months[str(i)] for i in range(1, 5))
                                    / monthly_total if monthly_total else None),
        "sample": sample,
        "urls": {
            "total": total_url, "coordinate": coord_url, "photo": photo_url,
            "human_photo": human_url, "month_facet": facet_url,
            "sample": sample_url,
        }
    }

def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", default="data/derived/gbif_photo_opportunity_snapshot.json")
    cli.add_argument("--sample-size", type=int, default=30)
    args = cli.parse_args()
    if not 0 <= args.sample_size <= 100:
        cli.error("sample-size must be between 0 and 100")
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "GBIF Occurrence API v1; live public query",
        "species": [],
        "sampling_note": (
            "Occurrence records, not individuals. Monthly peaks measure reporting "
            "effort and phenology together; they cannot identify flowering or sex. "
            "The first N records are a nonrandom screening queue. Repeated photos "
            "of one plant, cultivated records, and misidentifications may remain."
        ),
        "filters_photo": FILTERS
    }
    for name in SPECIES:
        species = one_species(name, args.sample_size)
        result["species"].append(species)
        short = {
            "name": name, "gbif_key": species["taxon"]["gbif_taxon_key"],
            "counts": species["counts"],
            "months": species["monthly_counts_in_photo_subset"],
            "month_facet_valid": species["monthly_sum_matches_count"],
        }
        print("GBIF_AUDIT_SPECIES=" + json.dumps(short, ensure_ascii=False, sort_keys=True), flush=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
                      + "\n", encoding="utf-8")
    print("GBIF_AUDIT_OUTPUT=" + str(output), flush=True)

if __name__ == "__main__":
    main()
