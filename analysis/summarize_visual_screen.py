#!/usr/bin/env python3
"""Check hand-screened first iNaturalist images; descriptive record counts only."""
import argparse
from collections import Counter
import csv
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "derived" / "inat_visual_screen_v01.csv"
VALID = {"fruit_visible","open_flower_visible","flower_or_bud_uncertain",
         "foliage_or_bud_only","not_reviewed_missing_license"}

def summarize(src=SRC):
    with src.open(newline="",encoding="utf-8") as fh:
        data=list(csv.DictReader(fh))
    if len(data)!=28 or len({r["observation_id"] for r in data})!=28:
        raise ValueError("Review-set size or unique ID count changed")
    if any(r["screen_class"] not in VALID for r in data):
        raise ValueError("Unknown screening category")
    if any(r["functional_sex"]!="unknown" for r in data):
        raise ValueError("Functional sex cannot be assigned here")
    if sum(r["user_annotation"]=="male" for r in data)!=1:
        raise ValueError("Source male annotation mismatch")
    if any(r["first_photo_only"]!="true" for r in data):
        raise ValueError("Unexpected screening scope change")
    for row in data:
        date.fromisoformat(row["observed_on"])
        if not row["source_url"].endswith("/"+row["observation_id"]):
            raise ValueError("Broken source URL")
    by_class=Counter(r["screen_class"] for r in data)
    flowers=[r for r in data if r["screen_class"]=="open_flower_visible"]
    fruits=[r for r in data if r["screen_class"]=="fruit_visible"]
    def months(rows):
        ctr=Counter(date.fromisoformat(r["observed_on"]).month for r in rows)
        return ", ".join(f"{m}: {ctr[m]}" for m in sorted(ctr))
    lines=[
        "# First-photo manual screening of Onishibari iNaturalist records",
        "",
        "**Screen v0.1, 2026-11-14. 28 research-grade observation events; one primary photograph per event. These are NOT population sex ratios, pollinator counts, or plant individual counts.**",
        "",
        "| Observed image state | Observation records |",
        "|---|---:|",
    ]
    labels={
        "open_flower_visible":"Open flower visibly present",
        "fruit_visible":"Fruit visibly present",
        "flower_or_bud_uncertain":"Open flower / bud uncertain",
        "foliage_or_bud_only":"Foliage/buds only in first photo",
        "not_reviewed_missing_license":"Not visually reviewed (missing/disallowed photo license)",
    }
    for klass in ("open_flower_visible","fruit_visible","flower_or_bud_uncertain",
                  "foliage_or_bud_only","not_reviewed_missing_license"):
        lines.append(f"| {labels[klass]} | {by_class[klass]} |")
    lines += [
        "",
        f"- Open-flower photo observation months (n={len(flowers)}): **{months(flowers)}**",
        f"- Fruit photo observation months (n={len(fruits)}): **{months(fruits)}**",
        "- One photo has a iNaturalist male community annotation: https://www.inaturalist.org/observations/339863389 ; zero have a female annotation. Not proof of functional sex.",
        "- **All 28 observations have functional sex unknown.** A fruit is evidence of seed-parent activity at some point but does not assign a permanent morph.",
        "- Other photos from the same observation can have useful evidence. Recheck all angles; calibrate with field-verified floral dissections.",
        "- Observation months reflect reporting effort and image availability; not a phenological population estimate.",
        "- The CSV has no coordinates or personally identifying data.",
        "",
        "### Reproducibility",
        "- Photo pages use iNaturalist source photos with license/attribution from analysis/build_inat_screening_sheet.py.",
        "- Manual classification: data/derived/inat_visual_screen_v01.csv",
        "- Public source: https://www.inaturalist.org/observations?taxon_id=447647&place_id=6737",
        "- GBIF + iNaturalist metadata audit: docs/2026-11-14_public_occurrence_audit.md",
        "",
    ]
    return "\n".join(lines)

if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--output",type=Path)
    args=cli.parse_args()
    report=summarize()
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(report,encoding="utf-8")
    else:
        print(report)
    print("PHOTO_SCREEN_VALIDATION=28 unique records, functional sex unknown")
