#!/usr/bin/env python3
"""Validate human photo triage against source iNaturalist observation metadata.

Outputs observation-level image phenophases and conservative flower morphology
opportunity: no mating-system, self-compatibility or sex-ratio inference.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

ALLOWED = {
    "open_flower_visible", "ripe_fruit_visible", "immature_fruit_visible",
    "bud_or_flower_uncertain", "unreviewable_license", "vegetative_only"
}

def validate(source_file, score_file):
    source = json.loads(Path(source_file).read_text(encoding="utf-8"))
    meta = {str(x["observation_id"]): x for x in source["observations"]}
    with open(score_file,encoding="utf-8",newline="") as stream:
        scores=list(csv.DictReader(stream))
    ids = [r["observation_id"] for r in scores]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate scored observation IDs")
    if set(ids) != set(meta):
        raise ValueError(f"Outdated manual scores: missing={sorted(set(meta)-set(ids))}, extra={sorted(set(ids)-set(meta))}")
    if source["total_observations"] != len(meta):
        raise ValueError("Inconsistent source record count")
    if source["total_photo_refs"] != sum(len(o["photos"]) for o in meta.values()):
        raise ValueError("Inconsistent photo-ref count")
    for row in scores:
        if row["image_state"] not in ALLOWED:
            raise ValueError("Invalid photo state " + row["image_state"])
        if row["sex_morph"] != "indeterminate":
            raise ValueError("Functional or presumed sex is not validated by imagery")
        if row["stigma_resolved"] != "unclear":
            raise ValueError("Stigma not verified from citizen-science contact sheets")
        if row["orange_central_structure"] not in ("anther_like_candidate", "not_scored_or_unclear"):
            raise ValueError("Unexpected anther morphology state")
        if row["review_status"] != "manual_v02":
            raise ValueError("Scoring protocol changed")
        if row["observation_url"] != meta[row["observation_id"]]["url"]:
            raise ValueError("Wrong source link " + row["observation_id"])
    return source, meta, scores

def summarise(source,meta,scores):
    cnt=Counter(r["image_state"] for r in scores)
    floral=[r for r in scores if r["image_state"]=="open_flower_visible"]
    ripe=[r for r in scores if r["image_state"]=="ripe_fruit_visible"]
    immature=[r for r in scores if r["image_state"]=="immature_fruit_visible"]
    orange=[r for r in scores if r["orange_central_structure"]=="anther_like_candidate"]
    def month_count(items):
        counts=Counter((meta[r["observation_id"]].get("observed_on") or "unknown")[5:7] for r in items)
        return ", ".join(f"{int(k)}月 {v}" for k,v in sorted(counts.items()))
    by_grade=Counter(o.get("quality_grade") or "unknown" for o in meta.values())
    flower_annotations = {
        oid for oid, o in meta.items()
        if any(a.get("term_id") == 12 and a.get("value_id") == 13
               for a in o.get("annotations", []))
    }
    fruit_annotations = {
        oid for oid, o in meta.items()
        if any(a.get("term_id") == 12 and a.get("value_id") == 14
               for a in o.get("annotations", []))
    }
    seen_flower = {r["observation_id"] for r in floral}
    seen_fruit = {r["observation_id"] for r in ripe+immature}
    flower_annotated_overlap = len(seen_flower & flower_annotations)
    fruit_annotated_overlap = len(seen_fruit & fruit_annotations)

    labels={
       "open_flower_visible":"Open flower discernible",
       "ripe_fruit_visible":"Ripe red fruit visible",
       "immature_fruit_visible":"Unripe green fruit visible",
       "bud_or_flower_uncertain":"Bud or flower opening unclear",
       "vegetative_only":"Vegetative only in supplied images",
       "unreviewable_license":"Unreviewed because image is not licensed/retrievable",
    }
    lines=[
       "# D. pseudomezereum: all-photo visual screening (v0.2)",
       "",
       "**This is an observation-level evidence audit of citizen-science photos, not an estimate of plant sex ratio, flowering prevalence, pollination, fecundity, or self-incompatibility.**",
       "",
       f"Source: iNaturalist Japan observations, 2016-01-01 to 2026-11-14. **{len(meta)} observation events, {source['total_photo_refs']} photographic references.** Photos are not independent plants.",
       "",
       "| Photo-observed state | Observation records |",
       "|---|---:|",
    ]
    for k in labels:
        lines.append(f"| {labels[k]} | {cnt[k]} |")
    lines+=["",
       f"- Research-grade observations: {by_grade['research']}; other grades: {len(meta)-by_grade['research']}. The image states were manually assessed across every licensed image of each observation, not the old first-photo-only queue.",
       f"- Dated open-flower evidence (n={len(floral)} observations): {month_count(floral)}.",
       f"- Ripe fruit evidence (n={len(ripe)} observations): {month_count(ripe)}.",
       f"- Unripe fruit evidence (n={len(immature)} observations): {month_count(immature)}.",
       f"- Yellow/orange anther-like central structures provisionally discernible in **{len(orange)}** observations; anatomical identity/pollen viability NOT verified. These are visual morphology candidates only.",
       "- Functional sexual morph is **indeterminate for every observation**. A male community annotation is not independent proof of viable pollen; fruit is evidence of female reproductive activity in that event, not necessarily a fixed lifetime morph.",
       "- An image with no visible open flower does not imply that the plant was not flowering elsewhere on the same shrub.",
       "- Unreviewed photos remain in the sample denominator but receive no reproductive stage score.",
       f"- iNaturalist flowering annotation (term 12/13) exists on {len(flower_annotations)} records, but only {flower_annotated_overlap} of {len(seen_flower)} visually confirmed open-flower records carry it; {len(seen_flower)-flower_annotated_overlap} visually confirmed flowering records lack this annotation.",
       f"- iNaturalist fruiting annotation (term 12/14) exists on {len(fruit_annotations)} records; {fruit_annotated_overlap} of {len(seen_fruit)} visually confirmed fruit records carry it; {len(seen_fruit)-fruit_annotated_overlap} visually confirmed fruit records lack this annotation. Inconclusive photographs are not false annotation proof.",
       "- iNaturalist and GBIF syndicate overlapping observations. Do not sum their record counts.",
       "- The recorded months are observation opportunities selected by photographers, NOT population flowering/fruiting curves; climate-date regressions require observation-effort and geographic controls.",
       "",
       "### Data provenance",
       "- Manual classification and explicit image-only uncertainty: data/derived/inat_visual_screen_v02.csv",
       "- Downloaded source metadata/photo IDs/licenses: data/derived/inat_all_photo_snapshot.json (dated artifact in the GitHub Actions run)",
       "- API retrieval: analysis/inat_all_photo_audit.py; gallery: analysis/inat_all_photo_gallery.py",
       "- Morphology scoring guidance: protocols/photo_sex_morph_scoring.md",
       "- Field proof of mating system still requires artificially crossed and selfed flowers plus pollen-tube/embryo and viable seed data.",
       "",
    ]
    return "\n".join(lines)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source", default="data/derived/inat_all_photo_snapshot.json")
    p.add_argument("--labels", default="data/derived/inat_visual_screen_v02.csv")
    p.add_argument("--output", default="docs/2026-11-14_all_photo_review.md")
    args=p.parse_args()
    src,meta,rows=validate(args.source,args.labels)
    txt=summarise(src,meta,rows)
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(txt,encoding="utf-8")
    print("ALL_PHOTO_MANUAL_CHECK=" + json.dumps({
        "observations":len(meta),
        "photos":src["total_photo_refs"],
        "image_states":dict(Counter(x["image_state"] for x in rows)),
        "sex_identified":0
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
