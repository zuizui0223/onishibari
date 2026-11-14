#!/usr/bin/env python3
"""Primary-source checks and historical-vs-modern photo phenology description.

This is *not* a test of changes in flowering phenology since 1955.
Hamaya (1955) section-wide sex-specific intervals must not be mislabelled
as quantified measurements for Daphne pseudomezereum alone.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
HIST=ROOT/"data/derived/primary_onishibari_morph_phenology_evidence.csv"
SCORES=ROOT/"data/derived/inat_visual_screen_v02.csv"
LIVE=ROOT/"data/derived/inat_all_photo_snapshot.json"
DOC="https://www.jstage.jst.go.jp/article/jjapbot/30/2/30_30_2_3828/_pdf/-char/en"

def read_csv(path):
    with open(path, encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file))

def validate_primary(rows):
    if len(rows)!=17:
        raise ValueError("Expected 17 bounded primary evidence entries")
    ids=[(r["source_id"],r["taxon_scope"],r["flower_morph"],r["variable"]) for r in rows]
    if len(set(ids))!=len(ids):
        raise ValueError("Duplicate primary evidence key")
    for r in rows:
        if r["source_id"] not in ("Hamaya1955","Lei2020"):
            raise ValueError("Unknown primary source")
        if r["source_id"]=="Hamaya1955" and r["source_url"]!=DOC:
            raise ValueError("Incorrect 1955 primary source")
        if r["source_id"]=="Lei2020" and r["source_url"]!="https://doi.org/10.1111/ppl.12972":
            raise ValueError("Incorrect ecophysiology DOI")
        if r["scope_type"]=="section_general" and "five taxa" not in r["taxon_scope"]:
            raise ValueError("Section-wide trait recoded as species-specific")
        if r["variable"] in ("calyx_tube_length","stamen_length"):
            if r["scope_type"]!="section_general":
                raise ValueError("Old sex-specific interval has no Dpm-only denominator")
            if not (r["min"] and r["max"] and float(r["min"])<=float(r["max"])):
                raise ValueError("Malformed range")
    known={ (r["taxon_scope"],r["variable"]):r for r in rows if r["scope_type"]=="taxon_specific" }
    for sp in ("Daphne pseudomezereum","Daphne jezoensis"):
        if (sp,"reported_flowering_month_window") not in known:
            raise ValueError("Missing taxon-specific flowering window")
    if known[("Daphne pseudomezereum","calyx_lobe_to_tube_ratio")]["min"]!="0.5":
        raise ValueError("Dpm calyx ratio was silently modified")
    if known[("Daphne jezoensis","calyx_lobe_to_tube_ratio")]["max"]!="1":
        raise ValueError("Djezo calyx ratio was silently modified")
    print("HAMAYA_PRIMARY_CHECK=17 entries; 4 section-wide sex ranges; historical species distinction retained")
    return known

def describe(known, snapshot, labels):
    obs={str(row["observation_id"]): row for row in snapshot["observations"]}
    scores={row["observation_id"]:row for row in labels}
    if len(obs)!=snapshot["total_observations"] or len(obs)!=len(scores):
        raise ValueError("Observation count mismatch")
    if set(obs)!=set(scores):
        raise ValueError("Review IDs differ from source IDs")
    if any(x["sex_morph"]!="indeterminate" for x in labels):
        raise ValueError("Unverified sex labels were introduced")
    if any(x["taxon_id_at_observation"] != 447647 for x in obs.values()):
        raise ValueError("Unexpected iNat taxon identification")
    flowers=[obs[k] for k,s in scores.items() if s["image_state"]=="open_flower_visible"]
    cnt=Counter(int(x["observed_on"][5:7]) for x in flowers if x.get("observed_on"))
    dpm=known[("Daphne pseudomezereum","reported_flowering_month_window")]
    first,last=int(dpm["min"]),int(dpm["max"])
    historical_months=set(list(range(first,13))+list(range(1,last+1))) if first>last else set(range(first,last+1))
    match=sum(n for m,n in cnt.items() if m in historical_months)
    grades=Counter(x["quality_grade"] for x in flowers)
    lines=[
        "# Historical taxonomic phenology vs photo-verified event months",
        "",
        "**DESCRIPTIVE ONLY, not a 1955–2026 time trend.** Flower observations are biased by species verification, photography choices, season, site and climate; Hamaya (1955) supplied a diagnostic window, not a measured temporal sample.",
        "",
        f"- Source: [Hamaya (1955)]({DOC}), original taxon key for *D. pseudomezereum*: **December–March**; *D. jezoensis*: **April–May**.",
        f"- Modern iNaturalist: **{len(flowers)} photo-confirmed open-flower observation events**, across {snapshot['total_observations']} photo-bearing observations. **{match}/{len(flowers)}** fall in the historical Dpm diagnostic month window.",
        f"- Visible-flower event months: "+", ".join(f"{m}月: {cnt[m]}件" for m in sorted(cnt))+".",
        f"- Quality of those photo-confirmed events: research grade {grades['research']}, needs ID {grades['needs_id']}. The latter require taxonomic re-review.",
        "",
        "## Scientific interpretation",
        "The overlap is compatible with the historical flowering-season description but does not demonstrate phenological stasis or species identity: a photograph can be misidentified, population sampling is uneven, and the 1955 key included taxonomic/calendar generalizations.",
        "Most importantly, phenology is not a proxy for the sex ratio, self-incompatibility or the genetic parentage of fruit.",
        "",
        "## Historical floral morphology",
        "- Section-general male calyx tubes 8–10 mm vs female 4–6 mm, without species-specific counts (Hamaya 1955).",
        "- Female ovaries relatively large and stigmas usually reach the height of lower stamens, but this does not establish sterility of any pistil.",
        "- Taxonomic cross-check: *D. pseudomezereum* calyx-lobe/tube approx 0.5, basal lateral leaf veins irregular; *D. jezoensis* ratio approx 1, veins more regular. These are fallible traits, not automatic labels.",
        "- Later molecular/SEM papers raise alternative taxon delimitations; Kew presently recognizes Dpm and Djezo as separate accepted species.",
        "",
        "## Seasonal allocation hypothesis",
        "Hamaya (1955) documented unusual overlap of fruits, old leaves and new inflorescences; Lei et al. (2020) experimentally connected summer leaf loss and autumn–winter carbon economy. That **does not establish sex-specific carbon costs**. Measure stored NSC and phenology separately in tagged male-biased versus female shrubs before hypothesizing causation.",
        "",
    ]
    return "\n".join(lines), dict(count=cnt, grades=grades, matched=match, total=len(flowers))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true",help="Only validate literature matrix; no network data")
    p.add_argument("--source",type=Path,default=LIVE)
    p.add_argument("--labels",type=Path,default=SCORES)
    p.add_argument("--output",type=Path,default=ROOT/"docs/hamaya_vs_photo_validation.md")
    args=p.parse_args()
    key=validate_primary(read_csv(HIST))
    if args.check:
        return
    if not args.source.exists():
        raise SystemExit("Source snapshot missing. Execute in iNaturalist all-photo audit workflow.")
    snap=json.loads(args.source.read_text(encoding="utf-8"))
    report,stats=describe(key,snap,read_csv(args.labels))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(report,encoding="utf-8")
    print("HISTORICAL_PHENOLOGY_AUDIT="+json.dumps(stats,default=dict,sort_keys=True))

if __name__=="__main__":
    main()
