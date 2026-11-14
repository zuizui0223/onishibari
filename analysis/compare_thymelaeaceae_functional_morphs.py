#!/usr/bin/env python3
"""Validate family-wide sexual-morph comparison and create auditable descriptors.

No phylogenetic-tree reconstruction or meta-analysis of overlapping publications.
Validate taxa and avoid treating nonfruiting functional males as physiological
self-incompatible hermaphrodites.
"""
import csv
import math
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"data/derived/thymelaeaceae_reproductive_function_cross_genera_v1.csv"
OUT=ROOT/"docs/2026-11-14_family_function_comparison_generated.md"
FIELDS=("reported_taxon","2026_taxon_context","comparison_tier","geographic_scope",
        "sex_system","morphological_evidence","female_function_evidence",
        "self_compatibility_or_SI","autonomous_selfing","genetic_selfing",
        "quantitative_trait","value","unit","n_structure","evidence_method",
        "primary_source","evidence_limit")
REQUIRED_TAXA={
    "Daphne pseudomezereum",
    "Daphne jezoensis",
    "Daphne kiusiana eastern lineage",
    "Daphne laureola",
    "Daphne arbuscula",
    "Daphnopsis filipedunculata",
    "Hawaiian Wikstroemia spp.",
    "Thymelina wikstroemiana",
    "Pimelea prostrata",
    "Pimelea traversii",
}
def load():
    with INPUT.open(encoding="utf-8",newline="") as stream:
        reader=csv.DictReader(stream)
        if tuple(reader.fieldnames or ())!=FIELDS:
            raise ValueError("Input schema drift")
        rows=list(reader)
    if len(rows)!=10 or {r["reported_taxon"] for r in rows}!=REQUIRED_TAXA:
        raise ValueError("Reviewed evidence set drift; requires human re-audit")
    if len(set(r["reported_taxon"] for r in rows))!=len(rows):
        raise ValueError("Duplicate biological comparison units")
    for r in rows:
        if not r["primary_source"].startswith("https://") or not r["evidence_limit"]:
            raise ValueError("Missing evidence source or limitation")
        if r["reported_taxon"]=="Daphnopsis filipedunculata":
            if r["comparison_tier"]!="other_family_branch":
                raise ValueError("Do not present this genus as a Daphne sister")
            if "NOT APPLICABLE" not in r["self_compatibility_or_SI"]:
                raise ValueError("Do not call functional male flower self-incompatible")
        if r["reported_taxon"]=="Daphne pseudomezereum":
            if r["self_compatibility_or_SI"]!="NOT_TESTED":
                raise ValueError("Target physiological SI not verified")
        if r["reported_taxon"]=="Thymelina wikstroemiana":
            if "2026" not in r["2026_taxon_context"]:
                raise ValueError("New Gnidia combination missing")
        if r["reported_taxon"]=="Hawaiian Wikstroemia spp.":
            if "unspecified" not in r["2026_taxon_context"]:
                raise ValueError("Do not imply species-level independent replication")
    return rows

def ratios(rows):
    pm={r["reported_taxon"]:r for r in rows}
    for name,value in [
        ("Pimelea prostrata",41.92/4.27),
        ("Pimelea traversii",83.54/40.64),
    ]:
        got=float(pm[name]["value"])
        if abs(got-value)>0.0015: raise ValueError("Fruit set ratio changed")
    assert pm["Daphnopsis filipedunculata"]["n_structure"].startswith("49 ")
    assert pm["Daphne arbuscula"]["genetic_selfing"]=="genetic selfing not estimated"
    assert pm["Daphne laureola"]["autonomous_selfing"]=="pollinator service required for successful selfing"
    return [(name,float(pm[name]["value"])) for name in ("Pimelea traversii","Pimelea prostrata")]

def generate(rows):
    pm={r["reported_taxon"]:r for r in rows}
    fs=ratios(rows)
    out=[
        "# Family-wide morph versus reproductive-function case inventory",
        "",
        "**Ten curated comparisons, not ten independent phylogenetic contrasts or a fitted model.**",
        "No new experiments on *D. pseudomezereum*, no asserted self-incompatibility for that species, no tree or branch lengths estimated.",
        "",
        "## The strongest direct functional comparison",
        "",
        "- *Daphnopsis filipedunculata*: male flowers retain a prominent pistillode, yet none of the multiple self/cross/open pollination treatments produced fruit; female open-pollination fruit percentage was 23.1%. A SNP sample identified 26 unique genets out of 49 visually sampled ramets. The mechanism is **cryptic functional dioecy**, NOT physiological self-incompatibility in a sexually functional hermaphrodite.",
        "- Hawaiian *Wikstroemia* complex: cryptic unisexual versus morphologically unisexual systems coexist across populations; crosses between different sex-type systems yielded many apparently perfect-flowered offspring. This is evidence of distinct genetic sex-expression architectures, not a quantitative selfing rate.",
        "- *Thymelina wikstroemiana* (formerly *Gnidia wikstroemiana*): female flowers smaller, styles longer, staminodes reduced; females in two populations were 36–37%, female open-fruit success 35%; females' total reproductive energy was 7.3% higher than hermaphrodites. Some hermaphrodites produced fruit after selfing, but no population genetic selfing rate was extracted.",
        "",
        "## Cross-population natural fruit success ratios for *Pimelea* (1960 historical results)",
        "",
        "| Species | Female fruit-set mean (%) | Hermaphrodite fruit-set mean (%) | Female / hermaphrodite |",
        "|---|---:|---:|---:|",
        f"| *P. traversii* | 83.54 | 40.64 | {fs[0][1]:.2f} × |",
        f"| *P. prostrata* | 41.92 | 4.27 | {fs[1][1]:.2f} × |",
        "",
        "**These are per-plant reported natural fruit-set means from Burrows (1960), not equal-pollen experimental estimates.** The original source reports some hermaphrodite self fruiting; female performance may reflect pollen quality, sex allocation, plant age and sampling. They are historical observations, not independent evolutionary origins.",
        "",
        "## Why similar small flowers do not imply selfing",
        "",
        "1. Female-specific smaller flowers in gynodioecy or cryptic dioecy reflect sexual allocation; females cannot autonomously self-fertilize because they produce no viable pollen.",
        "2. Hermaphrodite-specific reduced herkogamy, selfing genotype evidence and smaller flowers jointly qualify *D. kiusiana* eastern lineage as a genuine selfing syndrome.",
        "3. Outcross-self compatibility does not guarantee autonomous pollen deposition: *D. laureola* can produce mainly selfed seed genetically while requiring flower visitors.",
        "4. Showy, fragrant *D. arbuscula* failed self-fruit tests and was reported as self-incompatible, but is NOT a selfing-syndrome-morphology × SI counterexample.",
        "5. The pivotal target question is not simply 'self-incompatible?' but whether male-like pistils in *D. pseudomezereum* have any viable female function after matched outcross positive controls and contamination/pollen-quality checks.",
        "",
        "## Phylogenetic inference limits",
        "",
        "- An ITS study sampled 32 *Daphne* species and 238 samples, finding the genus broadly polyphyletic; some named *Daphne* taxa nested in *Wikstroemia* or *Edgeworthia*. Within *Daphne sensu stricto* the subgenus *Pseudomezereum* was monophyletic in that analysis but deep relationships remained unresolved. Do not apply a single fixed genus-level tree to all rows.",
        "- A 2026 family-level revision found *Gnidia* extensively polyphyletic and proposed *Thymelina wikstroemiana* as a new combination. The older 2006 ecology paper should be cited under its original name as well.",
        "- *Daphnopsis*, *Pimelea*, *Thymelina*, and Hawaiian *Wikstroemia* are independent-family **functional analogues** but not confirmed sister clades. Genetic distance, divergence time, and number of evolutionary transitions were not analyzed.",
        "- A novel phylogenetic meta-analysis will need independently rooted multilocus species trees, per-population functional sex data, synonym resolution, and effect-size sampling variance. The available table alone is not sufficient.",
        "",
        "## Source inventory",
    ]
    for row in rows:
        out.append(f"- *{row['reported_taxon']}*: {row['primary_source']} — {row['evidence_limit']}")
    out.append("")
    return "\n".join(out)

def check():
    rows=load()
    ratios(rows)
    assert "cryptic functional dioecy" in generate(rows)
    assert "10 independent" in generate(rows)
    print("FAMILY_REPRODUCTIVE_COMPARE_CHECK=SUCCESS 10 sourced cases; 4 genera; "
          "Pimelea fruit ratios validated; Dpm SI unknown; Daphnopsis functional male != SI",flush=True)

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    check()
    if not args.check:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(generate(load()),encoding="utf-8")
        print("FAMILY_REPRODUCTIVE_COMPARE_OUTPUT="+str(args.output))
