#!/usr/bin/env python3
"""Check Daphne phylogenetic/functional comparison without hallucinated traits.

Inputs are sourced literature assertions (NOT Onishibari raw biological data).
Do not make synonyms independent tips or conflate selfing, self compatibility,
autonomous pollen transfer, and female/androdioecious sex morphs.
"""
import csv
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]
EVIDENCE=ROOT/"data/derived/daphne_sex_selfing_phylogeny_evidence_v1.csv"
DIMS=ROOT/"data/derived/daphne_sex_morph_effect_size_source_means.csv"
UNKNOWN={"UNKNOWN","NOT EXPERIMENTALLY ESTABLISHED","NOT TESTED"}
SCHEMA={"taxon","evolutionary_comparison_tier","population",
        "floral_sex_system_reported","male_vs_female_floral_difference",
        "self_compatibility","autonomous_selfing","genetic_selfing_fraction",
        "estimate_method","key_evidence","primary_source","second_source",
        "uncertainty_or_conflict"}

def read(path):
    with path.open(encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle)
        if reader.fieldnames is None or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError("Missing or duplicated header")
        return set(reader.fieldnames),list(reader)

def check():
    names,rows=read(EVIDENCE)
    if names!=SCHEMA or len(rows)!=13:
        raise ValueError("Comparative evidence changed without manual re-audit")
    for row in rows:
        if not row["primary_source"].startswith("https://"):
            raise ValueError("Evidence lacks a verifiable primary source")
        if not row["uncertainty_or_conflict"]:
            raise ValueError("Each taxon must have uncertainty/limitation")
    keys=[(r["taxon"],r["population"]) for r in rows]
    if len(set(keys))!=len(keys):
        raise ValueError("Duplicate species/population strata")
    by={r["taxon"]:[] for r in rows}
    for row in rows:
        by[row["taxon"]].append(row)
    assert len(by)==12
    p=by["Daphne pseudomezereum"][0]
    assert p["self_compatibility"] in UNKNOWN
    assert p["genetic_selfing_fraction"]==""
    j=by["Daphne jezoensis"][0]
    assert "SELF COMPATIBLE" in j["self_compatibility"]
    assert "vector" in j["autonomous_selfing"]
    assert "1989" in j["uncertainty_or_conflict"]
    assert by["Daphne kamtschatica"][0]["self_compatibility"]=="NOT TESTED"
    syn=by["Daphne rebunensis"][0]
    assert syn["evolutionary_comparison_tier"].startswith("HISTORICAL SYNONYM")
    k=by["Daphne kamtschatica var. koreana"][0]
    assert "plastome" in k["key_evidence"]
    d=by["Daphne kiusiana"]
    assert len(d)==2
    east=next(r for r in d if r["population"].startswith("east"))
    west=next(r for r in d if r["population"].startswith("west"))
    assert (float(east["genetic_selfing_fraction"]),
            float(west["genetic_selfing_fraction"]))==(.920,.178)
    assert "Japanese field guides dioecy" in east["floral_sex_system_reported"]
    arb=by["Daphne arbuscula"][0]
    assert arb["genetic_selfing_fraction"]=="" # Zero selfed FRUITS is not population genetic selfing-rate estimate
    assert "168 hand-self" in arb["estimate_method"]
    assert "not demonstrated close sister" in arb["uncertainty_or_conflict"]
    mez=by["Daphne mezereum"][0]
    assert "UNRESOLVED" in mez["self_compatibility"]
    # A statement is not a new direct experimental SI assay in the target.
    assert not any("SELF INCOMPATIBLE" in r["self_compatibility"]
                   for r in by["Daphne pseudomezereum"])
    dims_cols,dims=read(DIMS)
    must={"taxon","study","trait","unit","female_mean","hermaphrodite_mean",
          "female_error","hermaphrodite_error","error_type",
          "female_n_flowers","female_n_plants","hermaphrodite_n_flowers",
          "hermaphrodite_n_plants","source_url","scope_note"}
    if dims_cols!=must or len(dims)!=5:
        raise ValueError("Five sex-dimorphism trait comparisons required")
    for row in dims:
        f=float(row["female_mean"])
        h=float(row["hermaphrodite_mean"])
        if not 0<f<h:
            raise ValueError("Female-vs-H floral contrast incorrect")
        if row["error_type"] not in {"SE","SD"}:
            raise ValueError("SD versus SE not identified")
        if row["taxon"] not in ("Daphne jezoensis","Daphne laureola"):
            raise ValueError("No numeric target sex data; do not claim Dpm measurements")
    assert sum(r["taxon"]=="Daphne laureola" for r in dims)==2
    assert len(set((r["taxon"],r["study"],r["trait"]) for r in dims))==5
    print("DAPHNE_COMPARISON_VALID=13 sourced rows; 12 taxon names (one synonym);"
          " 2 D. kiusiana populations; 5 non-target floral sex contrasts."
          " Dpm physiological SI UNTESTED. D. arbuscula no self fruit, genetic selfing fraction unknown, showy flowers.",
          flush=True)

if __name__=="__main__":
    check()
