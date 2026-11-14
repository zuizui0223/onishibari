#!/usr/bin/env python3
"""Effort-stratified GBIF overlap robustness checks for 3 Daphne taxa.

This is NOT a species distribution model, ecological niche estimate, or
proof of sympatry/speciation. Coordinates never leave this Python process.
It audits effects of sample size, basis of record, epoch and reported
coordinate uncertainty on 0.25-degree *observed* grid co-occurrence.
"""
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import argparse
import json
import math
from pathlib import Path
import random
from statistics import mean

from gbif_daphne_spatial_overlap import (
    SPECIES, download_taxon, valid, cell, GRID_DEGREES,
)

SEED = 20261114
REPLICATES = 400
CAP_RECORDS = 100
MIN_SUBSET = 10

def year_num(record):
    value = record.get("year")
    try:
        year = int(value)
    except (TypeError, ValueError):
        return None
    return year if 1800 <= year <= 2026 else None

def strata_predicates():
    return {
        "all_qc": lambda x: True,
        "preserved_only": lambda x: x.get("basisOfRecord") == "PRESERVED_SPECIMEN",
        "human_only": lambda x: x.get("basisOfRecord") == "HUMAN_OBSERVATION",
        "year_2000_onwards": lambda x: (year_num(x) or 0) >= 2000,
        "year_2016_onwards": lambda x: (year_num(x) or 0) >= 2016,
        "uncertainty_known_le_10km": lambda x:
            x.get("coordinateUncertaintyInMeters") is not None,
    }

def ecells(rows, spacing=GRID_DEGREES):
    return {cell(row, spacing) for row in rows}

def overlap(a, b):
    union = a | b
    return len(a & b) / len(union) if union else None

def pct(values, q):
    ordered = sorted(values)
    f = (len(ordered) - 1)*q
    i = math.floor(f)
    j = math.ceil(f)
    return (ordered[i]*(j-f) + ordered[j]*(f-i)) if i != j else ordered[i]

def calculate(rows_by_taxon, reps=REPLICATES, seed=SEED, cap=CAP_RECORDS):
    names=list(SPECIES)
    subset_sizes={k:len(rows_by_taxon[k]) for k in names}
    target=min(cap, min(subset_sizes.values()))
    all_grid={k:ecells(rows_by_taxon[k]) for k in names}
    all_pairs={}
    keys=[(a,b) for i,a in enumerate(names) for b in names[i+1:]]
    for a,b in keys:
        shared=len(all_grid[a]&all_grid[b])
        union=len(all_grid[a]|all_grid[b])
        all_pairs[f"{a} / {b}"]={
            "n_shared_cells":shared,"n_union_cells":union,
            "jaccard":round(shared/union,5) if union else None
        }
    result={
        "record_counts":subset_sizes,
        "occupied_cells":{k:len(v) for k,v in all_grid.items()},
        "all_records_pairwise":all_pairs,
        "balanced_target_records_per_species": target,
        "balance_draws":reps,
        "record_level_rarefaction":{},
    }
    if target < MIN_SUBSET:
        result["insufficient_for_rarefaction"]=True
        return result
    rng=random.Random(seed)
    draws={f"{a} / {b}":[] for a,b in keys}
    for iteration in range(reps):
        cells={
            k:ecells(rng.sample(rows_by_taxon[k], target)) for k in names
        }
        for a,b in keys:
            value=overlap(cells[a],cells[b])
            if value is not None:
                draws[f"{a} / {b}"].append(value)
    for key,vals in draws.items():
        if not vals:
            continue
        result["record_level_rarefaction"][key]={
            "median":round(pct(vals,.5),5),
            "p05":round(pct(vals,.05),5),
            "p95":round(pct(vals,.95),5),
            "mean":round(mean(vals),5),
        }
    return result

def get_data():
    rows={}
    fetch_meta={}
    for taxon, key in SPECIES.items():
        fetched, advertised=download_taxon(key)
        accepted=[]
        rejection=Counter()
        for record in fetched:
            reason=valid(record,key)
            if reason:
                rejection[reason]+=1
            else:
                accepted.append(record)
        rows[taxon]=accepted
        # Hash only IDs for deterministic provenance, never export identifiers.
        record_ids=sorted(str(x.get("key")) for x in accepted)
        ids_digest=sha256("\n".join(record_ids).encode()).hexdigest()
        fetch_meta[taxon]={
            "queried_count":advertised, "retained_qc":len(accepted),
            "excluded_by_reason":dict(rejection),
            "gbif_key":key, "record_ids_sha256":ids_digest,
        }
    return rows,fetch_meta

def self_test():
    # The core algorithm must not touch GBIF during a local check.
    rng=random.Random(123)
    def row(i,lat,lon,basis,year):
        return {"key":i, "decimalLatitude":lat,
                "decimalLongitude":lon,"basisOfRecord":basis, "year":year}
    obs={}
    for j,name in enumerate(SPECIES):
        obs[name]=[row(j*50+i,35+j+0.03*i,135+0.02*i,
                       "PRESERVED_SPECIMEN",2020) for i in range(20)]
    a=calculate(obs,reps=120,seed=9,cap=15)
    b=calculate(obs,reps=120,seed=9,cap=15)
    assert a==b
    assert a["balanced_target_records_per_species"]==15
    assert len(a["record_level_rarefaction"])==3
    assert abs(overlap({(1,1),(1,2)},{(1,2),(1,3)})-1/3)<1e-9
    assert year_num({"year":None}) is None and year_num({"year":"2020"})==2020
    assert year_num({"year":"unknown"}) is None
    print("EFFORT_SENSITIVITY_SELF_TEST=SUCCESS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="data/derived/gbif_daphne_effort_sensitivity.json")
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--reps",type=int,default=REPLICATES)
    args=p.parse_args()
    if args.self_test:
        self_test()
        return
    if not 100 <= args.reps <= 3000:
        p.error("--reps must be within 100..3000")
    clean,provenance=get_data()
    out={}
    for i,(name,condition) in enumerate(strata_predicates().items()):
        subset={sp:[r for r in clean[sp] if condition(r)] for sp in SPECIES}
        result=calculate(subset,reps=args.reps,seed=SEED+1009*i)
        out[name]=result
        print("EFFORT_STRATUM="+json.dumps({
            "stratum":name,
            "record_counts":result["record_counts"],
            "jaccard":{k:v["jaccard"] for k,v in result["all_records_pairwise"].items()},
            "balanced_n":result["balanced_target_records_per_species"],
            "balanced":result["record_level_rarefaction"],
        },sort_keys=True),flush=True)
    data={
        "retrieved_utc":datetime.now(timezone.utc).isoformat(),
        "source":"GBIF occurrence search API v1",
        "method":"quality-filtered 0.25-degree observed occurrence grid overlap",
        "coordinate_policy":"do not export points, locations, grids or GBIF record keys",
        "grid_degrees":GRID_DEGREES,
        "seed":SEED,"reps":args.reps,"record_cap":CAP_RECORDS,
        "source_provenance":provenance,
        "strata":out,
        "limitations":[
            "Jaccard of observed cells is NOT ecological range overlap or a niche model.",
            "Equal-record rarefaction controls unequal row totals, not unequal geographic sampling probabilities.",
            "Record-based draws can repeatedly represent the same plant or point; they are not independent individuals.",
            "Specimens and human observations have different spatial and temporal detectability.",
            "Coordinate uncertainty is unreported for many retained records; known-only subset is selected.",
            "Changing records by year changes source effort, so apparent temporal shifts are NOT range change.",
            "Taxonomy is imperfect; 2025 D. pseudomezereum var. koreana reclassification must be considered.",
        ],
    }
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("EFFORT_SENSITIVITY_OUTPUT="+str(path),flush=True)

if __name__=="__main__":
    main()
