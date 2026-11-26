#!/usr/bin/env python3
"""Garden-paired bootstrap contrast of two nested but distinct reproductive endpoints.

Original Reji Chacko et al. 2016 Zurich gardens, published 2026.
Reuses previously verified source joins; garden is the bootstrap/permutation unit.
No source garden identifiers or coordinates are written to outputs. Findings are
descriptive associations: not causal pollen transfer, fertilization or biomass.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, rankdata
from zurich2016_flower_date_link import prepare, load, STUDY
from zurich2016_display_to_seed import CATS, BASE, csv_from_bytes, num

OUT = Path("analysis_results_zurich_stage_contrast")
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261115
REPS = 10000


def aggregate_endpoint(raw, spec):
    by_site = defaultdict(lambda: [0.0, 0.0])
    for r in raw:
        site = str(r.get("Id", "")).strip()
        if site in spec["exclude"]:
            continue
        if str(r.get("exclude", r.get("excluded", ""))).lower() in ("1", "yes", "true"):
            continue
        success = num(r.get(spec["success"]))
        fail = num(r.get(spec["fail"]))
        if success is None or fail is None or success < 0 or fail < 0 or success + fail <= 0:
            continue
        by_site[site][0] += success
        by_site[site][1] += success + fail
    return {k: n/d for k,(n,d) in by_site.items() if d > 0}


def rho(x, y):
    if len(np.unique(x)) < 3 or len(np.unique(y)) < 3:
        return float("nan")
    return float(spearmanr(x,y).statistic)


def partial_rank(x,y,covs):
    a = np.column_stack([rankdata(x),rankdata(y),*[rankdata(c) for c in covs]])
    design = np.column_stack([np.ones(len(x)),a[:,2:]])
    rx=a[:,0]-design@np.linalg.lstsq(design,a[:,0],rcond=None)[0]
    ry=a[:,1]-design@np.linalg.lstsq(design,a[:,1],rcond=None)[0]
    if np.std(rx)==0 or np.std(ry)==0:return float("nan")
    return float(np.corrcoef(rx,ry)[0,1])


def interval(xs):
    z=np.array([x for x in xs if np.isfinite(x)])
    return [round(float(v),4) for v in np.quantile(z,[.025,.5,.975])] if len(z)>100 else None


def paired_assessment(x,y1,y2,covs, rng):
    """Paired garden bootstrap: compare strength of associations on same sites."""
    raw1,raw2=rho(x,y1),rho(x,y2)
    adj1=partial_rank(x,y1,covs)
    adj2=partial_rank(x,y2,covs)
    n=len(x)
    boot=np.empty((REPS,4),float)
    for k in range(REPS):
        ii=rng.integers(0,n,n)
        boot[k]=[rho(x[ii],y1[ii]),rho(x[ii],y2[ii]),
            partial_rank(x[ii],y1[ii],[c[ii] for c in covs]),
            partial_rank(x[ii],y2[ii],[c[ii] for c in covs])]
    leave_one_out=np.array([[rho(np.delete(x,j),np.delete(y1,j)),
                             rho(np.delete(x,j),np.delete(y2,j))] for j in range(n)])
    permutation_count=np.zeros(2,int)
    for k in range(REPS):
        ix=rng.permutation(n)
        permutation_count[0] += abs(rho(x[ix],y1))>=abs(raw1)
        permutation_count[1] += abs(rho(x[ix],y2))>=abs(raw2)
    valid1=np.isfinite(boot[:,0]);valid2=np.isfinite(boot[:,1])
    valid_a1=np.isfinite(boot[:,2]);valid_a2=np.isfinite(boot[:,3])
    return {
        "number_of_complete_gardens":n,
        "same_gardens_for_both_endpoints":True,
        "first_seeded_flower_rho":round(raw1,4),
        "ovules_with_seed_rho":round(raw2,4),
        "paired_rho_difference_first_minus_ovules":round(raw1-raw2,4),
        "paired_rho_difference_bootstrap_95":interval(boot[valid1&valid2,0]-boot[valid1&valid2,1]),
        "flower_endpoint_bootstrap_95":interval(boot[:,0]),
        "ovule_endpoint_bootstrap_95":interval(boot[:,1]),
        "flower_endpoint_permutation_p_two_sided":round((permutation_count[0]+1)/(REPS+1),5),
        "ovule_endpoint_permutation_p_two_sided":round((permutation_count[1]+1)/(REPS+1),5),
        "flower_endpoint_loo_rho_minmax":[round(float(v),4) for v in [np.min(leave_one_out[:,0]),np.max(leave_one_out[:,0])]],
        "ovule_endpoint_loo_rho_minmax":[round(float(v),4) for v in [np.min(leave_one_out[:,1]),np.max(leave_one_out[:,1])]],
        "urban500_flower_richness_partial_first_rho":round(adj1,4),
        "urban500_flower_richness_partial_ovule_rho":round(adj2,4),
        "partial_rho_difference":round(adj1-adj2,4),
        "partial_rho_difference_bootstrap_95":interval(boot[valid_a1&valid_a2,2]-boot[valid_a1&valid_a2,3]),
        "bootstrap_repetitions":REPS,
        "permutation_repetitions":REPS,
        "seed":SEED,
        "interpretation_guardrail":"Difference nonzero only if paired-bootstrap percentile interval excludes zero; this remains observational and does not isolate the fertilization stage."
    }


def main():
    report={"status":"failed",
       "source":"Reji Chacko et al. (2026) J Applied Ecology DOI 10.1111/1365-2664.70384; EnviDat DOI 10.16904/envidat.676",
       "question":"Is the flower-with-any-seed association statistically stronger than the ovule-completion association on identical gardens?",
       "original_analysis_not_replication":"The original paper already reported pollinator guild vs pollination success. This adds a paired strength-of-association uncertainty test, NOT a new independent experiment.",
       "privacy":"No garden identifiers, coordinates or original private records exported",
       "limitations":[
          "Non-significant ovule endpoint does not establish a late-stage bottleneck.",
          "Each garden's visitors and season-end seed outcomes may have shared resource/environment confounding.",
          "The two outcomes are not measured on identical individual flowers and are not a per-ovule longitudinal transition.",
          "Recorded insect landings are not measured compatible conspecific pollen deposition.",
          "The flower-number exposure is an approximate unit-hours sum and can be a mixed inflorescence measure.",
          "Even if the correlation contrast is positive, that is not experimental evidence for a block at fertilization."
       ]}
    try:
        field,qa,source = prepare()
        effort = defaultdict(lambda:[0.0,0.0])
        for r in field:
            if r["plant"]!="comfrey":continue
            site=r["site"];hours=r["effort_minutes"]/60
            effort[site][0]+=hours
            effort[site][1]+=hours*r["flower_units"]
        b=load(STUDY,80_000_000)
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            f=aggregate_endpoint(csv_from_bytes(z.read(CATS["comfrey"]["file"])), CATS["comfrey"])
            o=aggregate_endpoint(csv_from_bytes(z.read(CATS["comfrey_seed_ovules"]["file"])),CATS["comfrey_seed_ovules"])
        ab_raw=load(BASE+"cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt",500000)
        urban_raw=load(BASE+"raw_data/explanatory_variables.txt",500000)
        flowers_raw=load(BASE+"raw_data/plant_floristic_data.txt",500000)
        abundance={str(r["Id"]):r for r in csv_from_bytes(ab_raw,delimiter=";")}
        urban={str(r["Id"]):num(r["Urban_500"]) for r in csv_from_bytes(urban_raw,delimiter=";")}
        richness={str(r["Id"]):num(r["SR_all_insect_pollinated_May_August"]) for r in csv_from_bytes(flowers_raw,delimiter=";")}
        rows=[]
        for site in sorted(set(f)&set(o)):
            if site not in abundance or site not in effort:continue
            h,fh=effort[site]
            visitors=num(abundance[site].get("A_Bombus_Comfrey"))
            u=urban.get(site);z=richness.get(site)
            if not all(v is not None for v in [visitors,u,z]) or h<=0 or fh<=0:continue
            rows.append((visitors/h,visitors/fh,fh/h,f[site],o[site],u,z))
        if len(rows)!=23:raise ValueError("Expected 23 complete comfrey gardens; got "+str(len(rows))+"; source changed or join failed")
        a=np.array(rows,float)
        if not np.all(np.isfinite(a)):raise ValueError("Nonfinite in garden-level input")
        covs=[a[:,5],a[:,6]]
        rng=np.random.default_rng(SEED)
        report["per_observation_hour"]=paired_assessment(a[:,0],a[:,3],a[:,4],covs,rng)
        report["per_approx_flower_unit_hour"]=paired_assessment(a[:,1],a[:,3],a[:,4],covs,rng)
        report["floral_display_proxy_rho_with_flower_endpoint"]=round(rho(a[:,2],a[:,3]),4)
        report["floral_display_proxy_rho_with_ovule_endpoint"]=round(rho(a[:,2],a[:,4]),4)
        report["original_site_source_count"]=qa["unique_sites"]
        report["source_sha256"]={
            "envidat_zip":hashlib.sha256(b).hexdigest(),
            "source_guild_abundance":hashlib.sha256(ab_raw).hexdigest(),
            "urban500":hashlib.sha256(urban_raw).hexdigest(),
            "flowering_richness":hashlib.sha256(flowers_raw).hexdigest()
        }
        # Shared denominator: both endpoints are on same gardens but are neither
        # identical flowers nor causal transitions. Avoid fitting arrows V->F->O.
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
    out=OUT/"paired_comfrey_stage_associations.json"
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False),flush=True)
    if report["status"]!="success":raise RuntimeError(report["error"])


if __name__=="__main__":
    main()
