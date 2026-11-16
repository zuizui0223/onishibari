#!/usr/bin/env python3
"""v0.32: adversarial controls on original Daphne paternal flower display signal.

Controls planned before inspecting these fits:
1. Maternal-family spatial block holdout (10 & 20m grids, 0 and half-cell shifts).
2. Restrict to mothers coded exactly F in original adult table.
3. Exclude the single largest H flowering genotype or all top 10% flowering genets,
   condition on remaining confidently assigned H offspring and recompute beta.
Only 10m exponential distance decay for primary sensitivity; original 5/10/20m
source has already been run in v30.

Outputs only aggregated coefficients and heldout log scores, no coordinates,
source IDs or individual rows.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import (
   build,estimate,boot,RNG
)
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data

OUT=Path("analysis_results_daphne_fatherhood_robust_v32")
LAMBDA=10.

def v(x):
    return round(float(x),6) if np.isfinite(x) else None

def reduced(q, mother_mask=None, sire_keep=None):
    counts=q["counts"];dist=q["distance"];lf=q["logflowers"]
    if mother_mask is None:mother_mask=np.ones(len(counts),dtype=bool)
    if sire_keep is None:sire_keep=np.ones(counts.shape[1],dtype=bool)
    c=counts[np.asarray(mother_mask)][:,np.asarray(sire_keep)]
    d=dist[np.asarray(mother_mask)][:,np.asarray(sire_keep)]
    notzero=c.sum(axis=1)>0
    if not np.any(notzero):raise RuntimeError("no eligible father assignments remain after source sensitivity")
    return {"counts":c[notzero],"distance":d[notzero],
       "logflowers":lf[np.asarray(sire_keep)]}

def block_cv(q,block_m=20,offset=0):
    xy=q["private_mother_xy"]
    assert xy.shape[0]==q["counts"].shape[0]
    groups=np.floor((xy+offset)/block_m).astype(int)
    _,lab=np.unique(groups,axis=0,return_inverse=True)
    nblock=len(np.unique(lab))
    count=q["counts"]
    dist=q["distance"]
    lf=q["logflowers"]
    pred=[]
    all_exposure=0
    for b in np.unique(lab):
        keep=lab!=b
        hold=~keep
        if keep.sum()<5:continue
        wt=keep.astype(float)
        fit=estimate(q,LAMBDA,weights=wt)
        coef=float(fit["beta_hat"])
        base=-dist[hold]/LAMBDA
        U=base+coef*lf[None,:]
        U0=base
        diff=np.sum(count[hold]*(U-logsumexp(U,axis=1)[:,None]-
                                  U0+logsumexp(U0,axis=1)[:,None]))
        n=float(count[hold].sum())
        pred.append({"gain":float(diff),"n":int(n)})
        all_exposure+=n
    if not pred:return {"n_blocks":nblock,"status":"insufficient_block_coverage"}
    return {
      "n_blocks":nblock,
      "n_test_blocks":len(pred),
      "n_test_assigned_offspring":int(all_exposure),
      "spatial_block_holdout_log_score_gain_per_child":v(sum(p["gain"] for p in pred)/all_exposure),
      "n_blocks_positive_holdout_gain":sum(p["gain"]>0 for p in pred),
      "n_blocks_negative_holdout_gain":sum(p["gain"]<0 for p in pred),
      "blocks_not_independent_if_shifted_or_adjoining":True}

def scenario(name,q):
    r=estimate(q,LAMBDA,extra=False)
    b=boot(q,LAMBDA)
    return {"scenario":name,
      "remaining_sampled_maternal_families":len(q["counts"]),
      "remaining_source_H_genets":len(q["logflowers"]),
      "remaining_confident_H_paternity_assignments":int(q["counts"].sum()),
      "beta":r["beta_hat"],
      "beta_mother_boot_CI95":b["beta_mother_boot_ci95"],
      "per_child_gain_distance_only":r["per_child_in_sample_gain_over_distance_only"]}

def main():
    OUT.mkdir(exist_ok=True)
    files=data()
    outs=[]
    for year in (2018,2020):
        q=build(year,files)
        baseline=estimate(q,LAMBDA,extra=True)
        exactF=(np.asarray(q["private_mother_sex_codes"])=="F")
        assert exactF.sum()==(35 if year==2018 else 32)
        LF=q["logflowers"]
        largest=np.argmax(LF)
        keep_single=np.ones(len(LF),dtype=bool);keep_single[largest]=False
        top_n=int(np.ceil(len(LF)*.1))
        keep_top=~np.isin(np.arange(len(LF)),np.argsort(LF)[-top_n:])
        scenarios=[
          scenario("all_eligible_H_fathers",reduced(q)),
          scenario("only_exact_F_mothers",reduced(q,mother_mask=exactF)),
          scenario("exclude_one_max_flower_H_genet",reduced(q,sire_keep=keep_single)),
          scenario("exclude_top_10pct_flower_H_genets",reduced(q,sire_keep=keep_top))]
        blocks=[]
        for b in (10,20):
            for shift in (0,b/2):
                blocks.append({"block_cell_size_m":b,"shift_both_coords_m":shift,
                               **block_cv(q,b,shift)})
        outs.append({
            "year":year,
            "source_matched_H_children":int(q["counts"].sum()),
            "source_H_candidate_genets":int(len(LF)),
            "source_mothers":int(len(q["counts"])),
            "mothers_exact_F_code":int(exactF.sum()),
            "baseline_beta":baseline["beta_hat"],
            "baseline_mother_holdout_gain_per_child":baseline["mother_leave_one_out_gain_per_child"],
            "subset_sensitivity":scenarios,
            "spatial_block_holdout":blocks
            })
    out={"status":"v32_original_Daphne_male_flower_count_paternity_source_robustness",
      "original_paper":"10.1002/ajb2.16239",
      "original_data":"10.5281/zenodo.8253921",
      "lambda_metres_preselected":10,
      "years":outs,
      "limits":{
        "subsetting_candidates_changes_conditioned_risk_set":True,
        "spatial_blocks_do_not_remove_environmental_confounding":True,
        "many_nearby_mothers_share_pollinator_opportunities":True,
        "bootstrap_over_mothers_not_randomized_flower_size":True,
        "source_original_year_specific_natural_flower_numbers_not_daily_pollen":True,
        "not_full_population_paternal_seed_production":True,
        "not_focal_D_pseudomezereum":True}}
    (OUT/"daphne_H_flower_paternity_adversarial_v32.json").write_text(
       json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_MALE_PATERNITY_SENSITIVITY_V32",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
