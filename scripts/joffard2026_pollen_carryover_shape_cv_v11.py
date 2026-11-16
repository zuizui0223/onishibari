#!/usr/bin/env python3
"""A held-out-bee comparison of alternative original pollen-carryover shapes.

Primary public records: Joffard et al. (2026) AJB DOI 10.1002/ajb2.70276,
Dataverse DOI 10.57745/CF1VNW, pollen_counts_on_stigmas.tab.
Original data source SHA256 independently fixed in v0.9.
185 unique (bumblebee, visitation-position) cells from 186 source rows.
This is a conditional WITHIN-BEE shape analysis. It does NOT predict the
absolute pollen count of a new bee without measuring its pollen load. It does
NOT identify self-incompatibility, pollen germination, seed set or lifetime
fitness, and it does not join the separate female paternity experiment.
"""
from __future__ import annotations
import hashlib
import io
import json
import math
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

OUT=Path("analysis_results_pollen_carryover_shape_v11")
OUT.mkdir(exist_ok=True)
URL="https://entrepot.recherche.data.gouv.fr/api/access/datafile/765423"
EXPECTED_SHA="8291aea33a6c2894edccc55dd911f52a5c9c014049c20fb04ce7397f7e6dee68"
SEED=20261116
RNG=np.random.default_rng(SEED)
BOOT=7500
SHAPES={
    "exponential_position":lambda pos: np.column_stack([np.asarray(pos,float)]),
    "power_log_position":lambda pos: np.column_stack([np.log(np.asarray(pos,float))]),
    "early_step_plus_later_linear":lambda pos: np.column_stack(
        [(np.asarray(pos,float)>=2).astype(float),
         np.maximum(0,np.asarray(pos,float)-2)]),
    "continuous_kink_at_2":lambda pos: np.column_stack(
        [np.minimum(np.asarray(pos,float)-1,1),
         np.maximum(0,np.asarray(pos,float)-2)]),
    "continuous_kink_at_3":lambda pos: np.column_stack(
        [np.minimum(np.asarray(pos,float)-1,2),
         np.maximum(0,np.asarray(pos,float)-3)]),
}
# early_step+later_linear and continuous_kink_at_2 are exactly equivalent
# basis functions (first coefficient sets 1->2 log change). Keep only one
# candidate to prevent duplicate models being called separate evidence.
SHAPES.pop("early_step_plus_later_linear")

def scalar(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return round(float(x),7) if np.isfinite(x) else None
    if isinstance(x,(list,tuple,np.ndarray)):return [scalar(y) for y in x]
    if isinstance(x,dict):return {str(k):scalar(v) for k,v in x.items()}
    return x

def source():
    req=urllib.request.Request(URL,headers={"User-Agent":"Published-Silene-bee-holdout-validation-v11/0.1","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as r:b=r.read(100000)
    if hashlib.sha256(b).hexdigest()!=EXPECTED_SHA:raise ValueError("original public pollen table SHA mismatch")
    df=pd.read_csv(io.BytesIO(b),sep="\t",dtype=str)
    if list(df.columns)!=["pollinator_ID","position_of_the_flower_in_the_sequence",
                         "female_plant_ID","number_of_pollen_grains"]:
        raise ValueError("unexpected original data schema")
    df["pos"]=pd.to_numeric(df.position_of_the_flower_in_the_sequence,errors="coerce")
    df["grains"]=pd.to_numeric(df.number_of_pollen_grains,errors="coerce")
    if df[["pollinator_ID","pos","grains"]].isna().any().any():raise ValueError("missing original data")
    if not df.pos.between(1,10).all() or (df.grains<0).any():raise ValueError("invalid grains/position")
    if len(df)!=186 or df.pollinator_ID.nunique()!=25:raise ValueError("unexpected bee units")
    ndupe=int(df.duplicated(["pollinator_ID","pos"]).sum())
    if ndupe!=1:raise ValueError("unexpected duplicate bee-position grain")
    # Duplicate bee + ordinal position: use the average of two measured stigmas.
    # A unique bee-position is treated as the observation grain.
    a=df.groupby(["pollinator_ID","pos"],as_index=False).agg(grains=("grains","mean"))
    a["y"]=np.log1p(a.grains)
    if len(a)!=185:raise ValueError("unexpected bee-position units")
    return a,{"n_original_stigmas":186,"n_bees":25,"n_duplicate_bee_positions":ndupe,
              "n_unique_bee_position_cells":len(a),"original_data_sha256":EXPECTED_SHA}

def bee_demeaned(y,x,bee):
    """Within-bee intercept removal, preserving within-bee observation grain."""
    y=np.asarray(y,float)
    x=np.asarray(x,float)
    bee=np.asarray(bee)
    a=y.copy();z=x.copy()
    for b in np.unique(bee):
        mask=bee==b
        a[mask]-=a[mask].mean()
        z[mask,:]-=z[mask,:].mean(axis=0)
    return a,z

def fit_shape(frame,name):
    x=SHAPES[name](frame.pos.to_numpy(float))
    y=np.asarray(frame.y,float)
    bees=frame.pollinator_ID.to_numpy()
    yd,xd=bee_demeaned(y,x,bees)
    fit=sm.OLS(yd,xd).fit(cov_type="cluster",cov_kwds={"groups":bees})
    yh=xd@fit.params
    resid=yd-yh
    return {
        "fit":fit,
        "total_sse_centered":float(np.dot(resid,resid)),
        "coefficients":[float(v) for v in fit.params],
        "cluster_95":fit.conf_int().tolist(),
        "cluster_p":[float(v) for v in fit.pvalues],
        "source_n_bees":len(np.unique(bees)),
        "source_n_bee_positions":len(y),
        "empirical_within_bee_r2":float(1-np.sum(resid**2)/np.sum(yd**2)),
    }

def loo(frame,name):
    bybee=[]
    bees=sorted(frame.pollinator_ID.unique())
    for b in bees:
        tr=frame[frame.pollinator_ID!=b].copy()
        te=frame[frame.pollinator_ID==b].copy()
        tx=SHAPES[name](tr.pos.to_numpy(float));vy=tr.y.to_numpy(float)
        cy,cx=bee_demeaned(vy,tx,tr.pollinator_ID.to_numpy())
        beta=np.linalg.lstsq(cx,cy,rcond=None)[0]
        ux=SHAPES[name](te.pos.to_numpy(float))
        yy,xx=bee_demeaned(te.y.to_numpy(float),ux,np.repeat("heldout",len(te)))
        err=yy-xx@beta
        bybee.append({"bee_key_internal_only":b,"sse":float(np.dot(err,err)),
                      "mse":float(np.mean(err**2)),"n_positions":len(te),
                      "fold_predicted_with_test_bee_mean_subtracted":True})
    assert len(bybee)==25
    return bybee

def bootstrap_paired_rmse(a,b):
    """Compare two held-out bee models, paired on identical bee labels."""
    if len(a)!=len(b):raise ValueError("different validation groups")
    sa={x["bee_key_internal_only"]:x for x in a}
    sb={x["bee_key_internal_only"]:x for x in b}
    if set(sa)!=set(sb):raise ValueError("different held-out bees")
    labels=list(sa)
    aa=np.array([sa[x]["mse"] for x in labels])
    bb=np.array([sb[x]["mse"] for x in labels])
    diff=aa-bb
    ib=RNG.integers(0,len(diff),(BOOT,len(diff)))
    boot=diff[ib].mean(axis=1)
    return {"n_bees":len(labels),
        "mean_bee_mse_A_minus_B":scalar(diff.mean()),
        "paired_bee_bootstrap_95_mse_difference":[scalar(x) for x in np.quantile(boot,[.025,.975])],
        "n_bees_for_which_A_has_lower_mse":int((diff<0).sum()),
        "not_a_significance_test_of_a_prior_hypothesis":True}

def zeros_source(frame):
    original=frame.copy()
    original["zero"]=(original.grains==0).astype(int)
    # Correlated stigma zeros from same bee. Exchangeable working correlation.
    gee=smf.gee("zero ~ pos",groups="pollinator_ID",data=original,
        family=sm.families.Binomial(),cov_struct=sm.cov_struct.Exchangeable()).fit()
    b=float(gee.params["pos"]);ci=[float(x) for x in gee.conf_int().loc["pos"]]
    return {"n_bees":25,"n_bee_positions":len(original),
        "zero_grain_rate":scalar(original.zero.mean()),
        "gee_logit_or_per_additional_visit":scalar(math.exp(b)),
        "gee_robust_OR_95":[scalar(math.exp(x)) for x in ci],
        "p":scalar(gee.pvalues["pos"]),
        "gee_is_observation_model_not_absolute_mate_fitness":True}

def main():
    frame,qa=source()
    report={"status":"failed","study":"Joffard et al 2026 AJB doi:10.1002/ajb2.70276",
            "data_doi":"10.57745/CF1VNW","source_QA":qa,
            "models_compared":"within_bee_log1p_deposited_grains; observed ordinal positions 1..10",
            "cross_validation":"leave_one_bumblebee_out; evaluate demeaned within-bee shape, NOT prediction of unknown bee absolute pollen load",
            "bootstrap_reps":BOOT,"seed":SEED,
            "source_pollen_deposition_not_compatible_tubes_or_seeds":True,
            "held_out_bee_and_plant_ids_not_exported":True,
            "candidates":{},"pairwise_cv":{}}
    loores={}
    for name in SHAPES:
        res=fit_shape(frame,name)
        out=loo(frame,name);loores[name]=out
        report["candidates"][name]={
           "n_parameters":len(res["coefficients"]),
           "coefficient_per_basis":scalar(res["coefficients"]),
           "cluster95_by_basis":scalar(res["cluster_95"]),
           "cluster_p_by_basis":scalar(res["cluster_p"]),
           "within_bee_R2":scalar(res["empirical_within_bee_r2"]),
           "loo_bee_mean_MSE":scalar(np.mean([x["mse"] for x in out])),
           "loo_pooled_position_MSE":scalar(sum(x["sse"] for x in out)/sum(x["n_positions"] for x in out)),
           "n_bees":25,"n_unique_bee_positions":len(frame),
           "individual_bee_fold_performance_not_released":True}
    names=list(SHAPES)
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            a,b=names[i],names[j]
            report["pairwise_cv"][a+"_vs_"+b]=bootstrap_paired_rmse(loores[a],loores[b])
    chosen=min(names,key=lambda z:report["candidates"][z]["loo_bee_mean_MSE"])
    report["selected_by_exploratory_cv"]=chosen
    report["best_mse"]=report["candidates"][chosen]["loo_bee_mean_MSE"]
    report["ordinal_zero_grain_risk"]=zeros_source(frame)
    report["status"]="success_original_held_out_bee_shape_models"
    (OUT/"pollen_carryover_source_heldout_bee_shapes_v11.json").write_text(
        json.dumps(scalar(report),ensure_ascii=False,indent=2)+"\n")
    print("BEE_CARRYOVER_CV",json.dumps(scalar(report),ensure_ascii=False)[:30000],flush=True)

if __name__=="__main__":main()
