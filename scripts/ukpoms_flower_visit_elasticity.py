#!/usr/bin/env python3
"""Floral resource scaling: total insects versus per-floral-unit landings.

Hypothesis: For a fixed 10min/50cm PoMS FIT quadrat, visits scale as
    E[visits | X] = exp(alpha + beta * ln(F)) * exp(X gamma)
with F the source "floral_unit_count". beta<1 implies SUBPROPORTIONAL
number of landings per floral unit; beta=1 proportional; beta>1 concentration.
A "floral unit" can be a flower, flower head, umbel, or spike. Do not call
Y/F an effective pollen delivery probability or plant lifetime fitness.

This is an observational, one-national-survey association, not adaptation.
Poisson QMLE with 1km-square cluster robust errors; sensitivity with
within-square FE, exclude unnamed "Other" target flower category, drop top
1% of F, and negative-binomial GLM. No 1km square strings exported.
"""
from __future__ import annotations
import json
import math
import numpy as np
from pathlib import Path
import pandas as pd
import patsy
import scipy.stats as ss
import statsmodels.api as sm

from ukpoms_flower_vs_landcover_cv import prepare,matrix

OUT=Path("analysis_results_flower_scaling")
OUT.mkdir(exist_ok=True)
GUILDS=["all_insects_total","all_bees","hoverflies","bumblebees"]
FORMULA_ENV=(
    "log_floral_units + C(target_flower_corrected) + C(floral_unit)"
    " + C(flower_structure) + C(flower_cover) + C(flower_context)"
    " + C(year) + C(month) + hour + C(sunshine)"
    " + C(wind_speed) + C(cloud_cover) + C(land_cover) + C(habitat_type)"
)
FORMULA_WITHIN=(
    "log_floral_units + C(target_flower_corrected) + C(floral_unit)"
    " + C(flower_structure) + C(flower_cover) + C(flower_context)"
    " + C(year) + C(month) + hour + C(sunshine)"
    " + C(wind_speed) + C(cloud_cover) + C(square_id)"
)
SUMMARY_COLUMNS=["study","rows","squares","flower_types","min_flower_units",
                 "median_flower_units","max_flower_units"]

def fit_one(subset,design,y,variant,model_type):
    if model_type=="poisson":
        family=sm.families.Poisson()
    elif model_type=="negbin_alpha1":
        family=sm.families.NegativeBinomial(alpha=1)
    else:raise ValueError(model_type)
    ydata=subset[y].to_numpy(float)
    X=patsy.dmatrix(design,subset,return_type="dataframe",NA_action="raise")
    if not np.isfinite(X.to_numpy()).all():raise ValueError("Nonfinite design")
    groups=subset["square_id"].astype(str).to_numpy()
    model=sm.GLM(ydata,X,family=family)
    result=model.fit(maxiter=75,tol=1e-8,method="irls",
         cov_type="cluster",cov_kwds={"groups":groups})
    beta=float(result.params["log_floral_units"])
    se=float(result.bse["log_floral_units"])
    if not math.isfinite(beta) or not (math.isfinite(se) and se>0):
        raise ValueError("Invalid coefficient/SE "+str((beta,se)))
    z1=(beta-1)/se
    z0=beta/se
    q99=float(np.quantile(subset["floral_unit_count"],.99))
    overdisp=float(result.pearson_chi2/result.df_resid)
    return {
      "response":y,"variant":variant,"link":model_type,
      "model":"log E[10min visits] = beta * log floral units + controls",
      "n_visits_sum":int(np.round(ydata.sum())),
      "n_zero_count_rows":int((ydata==0).sum()),
      "n_observation_rows":len(subset),
      "n_distinct_1km_squares":len(set(groups)),
      "n_parameters_design":X.shape[1],
      "log_units_beta":round(beta,6),
      "cluster_se":round(se,6),
      "cluster_95pct_ci":list(map(lambda t:round(t,6),(beta-1.96*se,beta+1.96*se))),
      "p_vs_beta0_more_total_visits":round(float(2*ss.norm.sf(abs(z0))),8),
      "p_vs_beta1_proportional_units":round(float(2*ss.norm.sf(abs(z1))),8),
      "doubling_flower_units_total_visits_multiplier":round(2**beta,6),
      "doubling_flower_units_visits_per_unit_multiplier":round(2**(beta-1),6),
      "pseudopoisson_pearson_overdispersion":round(overdisp,4),
      "source_flower_count_99th_quantile":round(q99,2),
      "notes":"The beta=1 test is about scaling of recorded landings across flower units, not actual per-flower pollen deposition",
    }

def rowsummary(data, variant):
    x=data["floral_unit_count"]
    return {"variant":variant,"rows":len(data),
        "squares":int(data["square_id"].nunique()),
        "flower_targets":int(data["target_flower_corrected"].nunique()),
        "family_masked":"No raw coordinates or square IDs exported",
        "flower_unit_5_50_95_99pct":list(
            map(float,np.quantile(x,[.05,.5,.95,.99]))),
        "n_other_target":int(data["target_flower_corrected"].str.contains("Other",case=False).sum()),
        "count_units":"flower head, individual flower, flower spike, umbel" }

def main():
    report={"study":"UK PoMS FIT Counts 2017–2023; EIDC doi:10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "status":"failed",
      "hypothesis":"Quantitative floral display yields sub/super-linear visitor count scaling at fixed observation effort",
      "formal_unit":"10min / 0.25m2 flower-insect quadrat (not reproductive fitness)",
      "statistical_model":"Poisson quasi-likelihood with 1km-square clustered errors, log link and site/target/weather controls",
      "cluster_variance_not_evolutionary_inference":True,
      "interpretation":"beta<1 is a landings-per-floral-unit dilution proxy, not per-flower pollen delivery",
      "limitations":[
       "Observational flower-selection: PoMS observer selects a flowering target in quadrat",
       "Unmeasured visitor efficiency, clonal floral unit structure, mate compatibility and seed set",
       "Species-specific flower units may be head, spike, umbel, individual flower; category controls cannot fully harmonize units",
       "Site-blocked SE and site-fixed effects reduce spatial confounding, not establish density causality",
       "Square fixed effects plus category covariates can have nested factors and rank deficiency; inspect design",
       "High Poisson variance may cause overdispersion; NB alpha=1 is fixed as sensitivity, not inferred",
       "No direct biomass, phenology interactions, or independent country replication"]}
    try:
        original,rows,invalid,provenance=prepare()
        D=matrix(rows).copy()
        D["floral_unit_count"]=np.exp(D["log_floral_units"])
        D["hour"]=D["hour"].fillna(D["hour"].median())
        for col in ("target_flower_corrected","floral_unit","flower_structure",
                    "flower_cover","flower_context","sunshine","wind_speed",
                    "cloud_cover","land_cover","habitat_type"):
            D[col]=D[col].fillna("__MISSING__").astype(str)
        if len(D)<3000:raise ValueError("Missing valid input")
        q99=float(D["floral_unit_count"].quantile(.99))
        variants={
            "all_adjusted_habitat":(D,FORMULA_ENV),
            "all_within_1km_square":(D,FORMULA_WITHIN),
            "exclude_other_flower_target":(
                D.loc[~D["target_flower_corrected"].str.contains("Other",case=False)].copy(),
                FORMULA_WITHIN),
            "drop_top_1pct_flower_units":(D.loc[D["floral_unit_count"]<=q99].copy(),
                                           FORMULA_WITHIN),
        }
        report["original_rows"]=len(original)
        report["source"]=provenance
        report["source_rows_retained"]=len(D)
        report["source_record_drop_reason"]=dict(invalid)
        report["variants"]=list(map(lambda t:rowsummary(t[1][0],t[0]),variants.items()))
        outcomes=[]
        for variant,(data,form) in variants.items():
            targets=GUILDS if variant in ("all_adjusted_habitat","all_within_1km_square") else GUILDS[:3]
            for response in targets:
                if response not in data.columns:raise ValueError("Missing outcome "+response)
                results=fit_one(data,form,response,variant,"poisson")
                outcomes.append(results)
                print("POISSON",variant,response,results["log_units_beta"],results["cluster_95pct_ci"],
                      results["doubling_flower_units_visits_per_unit_multiplier"],flush=True)
        # NB sensitivity (fixed alpha); the same observed square clusters.
        for y in GUILDS[:3]:
            res=fit_one(D,FORMULA_ENV,y,"all_adjusted_habitat","negbin_alpha1")
            outcomes.append(res)
            print("NEG_BIN",y,res["log_units_beta"],res["cluster_95pct_ci"],flush=True)
        report["estimates"]=outcomes
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
        print("ERROR",report["error"],flush=True)
    (OUT/"flower_units_diminishing_return_summary.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    if report["status"]!="success":raise RuntimeError(report["error"])

if __name__=="__main__":main()
