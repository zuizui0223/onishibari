#!/usr/bin/env python3
"""Exploratory heterogeneity of floral-unit visitor saturation by accessibility.

PoMS allocated UK 1-km FIT 2017-2023, 10-minute fixed quadrat.
Predeclared moderators:
(A) target flower structure (open/closed) (BETWEEN flower types);
(B) surrounding flower context: isolated/same-flower/mixed-flower patch.
Fit Poisson quasi ML with site, target, climate/season/observer weather FE,
cluster-robust SE by 1-km square. No per-flower stigma pollen or seed data.
Posthoc ecological comparisons not designed to establish causation.
"""
from collections import Counter
import json,math
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import chi2
from ukpoms_flower_vs_landcover_cv import prepare,matrix
from ukpoms_flower_visit_elasticity import FORMULA_WITHIN

OUT=Path("analysis_results_flower_context_interactions")
OUT.mkdir(exist_ok=True)
GUILDS=("all_insects_total","all_bees","hoverflies")
MODERATORS={
 "flower_structure":("open","closed"),
 "flower_context":(
   "More or less isolated",
   "Growing in a larger patch of the same flower",
   "Growing in a larger patch of many different flowers"),
}

def slope_contrast(design_info,first,level,moderator):
    frame=first.copy()
    frame[moderator]=level
    frame["log_floral_units"]=0.
    x0=np.asarray(patsy.build_design_matrices([design_info],frame)[0])[0]
    frame["log_floral_units"]=1.
    x1=np.asarray(patsy.build_design_matrices([design_info],frame)[0])[0]
    return x1-x0

def fit_model(df,moderator,levels,response):
    D=df[df[moderator].isin(levels)].copy()
    if len(D)<500:raise ValueError("Too few rows")
    form=FORMULA_WITHIN+" + log_floral_units:C("+moderator+")"
    X=patsy.dmatrix(form,D,return_type="dataframe",NA_action="raise")
    y=D[response].to_numpy(float)
    model=sm.GLM(y,X,family=sm.families.Poisson())
    m=model.fit(maxiter=80,tol=1e-8,cov_type="cluster",
       cov_kwds={"groups":D["square_id"].astype(str)})
    betas=m.params.to_numpy()
    vc=m.cov_params().to_numpy()
    contrasts=[]
    estimates={}
    for level in levels:
        first=D.iloc[[0]].copy()
        contrast=slope_contrast(X.design_info,first,level,moderator)
        coef=float(contrast@betas)
        se=math.sqrt(max(0,float(contrast@vc@contrast)))
        if se<=0:raise ValueError("invalid contrast variance")
        estimates[level]={
            "n_source_cases":int(sum(D[moderator]==level)),
            "elasticity":round(coef,6),
            "cluster_95pct_CI":[round(coef-1.96*se,6),round(coef+1.96*se,6)],
            "doubling_floral_unit_visits_per_unit_multiplier":round(2**(coef-1),6)}
        contrasts.append(contrast)
    C=np.stack([c-contrasts[0] for c in contrasts[1:]])
    cov=C@vc@C.T
    est=C@betas
    stat=float(est@np.linalg.pinv(cov)@est)
    df_eff=int(np.linalg.matrix_rank(cov))
    pval=float(chi2.sf(stat,df_eff)) if df_eff else None
    return {"moderator":moderator,"outcome":response,
        "observations":len(D),"n_1km_squares":int(D["square_id"].nunique()),
        "levels":estimates,
        "joint_heterogeneity_chi2":round(stat,5),
        "joint_p_cluster_robust":round(pval,7) if pval is not None else None,
        "joint_effective_df":df_eff,
        "warnings":[
          "Flower structure/open/closed is often fixed within a flower species, so comparisons are among taxa, not within-taxon manipulations",
          "Flower context is survey categorical and may reflect unmeasured patch size, local flower density, site management",
          "Site and target fixed effects limit but do not remove sampling and ecosystem confounding",
          "The chi-square interaction test is exploratory across multiple moderators/outcomes"]}

def main():
    report={"source":"UK PoMS 2017–23, 1km-square allocated FIT observations",
      "doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "goal":"Does flower context or flower opening morphology moderate sublinear flower-unit visitor scaling?",
      "status":"failed","results":[],"test_family_n":6}
    try:
        _,cleaned,invalid,source=prepare()
        D=matrix(cleaned)
        D["hour"]=D["hour"].fillna(D["hour"].median())
        report["n_cases"]=len(D)
        for moderator,levels in MODERATORS.items():
            for outcome in GUILDS:
                item=fit_model(D,moderator,levels,outcome)
                report["results"].append(item)
                print("INTERACT",moderator,outcome,item["joint_p_cluster_robust"],
                      {k:v["elasticity"] for k,v in item["levels"].items()},flush=True)
        # BH adjusted across six planned modifier×guild tests.
        tests=sorted([(i,r["joint_p_cluster_robust"]) for i,r in enumerate(report["results"])],
                     key=lambda x:x[1] if x[1] is not None else 1)
        q=1
        for j in range(len(tests)-1,-1,-1):
            index,p=tests[j]
            if p is None:continue
            q=min(q,p*len(tests)/(j+1))
            report["results"][index]["joint_BH_q_across_6_tests"]=round(q,7)
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
        print("FAILED",report["error"],flush=True)
    (OUT/"floral_structure_patch_interaction.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
