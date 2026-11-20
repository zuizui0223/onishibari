#!/usr/bin/env python3
"""Does floral supply reduce ZERO-visitor PoMS FIT observation windows?

Uses the same UK 2017–2023 1km-square allocated PoMS dataset as elasticity
pilot, but response is P(at least one insect in a selected 50cm quadrat
over 10 minutes). This is NOT probability that every single flower receives
a visit nor effective compatible pollen deposition.
"""
from pathlib import Path
import json,math
import numpy as np
import patsy
import scipy.stats as ss
import statsmodels.api as sm
from ukpoms_flower_vs_landcover_cv import prepare,matrix
from ukpoms_flower_visit_elasticity import FORMULA_ENV

OUT=Path("analysis_results_zero_visit")
OUT.mkdir(exist_ok=True)
def fit(rows,y,design):
    X=patsy.dmatrix(design,rows,return_type="dataframe",NA_action="raise")
    outcome=(rows[y].to_numpy(float)>0).astype(int)
    groups=rows["square_id"].astype(str)
    z=sm.GLM(outcome,X,family=sm.families.Binomial()).fit(
        cov_type="cluster",cov_kwds={"groups":groups},maxiter=70)
    beta=float(z.params["log_floral_units"])
    se=float(z.bse["log_floral_units"])
    if not (math.isfinite(beta) and math.isfinite(se) and se>0):raise ValueError("invalid coefficient")
    return {"response":y,"n":len(rows),"n_squares":int(groups.nunique()),
      "proportion_windows_with_zero":round(float(np.mean(outcome==0)),5),
      "logistic_log_units_beta":round(beta,6),
      "cluster95_ci":[round(beta-1.96*se,6),round(beta+1.96*se,6)],
      "p_beta_eq_zero":round(float(2*ss.norm.sf(abs(beta/se))),8),
      "odds_ratio_any_visit_when_flower_units_double":round(2**beta,6),
      "outcome_unit":"at least one insect landing on selected flower target in fixed 10-min quadrat"}
def main():
    report={"source_doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "scope":"10 minute 50cm selected flower FIT observations in allocated UK squares 2017–2023",
      "hypothesis":"More local floral units reduces chance of all-zero timed observation even if visits per unit dilute",
      "status":"failure",
      "limitations":["Flowering quadrats are selectively sampled; do not infer scarcity in unobserved non-flowering habitat",
         "Binary total-visit does not imply pollinator effectiveness, pollen deposition or seed set",
         "Flower heads vs spikes/umbels differ; structure and target identity controlled but not mechanically equivalent",
         "Land cover, site type, date and weather modeled as covariates; residual confounding",
         "No site fixed effects due complete separation for many all-positive sites; errors clustered by 1km-square"]}
    try:
        _,rows,_,provenance=prepare()
        D=matrix(rows).copy()
        D["hour"]=D["hour"].fillna(D["hour"].median())
        report["source"]=provenance
        report["models"]={}
        for y in ("all_insects_total","all_bees","hoverflies","bumblebees"):
            report["models"][y]=fit(D,y,FORMULA_ENV)
            print("ANY_VISIT",y,report["models"][y],flush=True)
        report["status"]="success"
    except Exception as ex:
        report["error"]=type(ex).__name__+": "+str(ex)
        print("FAILED",report["error"],flush=True)
    (OUT/"zero_visitor_visit_probability.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    if report["status"]!="success":raise RuntimeError(report["error"])
if __name__=="__main__":main()
