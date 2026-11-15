#!/usr/bin/env python3
"""Does floral-unit abundance shift hoverfly-vs-bee visitor composition?

Each PoMS FIT count is 10min/0.25m2. Outcome is two-column binomial counts
[hoverflies, all bees] from the SAME observation, conditional on at least one
recorded bee or hoverfly. Positive slope = increasing hoverfly odds relative
to bee odds with more floral units (NOT any causal preference or effectiveness).
Poisson guild-specific slopes are studied separately.

Inferences: marginal association, not mechanistic competition or pollination.
"""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
import patsy
import scipy.stats as ss
import statsmodels.api as sm
from ukpoms_flower_vs_landcover_cv import prepare,matrix
from ukpoms_flower_visit_elasticity import FORMULA_ENV,FORMULA_WITHIN

OUT=Path("analysis_results_relative_guild")
OUT.mkdir(exist_ok=True)
def one(subset,formula,label):
    d=subset.copy()
    d["grouptotal"]=d["all_bees"]+d["hoverflies"]
    d=d[d["grouptotal"]>0].copy()
    if len(d)<400:raise ValueError("Too few mixed-guild rows")
    if (d["all_bees"]<0).any() or (d["hoverflies"]<0).any():
        raise ValueError("Negative visitor counts")
    X=patsy.dmatrix(formula,d,return_type="dataframe",NA_action="raise")
    y=np.column_stack((d["hoverflies"].to_numpy(float),d["all_bees"].to_numpy(float)))
    model=sm.GLM(y,X,family=sm.families.Binomial())
    fit=model.fit(maxiter=100,tol=1e-8,method="irls",
       cov_type="cluster",cov_kwds={"groups":d["square_id"].astype(str)})
    slope=float(fit.params["log_floral_units"])
    se=float(fit.bse["log_floral_units"])
    if not (math.isfinite(slope) and math.isfinite(se) and se>0):
        raise ValueError("Invalid slope/SE")
    return {"subset":label,"n_records_with_any_hoverfly_or_bee":len(d),
        "n_square_clusters":int(d["square_id"].nunique()),
        "source_unit":"relative hoverfly-versus-bee landing odds within same 10-minute quadrat record",
        "beta_logit_hoverfly_odds_vs_bees_per_ln_flower_units":round(slope,6),
        "cluster_95pct_ci":list(map(lambda x:round(x,6),
                          [slope-1.96*se,slope+1.96*se])),
        "cluster_p_two_sided_beta_eq0":round(float(2*ss.norm.sf(abs(slope/se))),8),
        "hoverfly_to_bee_odds_multiplier_if_floral_units_double":round(float(2**slope),6),
        "total_hoverfly_events":int(y[:,0].sum()),"total_bee_events":int(y[:,1].sum()),
        "controls":"target taxon/structure, weather, season + habitat OR 1km-site fixed effects",
        "warning":"A binomial share conditions on any counted bees/hoverflies and cannot identify habitat-driven total density."}

def main():
    report={"status":"failed","source_doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "scientific_question":"Do bees and hoverflies differ in their relative scaling with local flower units?",
      "target_pair":"hoverflies vs bumblebees+honeybees+solitary bees",
      "measurement":"plant visits per 10 minute 50-cm quadrat in allocated 101-square UK PoMS 2017-23",
      "limits":["Conditional group composition, not pollen deposition/plant seed set",
        "Different plant flower heads/umbels are not exchangeable with individual flowers",
        "Visitor classification confidence and habitat structure may covary with flower count",
        "Fixed-square/target adjustment and 101-cluster robust SE are observational, not causal",
        "Site data and raw locations never exported"]}
    try:
        _,cleaned,_,source=prepare()
        report["source"]=source
        D=matrix(cleaned).copy()
        D["floral_unit_count"]=np.exp(D["log_floral_units"])
        D["hour"]=D["hour"].fillna(D["hour"].median())
        for v in ("target_flower_corrected","floral_unit","flower_structure",
           "flower_cover","flower_context","sunshine","wind_speed",
           "cloud_cover","land_cover","habitat_type"):
            D[v]=D[v].fillna("__MISSING__").astype(str)
        q99=float(D["floral_unit_count"].quantile(.99))
        variants=[
          ("adjusted_habitat_all",D,FORMULA_ENV),
          ("within_square_all",D,FORMULA_WITHIN),
          ("exclude_other_target",D.loc[~D["target_flower_corrected"].str.contains("Other",case=False)].copy(),FORMULA_ENV),
          ("exclude_top_1pct",D.loc[D["floral_unit_count"]<=q99].copy(),FORMULA_WITHIN)]
        report["estimates"]=[]
        for label,subset,design in variants:
            try:
                res=one(subset,design,label)
                report["estimates"].append(res)
                print("GUILD",label,res["beta_logit_hoverfly_odds_vs_bees_per_ln_flower_units"],
                      res["cluster_95pct_ci"],flush=True)
            except Exception as ex:
                report.setdefault("failed_sensitivities",[]).append(
                    {"subset":label,"error":type(ex).__name__+": "+str(ex)[:350]})
                print("FAILED",label,type(ex).__name__,str(ex)[:160],flush=True)
        report["status"]="success" if any(x["subset"]=="adjusted_habitat_all" for x in report["estimates"]) else "failed"
    except Exception as ex:
        report["error"]=type(ex).__name__+": "+str(ex)
    (OUT/"hoverfly_vs_bee_share_sensitivity.json").write_text(
       json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error","No primary fit"))
if __name__=="__main__":main()
