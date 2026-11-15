#!/usr/bin/env python3
"""Preliminary two-part Zurich 2016 floral display–visitor rate audit.

Part A observed P(any captured insect on a flower-bearing sampling day).
Part B visitor count RATE only among observed-positive days.
This is NOT a fully specified zero-truncated hurdle likelihood; the
positive-only count fit is a descriptive sensitivity regression.

Separate plant species, 24-garden cluster robust uncertainty. Actual flowers
and observation effort are matched on the SAME site/date.
No individual-genet flower biomass, pollen efficacy or seed success inference.
"""
from collections import Counter
import json,math
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as ss
from zurich2016_flower_date_link import prepare,PLANTS

OUT=Path("analysis_results_zurich_two_part")
OUT.mkdir(exist_ok=True)
BINS=[(1,1),(2,3),(4,9),(10,999999)]

def binomial_interval(k,n):
    if not n:return None
    lo=0 if k==0 else ss.beta.ppf(.025,k,n-k+1)
    hi=1 if k==n else ss.beta.ppf(.975,k+1,n-k)
    return [round(float(lo),5),round(float(hi),5)]

def group_desc(rows):
    bins=[]
    for lo,hi in BINS:
        active=[r for r in rows if lo<=r["flower_units"]<=hi]
        n=len(active)
        zero=sum(r["all_visitor_count"]==0 for r in active)
        pos=n-zero
        bins.append({"flower_units_range":[lo,hi if hi<100000 else None],
            "garden_days":n,"n_zero_days":zero,
            "n_any_visit_days":pos,
            "any_visitor_fraction":round(pos/n,5) if n else None,
            "any_visitor_fraction_exact95":binomial_interval(pos,n),
            "median_visit_count_if_positive":round(float(np.median([r["all_visitor_count"] for r in active if r["all_visitor_count"]>0])),4) if pos else None})
    return bins

def fit(rows,day_adjust=False):
    d=pd.DataFrame(rows)
    d["logF"]=np.log(d["flower_units"].astype(float))
    d["doy"]=d["day_of_year"].astype(float)
    d["logEff"]=np.log(d["effort_minutes"].astype(float)/60)
    d["site"]=d["site"].astype(str)
    d["binary"]=np.where(d["all_visitor_count"]>0,1,0)
    d["doyc"]=(d["doy"]-d["doy"].mean())/10
    X=sm.add_constant(d[["logF"]+["doyc" if day_adjust else "logEff"]])
    if day_adjust:
        X["logEff"]=d["logEff"]
    if len(set(d["binary"]))==1:return None
    z=sm.GLM(d["binary"],X,family=sm.families.Binomial()).fit(
        maxiter=100,tol=1e-9,cov_type="cluster",cov_kwds={"groups":d["site"]})
    beta=float(z.params["logF"]);se=float(z.bse["logF"])
    if not (math.isfinite(beta) and math.isfinite(se) and se>0 and abs(beta)<15):
        # Only a single radish zero-visit date is present. Such logit
        # fits can experience quasi/complete separation and return an
        # enormous coefficient and NaN SE. Do not report that as data.
        e={"model":"any-visitor logit, not effective pollen delivery",
           "plant_source_days":len(d),"adjust_for_doy":day_adjust,
           "status":"NOT_IDENTIFIABLE_SEPARATION",
           "note":"Extreme logit parameter or nonfinite cluster SE; no odds ratio/p-value interpreted."}
    else:
        e={"model":"any-visitor logit, not effective pollen delivery",
           "plant_source_days":len(d),
           "adjust_for_doy":day_adjust,
           "flower_log_count_coefficient":round(beta,6),
           "cluster_se":round(se,6),
           "cluster95":[round(beta-1.96*se,6),round(beta+1.96*se,6)],
           "odds_any_visit_when_F_double":round(2**beta,6),
           "p_beta0":round(float(2*ss.norm.sf(abs(beta/se))),7),
           "status":"ESTIMABLE"}
    positive=d[d["binary"]==1].copy()
    if len(positive)<15:return e
    X2=sm.add_constant(positive[["logF"]+["doyc" if day_adjust else "logEff"]])
    if day_adjust:X2["logEff"]=positive["logEff"]
    # Deliberately a positive-only descriptive QMLE; not zero-truncated.
    p=sm.GLM(positive["all_visitor_count"],X2,family=sm.families.Poisson(),
         offset=positive["logEff"]).fit(
           maxiter=100,cov_type="cluster",cov_kwds={"groups":positive["site"]})
    b=float(p.params["logF"]);se2=float(p.bse["logF"])
    e["positive_only_count_QMLE"]={
       "n_positive_garden_days":len(positive),
       "flower_log_count_coefficient":round(b,6),
       "cluster95":[round(b-1.96*se2,6),round(b+1.96*se2,6)],
       "visitor_total_multiplier_when_F_double":round(2**b,6),
       "visits_per_unit_multiplier_when_F_double":round(2**(b-1),6),
       "caution":"Positive-only Poisson QMLE is not the zero-truncated count part of a proper hurdle likelihood."}
    return e

def main():
    report={"source":"2016 Zurich phytometer observer-day visits and flowers, Chacko et al 2026 JAE",
        "status":"failed","purpose":"Distinguish zero-visit occurrence from visit count on observed-positive flower-bearing days",
        "scope":"24 gardens, four standardized non-conspecific ornamental/wild phytometers",
        "cautions":[
          "Florally selected days only, small within-species site/date sample",
          "F>=1 and observer minutes>=1 already enforced, results not applicable to nonflowering populations",
          "Flower-unit structure differs: carrot umbels vs other species",
          "Binned descriptive thresholds post hoc, not prespecified ecological threshold",
          "The logit and positive-only Poisson fits are not a coherent fitted hurdle likelihood",
          "Garden cluster-robust SE with only 24 clusters; seasonal temperature/weather vary",
          "No causal manipulation or actual conspecific pollen success / floral dry mass / lifetime fitness"]}
    try:
        rows,qa,source=prepare()
        report["source_checks"]=qa
        report["source"]=source
        report["per_plant"]={}
        for plant in PLANTS:
            rr=[r for r in rows if r["plant"]==plant]
            result={"n_garden_date_observations":len(rr),
                "n_zero":sum(r["all_visitor_count"]==0 for r in rr),
                "binned_observational":group_desc(rr),
                "logistic_and_positive_count":{}}
            for mark in [False,True]:
                name="with_day_of_year" if mark else "without_day_of_year"
                try: result["logistic_and_positive_count"][name]=fit(rr,mark)
                except Exception as ex:
                    result["logistic_and_positive_count"][name]={"error":type(ex).__name__+": "+str(ex)[:400]}
            report["per_plant"][plant]=result
            print("TWO_PART",plant,result,flush=True)
        report["status"]="success"
    except Exception as ex:
        report["error"]=type(ex).__name__+": "+str(ex)
        print("FAILED",report["error"],flush=True)
    (OUT/"flower_resource_zero_vs_positive_descriptive.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
