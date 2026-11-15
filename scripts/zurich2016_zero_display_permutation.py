#!/usr/bin/env python3
"""Within-garden conditional randomization test: are larger flower displays
more likely to receive >=1 captured visitor on same garden×plant day?

24 Zurich gardens, 2016 phyto-experiment. Only gardens with at least one
zero day AND one positive-visit day are informative. Within each such garden,
keep exact number of zero visits; randomly shuffle day labels 10,000 times.
Statistic = equal-garden average of mean log flower units on positive days
minus mean log flower units on zero days. Compare raw log F and log F
detrended by simple global plant-specific DOY+DOY² OLS.

This does NOT test pollen efficiency, seed fitness, environmental causal
selection, or reliable multiyear visit availability.
"""
from collections import defaultdict
import json,math
from pathlib import Path
import numpy as np
from zurich2016_flower_date_link import prepare,PLANTS

OUT=Path("analysis_results_within_garden_any_visit")
OUT.mkdir(exist_ok=True)
N_RANDOM=10000

def prepare_study(rows,plant,use_doy_residual):
    d=[x for x in rows if x["plant"]==plant]
    F=np.log([x["flower_units"] for x in d])
    doy=np.array([x["day_of_year"] for x in d],float)
    if use_doy_residual:
        doy=(doy-doy.mean())/10.
        X=np.column_stack([np.ones(len(F)),doy,doy**2])
        fit=X@np.linalg.lstsq(X,F,rcond=None)[0]
        F=F-fit
    grouped=defaultdict(list)
    for r,x in zip(d,F):
        grouped[r["site"]].append((float(x),int(r["all_visitor_count"]>0)))
    mixed=[]
    for site,a in grouped.items():
        x=np.array([t[0] for t in a],float)
        y=np.array([t[1] for t in a],int)
        if (y==0).all() or (y==1).all():continue
        mixed.append((x,y))
    return mixed,len(grouped)

def test(groups,seed):
    RNG=np.random.default_rng(seed)
    if len(groups)<2:return {"status":"not_enough_informative_gardens","n_informative_gardens":len(groups)}
    def effect(x,y):
        return float(x[y==1].mean()-x[y==0].mean())
    observed=float(np.mean([effect(x,y) for x,y in groups]))
    null=np.zeros(N_RANDOM)
    for j in range(N_RANDOM):
        z=[]
        for x,y in groups:
            perm=RNG.permutation(y)
            z.append(effect(x,perm))
        null[j]=float(np.mean(z))
    return {
        "status":"computed","n_informative_gardens":len(groups),
        "n_informative_garden_date_observations":sum(len(x) for x,y in groups),
        "effect_observed_mean_positive_minus_zero_log_flower_units":round(observed,6),
        "exp_effect_positive_vs_zero_flower_unit_ratio":round(math.exp(observed),6),
        "null_mean":round(float(null.mean()),6),
        "null_randomization_95_percentile":list(map(lambda x:round(float(x),6),np.quantile(null,[.025,.975]))),
        "one_sided_p_positive_more_flower_units":round((1+int((null>=observed).sum()))/(1+N_RANDOM),6),
        "two_sided_p":round((1+int((np.abs(null)>=abs(observed)).sum()))/(1+N_RANDOM),6),
        "n_permutations":N_RANDOM,
        "null":"Exchange any-visit labels across flower-bearing sampled dates within same garden+plant, preserving visit-positive date count"}
def main():
    report={"source":"Zurich 2016 24 gardens 4 phytometer species, doi:10.16904/envidat.676",
      "status":"failed",
      "goal":"Separate within-garden flower display and same-day any-visit chance from site differences",
      "important_limitations":[
        "Conditional on sites with both zero and nonzero observation days; not a population prevalence estimate",
        "Observer minutes vary; this permutation does NOT adjust exposure beyond filtering positive recording minutes",
        "Date-adjusted detrending is approximate common seasonal quadratic, not garden-specific weather or flowering phase",
        "Small group count within each plant: permutation p values exploratory",
        "Visitor landings not pollen deposition; no biomass or reproduction modeled",
        "Do not call effect causation or a threshold in an individual plant"]}
    try:
        rows,qa,source=prepare()
        report["source"]=source
        report["by_plant"]={}
        for i,p in enumerate(PLANTS):
            gro,n=prepare_study(rows,p,False)
            det,m=prepare_study(rows,p,True)
            report["by_plant"][p]={
                "flower_active_gardens":n,
                "raw_within_garden":test(gro,seed=20261115+i),
                "after_global_quadratic_DOY_detrending":test(det,seed=20261121+i)}
            print("WITHIN",p,report["by_plant"][p],flush=True)
        report["status"]="success"
    except Exception as ex:
        report["error"]=type(ex).__name__+": "+str(ex)
        print("FAILED",report["error"],flush=True)
    (OUT/"any_visit_site_conditional_permutation.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
