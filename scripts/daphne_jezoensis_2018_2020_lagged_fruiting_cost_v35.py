#!/usr/bin/env python3
"""v0.35 | A two-year-lag diagnostic of the cost of female fruiting in same genets.

Daphne jezoensis Shibata & Kudo 2023 AJB source from Zenodo 8253921.
This does NOT causally measure a cost of reproduction: 2019 unsampled,
plant survivorship/size/resources unknown, HF outcome-conditioning, and
natural fruiting nonrandom. Stable exactly H and F genets are treated
separately to avoid conflating sex-code transitions with life history.

Outcomes natural 2020 flower number; predictors 2018 natural fruit presence,
2018 natural flower number (proxy for baseline size). Diagnostic:
log1p(flow2020) ~ intercept + log1p(flow2018) + I(fruit2018>0).
Genet bootstrap and Freedman-Lane residual permutation are exploratory.
Only group aggregates and coefficients, NO raw IDs, data rows or coordinates.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np,pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
OUT=Path("analysis_results_daphne_lagged_fruiting_v35")
RNG=np.random.default_rng(20261117)
BOOT=10000
PERM=7000

def s(v):
    x=float(v)
    return round(x,7) if np.isfinite(x) else None
def prep(ad,sex):
    d=ad[(ad["Sex_2018.g"]==sex)&(ad["Sex_2020.g"]==sex)].copy()
    d["fl18"]=pd.to_numeric(d["NaturalFlowerNo_2018.g"],errors="coerce")
    d["fl20"]=pd.to_numeric(d["NaturalFlowerNo_2020.g"],errors="coerce")
    d["fr18"]=pd.to_numeric(d["NaturalFruitNo_2018.g"],errors="coerce")
    d["fr20"]=pd.to_numeric(d["NaturalFruitNo_2020.g"],errors="coerce")
    a=(d.fl18.gt(0)&d.fl20.gt(0)&d.fr18.notna()&d.fr20.notna()&
       d.fr18.ge(0)&d.fr20.ge(0)&d.fr18.le(d.fl18)&d.fr20.le(d.fl20))
    d=d[a].copy()
    assert len(d)==(45 if sex=="H" else 19)
    d["fruited18"]=d.fr18.gt(0).astype(int)
    d["logf18"]=np.log1p(d.fl18)
    d["logf20"]=np.log1p(d.fl20)
    d["delta_log"]=d.logf20-d.logf18
    assert d.fruited18.nunique()==2
    return d

def ols(X,y):
    coef=np.linalg.lstsq(X,y,rcond=None)[0]
    return coef
def analyze(d,sex):
    y=d.logf20.to_numpy(float)
    x=d.logf18.to_numpy(float)
    a=d.fruited18.to_numpy(float)
    n=len(d)
    X=np.column_stack([np.ones(n),x,a])
    X0=np.column_stack([np.ones(n),x])
    bet=ols(X,y)
    residual=y-X@bet
    m0=ols(X0,y);fit0=X0@m0
    e0=y-fit0
    # With observational fruit status correlated with prior flower count, only a
    # conditional predictive diagnostic is described, not randomization of fruits.
    boot=[]
    for i in range(BOOT):
        idx=RNG.integers(0,n,size=n)
        if len(np.unique(a[idx]))<2:continue
        if np.linalg.matrix_rank(X[idx])<3:continue
        boot.append(ols(X[idx],y[idx])[2])
    assert len(boot)>BOOT*.90
    perms=[]
    for i in range(PERM):
        yp=fit0+e0[RNG.permutation(n)]
        perms.append(ols(X,yp)[2])
    coef=float(bet[2])
    p=(1+int((np.abs(perms)>=abs(coef)-1e-12).sum()))/(len(perms)+1)
    group=[]
    for f in (0,1):
        g=d[d.fruited18.eq(f)]
        group.append({
          "fruited_2018":bool(f),"n_genets":len(g),
          "mean_2018_flower_number":s(g.fl18.mean()),
          "median_2018_flower_number":s(g.fl18.median()),
          "mean_2020_flower_number":s(g.fl20.mean()),
          "median_2020_flower_number":s(g.fl20.median()),
          "mean_log_flower_change_2020_minus_2018":s(g.delta_log.mean()),
          "mean_2020_fruit_rate":s((g.fr20/g.fl20).mean()),
          "n_fruited_2020":int(g.fr20.gt(0).sum())})
    loo=[]
    # Six 2018 fruiting H can be highly influential.
    for i in np.where(a==1)[0]:
        ix=np.arange(n)!=i
        if np.linalg.matrix_rank(X[ix])==3:
            loo.append(float(ols(X[ix],y[ix])[2]))
    out={
      "sex_exact_code_both_years":sex,
      "n_genets":n,
      "n_fruited_2018":int((a>0).sum()),
      "baseline_adjusted_fruited2018_coefficient_log1p_flower2020":s(coef),
      "year_2018_log1p_flower_coefficient":s(bet[1]),
      "fruit2018_binary_coefficient_bootstrap_95%":[s(v) for v in np.quantile(boot,[.025,.975])],
      "bootstrap_n_valid":len(boot),
      "Freedman_Lane_residual_permutation_p_two_sided_EXPLORATORY":s(p),
      "permutation_n":len(perms),
      "leave_one_fruited2018_genet_out_coefficient_range":
           [s(min(loo)),s(max(loo))] if loo else None,
      "within_group_descriptive":group,
      "plant_flower_change_logratio_unadjusted_fruited_minus_nonfruited":
          s(d.loc[d.fruited18.eq(1),"delta_log"].mean()-d.loc[d.fruited18.eq(0),"delta_log"].mean()),
      "limits":{
        "two_year_interval_has_unobserved_2019_fruiting_and_floral_output":True,
        "baseline_flower_number_is_imperfect_size_and_resource_proxy":True,
        "flower_number_can_regress_to_mean_between_years":True,
        "survivorship_conditioned_on_flowering_both_years":True,
        "not_random_assignment_of_fruit_production":True,
        "not_total_floral_construction_cost_or_paternal_fitness":True,
        "not_a_proof_of_no_cost_even_if_interval_includes_zero":True}
    }
    return out

def main():
    OUT.mkdir(exist_ok=True)
    files=data()
    ad=files["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(ad)==492
    h=prep(ad,"H");f=prep(ad,"F")
    out={"status":"original_Daphne_lagged_fruit_to_future_flower_diagnostic_v35",
      "source_doi":"10.5281/zenodo.8253921",
      "paper_doi":"10.1002/ajb2.16239",
      "source_2022_experimental_prior_doi":"10.1111/plb.13458",
      "sex_stratified":[analyze(h,"H"),analyze(f,"F")],
      "do_not_claim_reproductive_costs_causally":True,
      "source_study_same_as_v34_v33":True,
      "no_D_pseudomezereum_observations":True}
    (OUT/"daphne_H_F_2018_fruit_2020_flower_cost_sensitivity_v35.json").write_text(
       json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_LAGGED_FEMALE_COST_V35",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
