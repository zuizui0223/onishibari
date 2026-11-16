#!/usr/bin/env python3
"""v0.36: common-support and reverse-time diagnostics of an apparent H fruit-cost.

Builds on v35 original-genet 2018/2020 fruit/flower data.
Never interpret observational intervals as a causal reproductive cost.
A reverse-time placebo is a diagnostic for static genotype/resource confounding,
not a valid randomized negative control. Shared-source independent IDs only.
"""
import json
from pathlib import Path
import numpy as np,pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
from daphne_jezoensis_2018_2020_lagged_fruiting_cost_v35 import prep,ols
OUT=Path("analysis_results_daphne_lagged_fruiting_identification_v36")
RNG=np.random.default_rng(20261117)
BOOTS=5000

def val(x):
    x=float(x)
    return round(x,6) if np.isfinite(x) else None
def desc(x):
    x=np.asarray(x,float)
    return {"n":len(x),"min":val(np.min(x)),"q25":val(np.quantile(x,.25)),
      "median":val(np.median(x)),"q75":val(np.quantile(x,.75)),"max":val(np.max(x)),
      "mean":val(np.mean(x))}
def coef(y,baseline,fruit):
    x=np.column_stack([np.ones(len(y)),np.log1p(baseline),fruit])
    return float(ols(x,np.log1p(y))[2])
def analyze(d,sex):
    a=d[d.fruited18.eq(1)].copy()
    b=d[d.fruited18.eq(0)].copy()
    assert len(a)>0 and len(b)>0
    lower=max(a.fl18.min(),b.fl18.min())
    upper=min(a.fl18.max(),b.fl18.max())
    feasible=(d.fl18>=lower)&(d.fl18<=upper)
    common=d[feasible].copy()
    a2=common[common.fruited18.eq(1)]
    b2=common[common.fruited18.eq(0)]
    # Within observed common-support range; no extrapolation for the other genets.
    regression_in_support=coef(common.fl20.to_numpy(float),common.fl18.to_numpy(float),
       common.fruited18.to_numpy(float)) if len(a2)>=2 and len(b2)>=3 else None
    # Greedy rank nearest control with replacement on log1p baseline number.
    train=b.logf18.to_numpy(float)
    matches=[]
    for _,r in a.iterrows():
        delta=np.abs(train-float(r.logf18))
        idx=int(np.argmin(delta))
        control=b.iloc[idx]
        matches.append({
          "abs_log_flower_baseline_difference":float(delta[idx]),
          "treated_minus_control_log_2020_flowers":float(r.logf20-control.logf20)
        })
    matched={"n_fruiting_genets":len(a),
      "n_with_nearest_control_differing_by_lt_20pct_in_baseline_flower_number":
       sum(m["abs_log_flower_baseline_difference"]<=np.log(1.2) for m in matches),
      "n_with_nearest_control_differing_by_lt_50pct":
       sum(m["abs_log_flower_baseline_difference"]<=np.log(1.5) for m in matches),
      "max_nearest_log_flower_number_gap":val(max(m["abs_log_flower_baseline_difference"] for m in matches)),
      "mean_matched_difference_log_flower2020_with_replacement_no_uncertainty":
       val(np.mean([m["treated_minus_control_log_2020_flowers"] for m in matches]))
    }
    # "Reverse time" regression: 2020 fruit state predictor of 2018 flowers,
    # controlling 2020 flower number. Any association cannot arise as an effect
    # of fruit2020 on flowers2018; still observational/cross-time selection.
    revfruit=d.fr20.gt(0).to_numpy(float)
    reverse=coef(d.fl18.to_numpy(float),d.fl20.to_numpy(float),revfruit)
    # Future flower count effect from 2018 fruit with prior size offset conditional.
    forward=coef(d.fl20.to_numpy(float),d.fl18.to_numpy(float),d.fruited18.to_numpy(float))
    # Paired label switches counted, not a true stochastic comparison.
    result={"sex_exact_both_years":sex,
      "n_genets":len(d),
      "2018_fruit_positive_n":len(a),
      "2020_fruit_positive_n":int(d.fr20.gt(0).sum()),
      "2018_flower_number_if_fruited":desc(a.fl18),
      "2018_flower_number_if_not_fruited":desc(b.fl18),
      "2020_flower_number_if_fruited2018":desc(a.fl20),
      "2020_flower_number_if_not_fruited2018":desc(b.fl20),
      "common_support_2018_flower_range":[int(lower),int(upper)] if lower<=upper else None,
      "n_genets_common_support":len(common),
      "n_fruiting_genets_common_support":len(a2),
      "n_nonfruiting_genets_common_support":len(b2),
      "common_support_regression_coef":val(regression_in_support) if regression_in_support is not None else None,
      "nearest_baseline_flower_match_diagnostic":matched,
      "baseline_adjusted_forward_fruit2018_to_flower2020_coef":val(forward),
      "reverse_time_fruit2020_to_prior_flower2018_coef_not_causal":val(reverse),
      "n_years_observed":2,
      "n_calendar_years_between":2,
      "n_missing_intermediate_2019_census":1,
      "cannot_infer_effect_on_fatherhood_or_genet_survival":True}
    return result
def main():
    OUT.mkdir(exist_ok=True)
    ad=data()["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    out={"status":"source_original_Daphne_cost_identification_common_support_v36",
      "source_2023_paper":"10.1002/ajb2.16239",
      "author_data":"10.5281/zenodo.8253921",
      "resource_experiment_2022":"10.1111/plb.13458",
      "by_sex":[analyze(prep(ad,sex),sex) for sex in ("H","F")],
      "interpretation":{
        "observational_forward_association_not_direct_resource_cost":True,
        "reverse_year_diagnostic_not_true_randomized_negative_control":True,
        "unobserved_2019_and_survivor_selection_cannot_be_fixed_by_matching":True,
        "same_adult_table_study_as_v34_v35":True,
        "no_compatible_seeds_or_lifetime_male_fitness":True,
        "not_focal_D_pseudomezereum":True}}
    (OUT/"daphne_2018fruit_2020flower_common_support_v36.json").write_text(
      json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_FRUIT_COST_IDENTIFICATION_V36",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
