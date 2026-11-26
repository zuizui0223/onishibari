#!/usr/bin/env python3
"""Independent specification stress test for Clarkia study-5 patch × donor-pollen.

The primary authors (Miller & Kay, Madrono 2025) reported an omnibus 2-df
display×patch interaction of p~0.06 in their ANOVA, so our post-hoc
one-df contrast (six intact male flowers vs six emasculated display flowers,
difference between small and large patches) must NOT be presented as a
confirmed biological mechanism based on one model.

Re-fit OLS, categorical original location as sensitivity, NB count model,
winsorized outcomes, date-block bootstrap and within-date/patch label
permutation for the one-df difference-in-differences contrast.

Source data stable SHA-checked via clarkia_concinna_five_stage_refit_v13.py.
No original plant IDs, location/site names or row-level results exported.
"""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm
import statsmodels.formula.api as smf
from clarkia_concinna_five_stage_refit_v13 import fetch_dataset,scalar

OUT=Path("analysis_results_clarkia_concinna_robustness_v13");OUT.mkdir(exist_ok=True)
SEED=20261116;RNG=np.random.default_rng(SEED);N_BOOT=5500;N_PERM=8000
TR="C(treatment)[T.6]"
TE="C(treatment)[T.6e]"
PATCH="C(patch)[T.small]"
HYP=TR+":"+PATCH+" - "+TE+":"+PATCH+" = 0"

def get_frame():
    _,allstudies=fetch_dataset()
    x=allstudies[5].copy().rename(columns={
        "display treatment":"treatment","patch size":"patch","st dia":"stem_diameter"})
    x["treatment"]=pd.Categorical(x.treatment.astype(str),categories=["1","6e","6"])
    x["patch"]=pd.Categorical(x.patch.astype(str),categories=["large","small"])
    x["date"]=x.date.astype(str)
    x["location"]=x.location.astype(str).str.lower().str.strip()
    x["grains"]=pd.to_numeric(x.grains,errors="coerce")
    x["stem_diameter"]=pd.to_numeric(x.stem_diameter,errors="coerce")
    if len(x)!=204 or x[["treatment","patch","date","grains","stem_diameter"]].isna().any().any():
        raise ValueError("study5 source unexpectedly changed")
    x["y"]=np.log1p(x.grains)
    return x

def diff_in_diff(z,y="y"):
    # Difference: (6-6e)_small - (6-6e)_large
    out={}
    for patch in ("large","small"):
        for tr in ("6","6e"):
            d=z[(z.patch==patch)&(z.treatment==tr)][y]
            if d.empty:raise ValueError("empty observed group")
            out[(patch,tr)]=float(d.mean())
    return (out["small","6"]-out["small","6e"])-(out["large","6"]-out["large","6e"])

def fit(x,formula,family=None):
    if family is None:
        m=smf.ols(formula,data=x).fit(cov_type="HC3")
    else:
        m=smf.glm(formula,data=x,family=family).fit(cov_type="HC3")
    tt=m.t_test(HYP)
    beta=float(np.squeeze(tt.effect))
    se=float(np.squeeze(tt.sd))
    p=float(np.squeeze(tt.pvalue))
    return {"contrast":"(intact6_minus_emasculated6)_small minus same_large",
      "estimate_log_scale":scalar(beta),"HC3_95":[scalar(beta-1.96*se),scalar(beta+1.96*se)],
      "p_two_sided_exploratory":scalar(p),
      "comparison_is_post_hoc_one_degree_of_freedom_not_original_omnibus":True,
      "n_source_plants":int(m.nobs)}

def block_bootstrap(x):
    original=diff_in_diff(x)
    dates=list(x.date.unique())
    chunks={d:x[x.date==d] for d in dates}
    reps=[]
    for i in range(N_BOOT):
        sampled=RNG.choice(dates,len(dates),replace=True)
        dat=pd.concat([chunks[d] for d in sampled],ignore_index=True)
        try:reps.append(diff_in_diff(dat))
        except ValueError:pass
    if len(reps)<N_BOOT-50:raise RuntimeError("unstable resampled study groups")
    return {"observed_unadjusted_difference_in_differences_log1p":scalar(original),
      "calendar_date_block_count":len(dates),
      "n_valid_date_bootstrap":len(reps),
      "date_resampling95":[scalar(t) for t in np.quantile(reps,[.025,.975])],
      "not_an_original_preregistered_contrast":True}

def permutation_within_date_patch(x):
    # Respect calendar day, environment patch, and observed allocation counts.
    # Does not remove site-driven confounding or exactly reconstruct field randomization.
    groups=[g for _,g in x.groupby(["date","patch"],observed=True)]
    original=diff_in_diff(x)
    absvals=np.empty(N_PERM)
    for i in range(N_PERM):
        shuff=[]
        for g in groups:
            gg=g.copy()
            gg["treatment"]=RNG.permutation(g.treatment.to_numpy())
            shuff.append(gg)
        new=pd.concat(shuff,ignore_index=True)
        absvals[i]=abs(diff_in_diff(new))
    p=float((1+(absvals>=abs(original)-1e-12).sum())/(N_PERM+1))
    return {"n_permutations":N_PERM,"source_strata":"date_by_patch","observed_DiD":scalar(original),
        "two_sided_stratified_label_p":scalar(p),
        "does_not_remove_natural_patch_assignment_or_site_confounds":True}

def main():
    x=get_frame()
    baseline="y ~ C(treatment)*C(patch) + stem_diameter + C(date)"
    withsite=baseline+" + C(location)"
    count="grains ~ C(treatment)*C(patch) + stem_diameter + C(date)"
    report={"status":"failed",
      "original_source_DOI":"10.5061/dryad.8pk0p2nz2",
      "reported_original_omnibus_interaction_p_approx":.06,
      "source_randomized_display_but_patch_partly_natural_and_thinned":True,
      "n_source_plants":len(x),"n_observation_days":int(x.date.nunique()),
      "n_locations":int(x.location.nunique()),"source_six_arm_n":{
         str(a)+"_"+str(b):int(len(g)) for (a,b),g in x.groupby(["treatment","patch"],observed=True)},
      "site_patch_cell_count_matrix":{
        str(k):{str(j):int(t) for j,t in g.patch.value_counts().items()}
         for k,g in x.groupby("location")},
      "models":{}}
    report["models"]["original_scale_log1p_date_covariate"]=fit(x,baseline)
    report["models"]["original_scale_log1p_date_and_location_covariate"]=fit(x,withsite)
    report["models"]["raw_counts_NB_alpha_1_date_covariate"]=fit(x,count,sm.families.NegativeBinomial(alpha=1))
    w=x.copy()
    high=float(w.grains.quantile(.95))
    w["y"]=np.log1p(np.minimum(w.grains,high))
    report["source_95pct_grain_winsorization_threshold"]=scalar(high)
    report["models"]["winsor95_log1p_date_covariate"]=fit(w,baseline)
    report["date_cluster_bootstrap_on_unadjusted_one_df_contrast"]=block_bootstrap(x)
    report["within_day_patch_treatment_permutation"]=permutation_within_date_patch(x)
    report["conclusion_boundaries"]=[
      "One-df post-hoc patch effect can disagree with original 2-df omnibus ANOVA p~0.06",
      "Patch treated environment was not randomized among all source sites",
      "Focal stigma grains are not donor-identified self grains, compatible pollen tubes or offspring",
      "Different specification sensitivity not independent source replications",
      "Avoid declaring evolutionary compensation or a new universal transition threshold from one study"]
    report["status"]="success_secondary_specification_robustness_audit"
    (OUT/"clarkia_study5_pollen_supply_patch_interaction_sensitivity.json").write_text(
      json.dumps(scalar(report),ensure_ascii=False,indent=2)+"\n")
    print("CLARKIA_STUDY5_SENSITIVITY",json.dumps(scalar(report),ensure_ascii=False)[:19000],flush=True)
if __name__=="__main__":main()
