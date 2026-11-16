#!/usr/bin/env python3
"""v0.30 | Does floral display predict *observed* fatherhood, conditional on distance?

Reanalysis of ORIGINAL Shibata & Kudo (2023), AJB doi:10.1002/ajb2.16239,
public author Zenodo record 8253921.

Only directly recorded one-letter H plants with positive same-year *NaturalFlowerNo*
and known adult coordinates are eligible as father candidates. We condition on
offspring with a >=0.95 SSR-assigned father *in the eligible H risk set*; all
excluded offspring and reasons are counted, not filled with zero.

Conditional mother-choice likelihood:
p(sire_i | sampled mother_j, year, lambda, beta) propto
   (number_of_naturally_recorded_H_flowers_i)^beta * exp(-d_ij/lambda)
Predeclared lambda = 5/10/20 m. Beta=0 is distance-only; beta=1 is a
proportional-to-flowers heuristic. Same-year fatherhood is a sampled relative
paternity opportunity, NOT lifetime male fitness, causation, floral mass, or
pollination by any particular insect.

Cluster bootstrap samples maternal families. Leave-one-mother-out log-score
reports out-of-mother predictive advantage over beta=0 to reduce overfitting.
No source individual IDs, raw genotypes or coordinates are logged/artifacted.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import spearmanr
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data as source, stringify

OUT=Path("analysis_results_daphne_flower_to_paternity_v30")
FILENAMES={
 2018:"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
 2020:"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"}
N_BOOT=800
RNG=np.random.default_rng(20261116)
LAMBDAS=(5.,10.,20.)

def val(x):
    x=float(x)
    return round(x,6) if np.isfinite(x) else None
def build(year,files,include_HF=False):
    raw=files[FILENAMES[year]]
    adult=files["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(raw)==(257 if year==2018 else 327)
    assert len(adult)==492
    sexcol=f"Sex_{year}.g"
    flowercol=f"NaturalFlowerNo_{year}.g"
    assert {"GenetID2","East.g","North.g",sexcol,flowercol}<=set(adult)
    ad=adult[["GenetID2","East.g","North.g",sexcol,flowercol]].copy()
    ad["genet"]=stringify(ad.GenetID2)
    assert not ad.genet.duplicated().any()
    ad["x"]=pd.to_numeric(ad["East.g"],errors="coerce")
    ad["y"]=pd.to_numeric(ad["North.g"],errors="coerce")
    ad["flower"]=pd.to_numeric(ad[flowercol],errors="coerce")
    ad["sexcode"]=ad[sexcol].astype(str)
    sex_candidate_codes=["H","HF"] if include_HF else ["H"]
    eligible=ad[(ad.sexcode.isin(sex_candidate_codes))&ad.flower.gt(0)&ad.x.notna()&ad.y.notna()].copy()
    assert len(eligible)==((108 if year==2018 else 79) if include_HF else (105 if year==2018 else 71))
    eligible=eligible.reset_index(drop=True)
    assert eligible.genet.nunique()==len(eligible)
    raw=raw.copy()
    raw["Father_key"]=stringify(raw.Father)
    raw["Mother_key"]=stringify(raw.Mother)
    raw["Probability"]=pd.to_numeric(raw.Probability,errors="coerce")
    assigned=raw[raw.Father.notna()&raw.Probability.ge(.95)].copy()
    assert len(assigned)==(209 if year==2018 else 261)
    mother_ad=ad[["genet","x","y","sexcode"]].rename(columns={
      "genet":"Mother_key","x":"mx","y":"my","sexcode":"mother_sex"})
    father_ad=ad[["genet","sexcode","flower","x","y"]].rename(columns={
      "genet":"Father_key","sexcode":"father_sex","flower":"father_flowers",
      "x":"fx","y":"fy"})
    j=assigned.merge(mother_ad,on="Mother_key",how="left",validate="m:1")
    j=j.merge(father_ad,on="Father_key",how="left",validate="m:1")
    assert len(j)==len(assigned)
    assert j[["mx","my","fx","fy","father_sex"]].notna().all().all()
    father_H=j[j.father_sex.isin(sex_candidate_codes)].copy()
    assert len(father_H)==((183 if year==2018 else 235) if include_HF else (180 if year==2018 else 227))
    father_H["eligible"]=father_H.Father_key.isin(eligible.genet)
    rows=father_H[father_H.eligible].copy()
    # Exact unit: each row is confidently assigned offspring; repeats within
    # a maternal family are correlated and bootstrapped as one unit.
    assert rows["father_flowers"].gt(0).all()
    assert rows["mx"].notna().all()
    father_to_index={g:i for i,g in enumerate(eligible.genet)}
    rows["sire_index"]=rows.Father_key.map(father_to_index)
    assert rows.sire_index.notna().all()
    # Mother gene IDs are labels only for grouping; no transformation to space.
    mum=rows.groupby("Mother_key",sort=True)
    moms=[(str(k),g) for k,g in mum]
    assert len(moms)>=12, "too few independent recipient families"
    mother_xy=np.array([[g.mx.iloc[0],g.my.iloc[0]] for _,g in moms],dtype=float)
    assert np.all(np.isfinite(mother_xy))
    # Guard mother coordinate consistency within family
    assert all(g.mx.nunique()==1 and g.my.nunique()==1 for _,g in moms)
    xy=eligible[["x","y"]].to_numpy(dtype=float)
    distance=np.hypot(mother_xy[:,None,0]-xy[None,:,0],
                      mother_xy[:,None,1]-xy[None,:,1])
    counts=np.zeros((len(moms),len(eligible)),dtype=float)
    for k,(_,g) in enumerate(moms):
        ix=g.sire_index.to_numpy(dtype=int)
        np.add.at(counts[k],ix,1)
    assert int(counts.sum())==len(rows)
    assert np.all(counts.sum(axis=1)>=1)
    logflower=np.log(eligible.flower.to_numpy(float))
    qc={
      "year":year,
      "candidate_sex_codes":sex_candidate_codes,
      "original_seed_rows":len(raw),
      "original_95percent_assigned_rows":len(assigned),
      "assigned_to_H_sex_code":len(father_H),
      "assigned_H_excluded_due_to_no_positive_flower_or_coordinate":len(father_H)-len(rows),
      "assigned_H_eligible_for_conditional_model":len(rows),
      "eligible_source_H_flowering_genets":len(eligible),
      "eligible_H_fathers_with_at_least_one_sampled_offspring":int((counts.sum(axis=0)>0).sum()),
      "eligible_H_fathers_with_no_sampled_offspring":int((counts.sum(axis=0)==0).sum()),
      "independent_sampled_mothers_with_eligible_H_offspring":len(moms),
      "sampled_mother_sex_code_counts":{str(k):int(v) for k,v in rows.groupby("Mother_key").mother_sex.first().value_counts().items()},
      "same_sex_year_recorded_flower_count_available":True,
      "candidate_father_genet_zero_assignment_is_only_sampled_observation":True,
      "n_father_eligible_risk_set_matching_genet_source":int(len(rows)),
      "mother_offspring_range":[int(counts.sum(axis=1).min()),int(counts.sum(axis=1).max())],
    }
    # Descriptive raw count association not controlling mother opportunity.
    sire_count=counts.sum(axis=0)
    qc["naive_H_flower_count_vs_sire_offspring_spearman"]=val(spearmanr(
       eligible.flower.to_numpy(float),sire_count).statistic)
    most=np.argsort(eligible.flower.to_numpy(float))
    top=most[-max(1,int(np.ceil(len(eligible)*.1))):]
    qc["top10pct_by_flower_genet_count"]=len(top)
    qc["top10pct_H_genets_share_of_candidate_flower_counts"]=val(
      eligible.flower.to_numpy(float)[top].sum()/eligible.flower.sum())
    qc["top10pct_H_genets_share_of_assigned_H_offspring"]=val(sire_count[top].sum()/sire_count.sum())
    return {"year":year,"counts":counts,"distance":distance,
            "logflowers":logflower,"qc":qc,
            "private_mother_xy":mother_xy,
            "private_candidate_father_xy":xy,
            "private_mother_sex_codes":[str(g.mother_sex.iloc[0]) for _,g in moms]}

def estimate(q,lam,weights=None,extra=False):
    counts=q["counts"]
    n=counts.sum(axis=1)
    distance=q["distance"]
    flowers=q["logflowers"]
    w=np.ones(len(counts)) if weights is None else weights
    assert len(w)==len(counts)
    assert w.sum()>0
    baseline=-distance/lam
    def loglik(b):
        logits=baseline+b*flowers[None,:]
        return float(np.sum(w[:,None]*counts*logits)-
                     np.sum(w*n*logsumexp(logits,axis=1)))
    opt=minimize_scalar(lambda b:-loglik(b),bounds=(-2,3.5),
                        method="bounded",options={"xatol":0.004})
    if not opt.success:raise RuntimeError("one-parameter source sire fit failed")
    beta=float(opt.x)
    result={"beta_hat":val(beta),"fit_at_beta_boundary":bool(beta<-1.98 or beta>3.48),
      "per_child_in_sample_gain_over_distance_only":val((loglik(beta)-loglik(0))/(w@n))}
    if extra:
        result.update({"distance_only_log_score_per_child":val(loglik(0)/(w@n)),
          "flower_proportional_beta1_gain_per_child":val((loglik(1)-loglik(0))/(w@n)),
          "observed_H_child_rows":int(counts.sum()),
          "possible_sire_genets":len(flowers),
          "recipient_mothers":len(counts)})
        # LOO predictive log score; retain no per-family values.
        deltas=[]
        for i in range(len(counts)):
            wt=np.ones(len(counts));wt[i]=0
            e=estimate(q,lam,weights=wt,extra=False)
            held_idx=i
            U=baseline[held_idx]+e["beta_hat"]*flowers
            logp=U-logsumexp(U)
            U0=baseline[held_idx]
            logp0=U0-logsumexp(U0)
            delta=float(np.dot(counts[held_idx],logp-logp0))
            deltas.append(delta)
        result["mother_leave_one_out_gain_per_child"]=val(sum(deltas)/counts.sum())
        result["n_holdout_mother_families"]=len(deltas)
        result["loo_families_with_positive_log_score_gain"]=int(sum(v>0 for v in deltas))
        result["loo_families_with_negative_log_score_gain"]=int(sum(v<0 for v in deltas))
        result["loo_total_score_gain_over_distance_only"]=val(sum(deltas))
    return result

def boot(q,lam):
    b=[]
    n=len(q["counts"])
    for i in range(N_BOOT):
        ids=RNG.integers(0,n,n)
        weights=np.bincount(ids,minlength=n).astype(float)
        beta=estimate(q,lam,weights=weights)["beta_hat"]
        if beta is not None:b.append(beta)
    assert len(b)>=N_BOOT*.95
    return {"n_family_resamples":len(b),
        "beta_mother_boot_ci95":[val(x) for x in np.quantile(b,[.025,.975])],
        "fraction_boot_beta_positive":val(np.mean(np.array(b)>0)),
        "fraction_boot_beta_less_than_one":val(np.mean(np.array(b)<1))}

def main():
    OUT.mkdir(exist_ok=True)
    files=source() # pinned SHA256 original author source via v28
    yearly=[]
    for year in (2018,2020):
        q=build(year,files)
        results={"year":year,"quality_control":q["qc"],"kernel_sensitivity":[]}
        for lam in LAMBDAS:
            fit=estimate(q,lam,extra=True)
            fit.update(boot(q,lam))
            fit["distance_kernel_decay_metres"]=int(lam)
            fit["mother_family_cluster_bootstrap_not_independent_seeds"]=True
            results["kernel_sensitivity"].append(fit)
        yearly.append(results)
    out={"status":"original_Daphne_H_flower_to_observed_paternity_conditional_v30",
       "original_paper_doi":"10.1002/ajb2.16239",
       "original_public_data":"10.5281/zenodo.8253921",
       "source_integrity_file":"scripts/daphne_jezoensis_2018_2020_HF_fathers_v28.py",
       "years":yearly,
       "strict_assumptions":{
         "H_candidates_only_verified_single_code_and_natural_flower_abundance":True,
         "offspring_with_missing_same_year_father_flower_counts_excluded_not_zero_filled":True,
         "assigned_offspring_conditioned_on_threshold_095_and_eligible_H_fathers":True,
         "distance_decay_m_preselected_5_10_20_sensitivity_only":True,
         "maternal_family_bootstrap_not_seed_level_independence":True,
         "male_flower_count_not_dry_mass_biomass":True,
         "sire_counts_are_sampled_not_population_wide_lifetime_siring":True,
         "mother_holdout_predicts_sampled_sire_choice_not_offspring_survival":True,
         "paternal_flower_number_phenology_at_pollination_unmeasured":True,
         "source_one_population_two_years_not_independent_species_replication":True,
         "no_focal_D_pseudomezereum_data":True},
       }
    path=OUT/"daphne_2018_2020_H_flower_numbers_conditional_paternity_v30.json"
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_MALE_DISPLAY_PATERNITY_V30",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
