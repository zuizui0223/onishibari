#!/usr/bin/env python3
"""v0.38 | Retrospective cross-year father-choice model transfer & sire influence.

Source: author-original Shibata & Kudo 2023, Zenodo 8253921, source-sha
verified by v28. Same population, distinct offspring cohorts in 2018/2020.

Train distance-conditioned flowering coefficient beta on ALL mother families
from one year and apply without re-fitting to offspring of other year.
Report log-probability gain per confidently assigned offspring vs beta0.
Recipient-mother bootstrap conditions on the training-year beta and cannot
quantify training-year, sire-family or general geographical uncertainty.

Separately refit after removing one father candidate and all offspring
attributed to him, conditioning on the reduced risk set. This checks whether
the entire positive association is driven by ONE disproportionate sire,
but is not an independent father bootstrap or test of causal floral selection.

No raw identities, genotype data, coordinates, father-level results,
or individual paternity rows are output.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import build,estimate
from daphne_jezoensis_2018_2020_paternal_adversarial_v32 import reduced

OUT=Path("analysis_results_daphne_crossyear_paternity_v38")
RNG=np.random.default_rng(20261117)
LAMBDA=10.0
N_MOTHER_BOOT=10000

def z(v):
    return round(float(v),6) if np.isfinite(v) else None

def evaluation(q,beta,bootstrap=True):
    c=q["counts"]
    lf=q["logflowers"]
    base=-q["distance"]/LAMBDA
    a=base+beta*lf[None,:]
    alt=np.sum(c*(a-logsumexp(a,axis=1)[:,None]),axis=1)
    nul=np.sum(c*(base-logsumexp(base,axis=1)[:,None]),axis=1)
    bymother=alt-nul
    n=c.sum(axis=1)
    result={"n_assigned_offspring":int(n.sum()),
       "n_heldout_maternal_families":len(n),
       "gain_over_distance_only_nats_per_child":z(bymother.sum()/n.sum()),
       "families_with_positive_log_score_gain":int((bymother>0).sum()),
       "families_with_negative_log_score_gain":int((bymother<0).sum()),
       "positive_gain_total":bool(bymother.sum()>0),
       "not_an_independent_site_or_species_holdout":True}
    if bootstrap:
        ix=RNG.integers(0,len(n),size=(N_MOTHER_BOOT,len(n)))
        g=bymother[ix].sum(axis=1)/n[ix].sum(axis=1)
        result["fixed_train_beta_test_mother_boot_95%"]=[
          z(v) for v in np.quantile(g,[.025,.975])]
    return result

def donor_influence(q):
    v=q["counts"].sum(axis=0)
    active=np.where(v>0)[0]
    allcoef=[]
    losses=[]
    for idx in active:
        keep=np.ones(len(v),dtype=bool)
        keep[idx]=False
        r=reduced(q,sire_keep=keep)
        assert r["counts"].sum()<q["counts"].sum()
        fit=estimate(r,LAMBDA)
        allcoef.append(float(fit["beta_hat"]))
        losses.append(int(v[idx]))
    allcoef=np.asarray(allcoef,float)
    assert len(allcoef)==len(active)
    return {
       "n_influential_source_sire_genets_tested":len(active),
       "n_sire_genets_with_zero_observed_offspring_not_excluded":int((v==0).sum()),
       "max_assigned_offspring_from_one_source_sire":int(max(losses)),
       "n_different_sire_exclusion_estimates_positive":int((allcoef>0).sum()),
       "lowest_beta_after_dropping_one_sire":z(min(allcoef)),
       "highest_beta_after_dropping_one_sire":z(max(allcoef)),
       "median_beta_after_dropping_one_sire":z(np.median(allcoef)),
       "lower_beta_exclusion_vs_main":z(min(allcoef)-float(estimate(q,LAMBDA)["beta_hat"])),
       "not_independent_sire_resampling":True,
       "sire_exclusion_changes_conditioned_risk_set":True
    }

def do_group(fs,includeHF):
    samples={y:build(y,fs,include_HF=includeHF) for y in (2018,2020)}
    coeff={y:float(estimate(q,LAMBDA)["beta_hat"]) for y,q in samples.items()}
    transfers=[]
    for train,test in ((2018,2020),(2020,2018)):
        case=evaluation(samples[test],coeff[train])
        case.update({
          "train_year":train,"test_year":test,
          "beta_trained_other_year":z(coeff[train]),
          "in_test_year_best_fit_beta_for_context_only":z(coeff[test]),
          "train_only_model_applied_without_test_year_refit":True,
        })
        transfers.append(case)
    influence=[]
    for year in (2018,2020):
        influence.append({"year":year,
          "n_candidate_genets":len(samples[year]["logflowers"]),
          "n_ssr_assigned_sire_children":int(samples[year]["counts"].sum()),
          "main_beta":z(coeff[year]),
          **donor_influence(samples[year])})
    return {"donor_sex_codes":["H","HF"] if includeHF else ["H"],
      "per_year_betas":[{"year":y,"beta":z(coeff[y])} for y in (2018,2020)],
      "cross_year_transfers":transfers,
      "one_observed_sire_genet_leaveout_sensitivity":influence}

def main():
    OUT.mkdir(exist_ok=True)
    fs=data()
    rows=[do_group(fs,hf) for hf in (False,True)]
    output={
      "status":"v38_Daphne_2018_2020_crossyear_father_choice_and_sire_influence",
      "original_2023_paper_doi":"10.1002/ajb2.16239",
      "source_doi":"10.5281/zenodo.8253921",
      "distance_lambda_m":int(LAMBDA),
      "groups":rows,
      "limits":{
        "cross_year_transfer_is_retrospective_not_pre_registered_forecast":True,
        "same_local_population_both_years_not_external_biogeographic_replication":True,
        "observed_sire_drop_is_conditioned_on_sire_having_assigned_offspring":True,
        "candidate_father_risk_set_limited_to_exact_H_or_HF_positive_flowers":True,
        "seed_paternity_sample_not_full_population_fitness":True,
        "flower_number_not_true_pollen_export_or_construction_biomass":True,
        "individual_plants_may_recur_between_years":True,
        "unmeasured_sire_spatial_resource_and_compatibility_confounding":True,
        "no_original_father_genet_identifiers_or_coords_emitted":True,
        "not_focal_Daphne_pseudomezereum":True}}
    path=OUT/"daphne_fatherhood_year_transfer_and_sire_exclusion_v38.json"
    path.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_YEAR_TRANSFER_SIRE_INFLUENCE_V38",json.dumps(output,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
