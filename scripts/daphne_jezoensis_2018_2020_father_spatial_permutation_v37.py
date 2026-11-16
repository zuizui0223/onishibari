#!/usr/bin/env python3
"""v0.37 | Father-site constrained flower-count permutation nulls.

Purpose: v0.30-33 showed male flower number predicts maternal-conditioned
SSR father assignment when 10m distance is used. But pollen-donor size may
be spatially sorted. Here the observed genotype-specific flower counts are
shuffled among male candidates ONLY inside predeclared 10/20/40m father
coordinate grids (two origins each), with parentage, mother/father coords,
sex risk set, offspring count and distance kernel all fixed. A global
unconstrained shuffle is included as an easier, weak negative control.

Data: original author source of Shibata & Kudo (2023), AJB doi 10.1002/ajb2.16239,
Zenodo 10.5281/zenodo.8253921. The source is SHA256 pinned by the v28 source
reader, with exact 2018/2020 SSR and adult data source assertions in v30.

This is an *observational conditional-label diagnostic*, not a randomized
experiment: donor genets are not exchangeable in reality and ecological,
genetic, resource and compatibility differences can still confound. It cannot
establish selection on flower size, pollen amount, pistil cost or Onishibari
fitness, and it is not an independent published species replication.

No original father or mother IDs, individual positions, source genotypes,
individual assignment rows or random permutations are emitted.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp,softmax
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import (
    build,estimate
)
OUT=Path("analysis_results_daphne_spatial_donor_permutation_v37")
N_PERM=1000
LAMBDA=10.
RNG=np.random.default_rng(20261117)
GRIDS=(10.,20.,40.)
SHIFTS=(0.0,0.5)
BETA_LOW=-2.0
BETA_HIGH=3.5

def z(x):
    v=float(x)
    return round(v,6) if np.isfinite(v) else None

def score_model(q, logflower):
    """Maximize maternal conditional multinomial log likelihood, one beta."""
    d=q["distance"]
    count=q["counts"]
    n=count.sum(axis=1)
    assert len(logflower)==count.shape[1]
    base=-d/LAMBDA
    obslog=float(np.sum(count*logflower[None,:]))
    def grad(beta):
        probs=softmax(base+beta*logflower[None,:],axis=1)
        return obslog-float(np.sum(n*(probs@logflower)))
    gl=grad(BETA_LOW)
    gh=grad(BETA_HIGH)
    if gl<=0:beta=BETA_LOW
    elif gh>=0:beta=BETA_HIGH
    else:beta=float(brentq(grad,BETA_LOW,BETA_HIGH,xtol=1e-6))
    def ll(b):
        return float(np.sum(count*(base+b*logflower[None,:]))-
                     np.dot(n,logsumexp(base+b*logflower[None,:],axis=1)))
    gain=(ll(beta)-ll(0))/float(count.sum())
    return beta,gain

def split_grid(xy, width, shift_fraction):
    assert xy.ndim==2 and xy.shape[1]==2
    # The father coordinates are temporary within this analysis.
    labels=np.floor(xy/width+shift_fraction).astype(int)
    _,inv=np.unique(labels,axis=0,return_inverse=True)
    groups=[np.where(inv==i)[0] for i in np.unique(inv)]
    shuffle=[x for x in groups if len(x)>=2]
    sizes=[len(x) for x in groups]
    return shuffle,{"n_grid_cells_with_candidate_H":len(groups),
       "n_cells_eligible_for_label_shuffle":len(shuffle),
       "n_candidate_H_in_multigenet_cells":int(sum(len(x) for x in shuffle)),
       "fraction_candidate_H_shufflable":z(sum(len(x) for x in shuffle)/len(xy)),
       "max_candidate_genets_in_one_cell":int(max(sizes)),
       "n_single_genet_cells":int(sum(x==1 for x in sizes))}

def spatial_perm(q, observed, label, groups, extra):
    betas=np.empty(N_PERM)
    gains=np.empty(N_PERM)
    lf=q["logflowers"]
    n_reorders=0
    for i in range(N_PERM):
        shuffled=lf.copy()
        for idx in groups:
            shuffled[idx]=lf[RNG.permutation(idx)]
        if not np.array_equal(shuffled,lf):n_reorders+=1
        beta,gain=score_model(q,shuffled)
        betas[i]=beta;gains[i]=gain
    obeta,ogain=observed
    # These tails are conditional-on-geometry diagnostics; not formal design p.
    p_gain=(1+int(np.sum(gains>=ogain-1e-9)))/(N_PERM+1)
    p_beta=(1+int(np.sum(betas>=obeta-1e-9)))/(N_PERM+1)
    return {"permutation_scenario":label,**extra,
        "n_permutations":N_PERM,
        "permutations_changing_a_flower_label":int(n_reorders),
        "observed_conditional_beta":z(obeta),
        "observed_per_offspring_log_likelihood_gain_over_beta0":z(ogain),
        "null_beta_median":z(np.median(betas)),
        "null_beta_ci95":[z(v) for v in np.quantile(betas,[.025,.975])],
        "null_ll_gain_95th_quantile_per_offspring":z(np.quantile(gains,.95)),
        "null_ll_gain_99th_quantile_per_offspring":z(np.quantile(gains,.99)),
        "tail_fraction_gain_at_least_observed_conditional":z(p_gain),
        "tail_fraction_beta_at_least_observed_conditional":z(p_beta),
        "not_exact_randomization_test":True}

def run_year(year,files,with_hf):
    q=build(year,files,include_HF=with_hf)
    observed=score_model(q,q["logflowers"])
    old=estimate(q,LAMBDA,extra=True)
    assert abs(observed[0]-old["beta_hat"])<0.003, "Independent optimizer must recover v30/v33 beta"
    assert abs(observed[1]-old["per_child_in_sample_gain_over_distance_only"])<2e-5
    xy=q["private_candidate_father_xy"]
    n=len(xy)
    designs=[]
    designs.append(spatial_perm(q,observed,"global_shuffle_no_spatial_constraint",
        [np.arange(n)],{"n_candidate_H_in_multigenet_cells":n,
          "fraction_candidate_H_shufflable":1.0,
          "n_grid_cells_with_candidate_H":1,"n_cells_eligible_for_label_shuffle":1,
          "max_candidate_genets_in_one_cell":n,"n_single_genet_cells":0}))
    for m in GRIDS:
        for shift in SHIFTS:
            groups,info=split_grid(xy,m,shift)
            designs.append(spatial_perm(q,observed,
               f"donor_genet_within_grid_{int(m)}m_origin_{shift}",
               groups,{"father_grid_width_m":int(m),
                   "father_grid_origin_shift_fraction":shift,**info}))
    return {"year":year,"donor_sex_codes":["H","HF"] if with_hf else ["H"],
        "n_candidate_genets":n,
        "n_ssr_assigned_children_in_conditioned_father_risk_set":int(q["counts"].sum()),
        "n_recipient_mothers":len(q["counts"]),
        "lambda_metres":int(LAMBDA),
        "beta_observed":z(observed[0]),
        "source_original_mother_leave_one_out_gain_per_child":old["mother_leave_one_out_gain_per_child"],
        "observed_ll_gain_per_child":z(observed[1]),
        "permutation_results":designs,
        "source_eligible_and_sampled_conditioning_unaltered":True}

def main():
    OUT.mkdir(exist_ok=True)
    fs=data()
    allrows=[]
    for year in (2018,2020):
        # main strict H; HF serves as sensitivity, not unambiguously prespecified
        for hf in (False,True):
            allrows.append(run_year(year,fs,hf))
    res={"status":"source_verified_v37_conditional_donor_spatial_label_permutations",
       "article_doi":"10.1002/ajb2.16239",
       "source_doi":"10.5281/zenodo.8253921",
       "years_and_sex_definitions":allrows,
       "limits":{
           "null_shuffles_flower_size_not_paternity_distance_or_fruit":True,
           "spatial_candidate_label_permutations_not_true_random_assignment":True,
           "block_width_sensitive_predeclared_10_20_40m":True,
           "source_genotype_and_flower_size_may_share_unmeasured_environment":True,
           "paternity_sampling_and_95percent_assignment_conditioned":True,
           "repeated_offspring_within_sire_and_maternal_families_not_independent":True,
           "father_coord_label_permutation_preserves_source_model_risk_set":True,
           "one_original_research_population_not_independent_replicates":True,
           "no_evidence_for_focal_Daphne_pseudomezereum":True}}
    path=OUT/"daphne_2018_2020_H_flower_father_spatial_label_permutation_v37.json"
    path.write_text(json.dumps(res,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_FATHER_SPATIAL_LABEL_NULL_V37",json.dumps(res,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
