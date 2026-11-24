#!/usr/bin/env python3
"""v0.54: Within-mother same-father sibling clustering beyond distance,
male flower display, and observed sire productivity in Daphne jezoensis.

Published Shibata & Kudo (2023), doi:10.1002/ajb2.16239, data from
Zenodo 8253921, pinned SHA by v28. This is the SAME source cohort as
v30-v38/v53, not independent data or novel multiple-paternity discovery.

Unit: SSR-assigned offspring to strictly H-coded fathers with positive
natural flower counts (>=.95 assignment), conditional on mother families.
- Null A: independent allocations within each sampled mother using
           distance+flower weighted father probabilities, fixing mother n.
- Null B: shuffle existing paternal labels over offspring slots, fixing
           exact mother n and global father totals (no spatial weights).
- Null C: fixed mother n and global sire totals; Metropolis edge swaps
           weighted by exp(-mother-father distance/lambda).
  The flower-number part of the earlier donor weight cancels exactly
  when per-father offspring totals are held fixed.
Use family-weighted full-sibling pair fraction, not a fictitious number
of independent parents. Source IDs, coordinates, assignments, individual
family statistics NEVER written to output or logs.

Monte Carlo nulls include the *sampled* parentage data only. They do
not identify a single pollen-dispersal event, a visitor identity, complete
male fitness, pollen competition, or the cause of any excess clustering.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.special import softmax
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data, ZIP_SHA
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import build, estimate

RNG_SEED=20261124
LAMBDA_PRIMARY=10.0
INDEPENDENT_NULL_DRAWS=10000
FIXED_MARGIN_DRAWS=10000
MCMC_DRAWS=1500
MCMC_CHAINS=2
MCMC_BURNIN=18000
MCMC_THIN=260
LAMBDA_SENS=(5.,10.,20.)
OUT=Path("analysis_results_daphne_mother_sire_clustering_v54")

def pair_fraction(counts):
    n=counts.sum(axis=1)
    pairs_total=int(np.sum(n*(n-1)//2))
    assert pairs_total>0
    same=int(np.sum(counts*(counts-1)//2))
    return same/pairs_total,pairs_total,same

def rate_samples(samples,obs):
    lo,hi=np.quantile(samples,[.025,.975])
    return {
     "null_mean":float(np.mean(samples)),
     "null_sd":float(np.std(samples,ddof=1)),
     "null_95":[float(lo),float(hi)],
     "one_sided_excess_p":float((1+np.count_nonzero(samples>=obs))/(len(samples)+1)),
     "two_sided_centered_p":float((1+np.count_nonzero(np.abs(samples-np.mean(samples))>=abs(obs-np.mean(samples))))/(len(samples)+1)),
     "observed_minus_null_mean":float(obs-np.mean(samples)),
     "draws":int(len(samples))
    }

def member_slots(counts):
    m,k=counts.shape
    mother=np.repeat(np.arange(m),counts.sum(axis=1))
    father=np.concatenate([np.repeat(np.arange(k),counts[j]) for j in range(m)])
    assert len(mother)==len(father)==int(counts.sum())
    assert np.array_equal(np.bincount(mother,minlength=m),counts.sum(axis=1))
    assert np.array_equal(np.bincount(father,minlength=k),counts.sum(axis=0))
    return mother.astype(np.int32),father.astype(np.int32)

def equal_mother_pair_fraction(counts):
    n=counts.sum(axis=1)
    family_pairs=n*(n-1)//2
    i=family_pairs>0
    assert i.sum()>0
    numerator=np.sum(counts*(counts-1)//2,axis=1)
    return float(np.mean(numerator[i]/family_pairs[i]))

def family_pair_fraction_slots(mother,father,n_moms,n_fathers,total_pairs):
    k=np.bincount(mother*n_fathers+father,minlength=n_moms*n_fathers).reshape(n_moms,n_fathers)
    n=k.sum(axis=1)
    family_pairs=n*(n-1)//2
    same_by_mother=np.sum(k*(k-1)//2,axis=1)
    i=family_pairs>0
    return (float(same_by_mother.sum()/total_pairs),
            float(np.mean(same_by_mother[i]/family_pairs[i])))

def random_margins_null(counts,draws,rng,pairs_total):
    m,f=member_slots(counts)
    vals=np.empty((draws,2))
    for i in range(draws):
        vals[i]=family_pair_fraction_slots(m,rng.permutation(f),*counts.shape,pairs_total)
    return vals

def independent_null(counts,p,rng,pairs_total):
    nrows=counts.sum(axis=1)
    family_pairs=nrows*(nrows-1)//2
    eligible=family_pairs>0
    draws=np.zeros((INDEPENDENT_NULL_DRAWS,2))
    for n,npair,prob in zip(nrows,family_pairs,p):
        if n<2:continue
        k=rng.multinomial(int(n),prob,size=INDEPENDENT_NULL_DRAWS)
        paired=np.sum(k*(k-1)//2,axis=1)
        draws[:,0]+=paired/pairs_total
        draws[:,1]+=paired/npair/int(eligible.sum())
    return draws

def mcmc_fixed_margins(counts,dist,lam,seed,draws=MCMC_DRAWS):
    rng=np.random.default_rng(seed)
    mother,father=member_slots(counts)
    n=len(mother);m,k=counts.shape
    pairs_total=int(np.sum(counts.sum(axis=1)*(counts.sum(axis=1)-1)//2))
    sample=np.empty((draws,2))
    # Start one chain at the observed configuration (potentially highly
    # clustered), one at the uniform shuffled configuration; stationarity
    # and sensitivity to starts can be inspected from chain means.
    chain_id=(seed//10)%2
    if chain_id==1:father=rng.permutation(father)
    attempts=0;accepted=0;new_pair=0
    for step in range(MCMC_BURNIN+draws*MCMC_THIN):
        ix=rng.integers(n,size=2)
        a,b=int(ix[0]),int(ix[1])
        if a!=b:
            ma,mb=int(mother[a]),int(mother[b])
            fa,fb=int(father[a]),int(father[b])
            if ma!=mb and fa!=fb:
                attempts+=1
                dd=(dist[ma,fb]+dist[mb,fa]-dist[ma,fa]-dist[mb,fb])/lam
                if dd<=0 or rng.random()<np.exp(-dd):
                    father[a],father[b]=fb,fa
                    accepted+=1
        if step>=MCMC_BURNIN and (step-MCMC_BURNIN+1)%MCMC_THIN==0:
            sample[new_pair]=family_pair_fraction_slots(mother,father,m,k,pairs_total)
            new_pair+=1
    assert new_pair==draws
    assert np.array_equal(np.bincount(father,minlength=k),counts.sum(axis=0))
    assert np.array_equal(np.bincount(mother,minlength=m),counts.sum(axis=1))
    return sample,{"attempted_valid_label_swaps":attempts,"accepted_swaps":accepted,
                   "conditional_swap_acceptance":float(accepted/max(1,attempts)),
                   "chain_pair_mean":float(sample[:,0].mean()),
                   "chain_pair_95":np.quantile(sample[:,0],[.025,.975]).tolist(),
                   "chain_mother_equal_mean":float(sample[:,1].mean()),
                   "chain_mother_equal_95":np.quantile(sample[:,1],[.025,.975]).tolist(),
                   "chain_start":"observed" if chain_id==0 else "uniform_randomized"}

def analyse_year(q,year):
    counts=q["counts"].astype(np.int32)
    m,k=counts.shape
    assert int(counts.sum())==(180 if year==2018 else 206)
    assert k==(105 if year==2018 else 71)
    obs,pairs_total,pairs_same=pair_fraction(counts)
    equal_obs=equal_mother_pair_fraction(counts)
    mother_sizes=counts.sum(axis=1)
    donor_sizes=counts.sum(axis=0)
    beta=float(estimate(q,LAMBDA_PRIMARY)["beta_hat"])
    prob=softmax(-q["distance"]/LAMBDA_PRIMARY+beta*q["logflowers"][None,:],axis=1)
    rng=np.random.default_rng(RNG_SEED+year)
    nullA=independent_null(counts,prob,rng,pairs_total)
    nullB=random_margins_null(counts,FIXED_MARGIN_DRAWS,rng,pairs_total)
    lam_results={}
    for il,lam in enumerate(LAMBDA_SENS):
        chains=[];meta=[]
        for c in range(MCMC_CHAINS):
            # seed ending 0=>observed; seed ending 1=>uniform
            seed=RNG_SEED+year*100+il*20+c
            # seed parity in tens set so chain_id below truly matches
            seed=(seed//10)*20+c*10
            vals,stat=mcmc_fixed_margins(counts,q["distance"],lam,seed)
            chains.append(vals);meta.append(stat)
        combined=np.concatenate(chains)
        lam_results[str(int(lam))]={
            **rate_samples(combined[:,0],obs),
            "mother_equal_metric":rate_samples(combined[:,1],equal_obs),
            "chains":meta,
            "chain_mean_absolute_difference":float(abs(meta[0]["chain_pair_mean"]-meta[1]["chain_pair_mean"])),
            "chain_mother_equal_mean_absolute_difference":float(abs(meta[0]["chain_mother_equal_mean"]-meta[1]["chain_mother_equal_mean"])),
            "source_data_observed_start_and_random_start_in_same_window":bool(abs(meta[0]["chain_pair_mean"]-meta[1]["chain_pair_mean"])<0.05 and abs(meta[0]["chain_mother_equal_mean"]-meta[1]["chain_mother_equal_mean"])<0.05)
        }
    # F sex mother subgroup source descriptive only; mothers, not offspring,
    # are the relevant biological clusters and the source fraction can be
    # conditioning on maternally sampled fruit set.
    groups={}
    for code in ["F","H","HF"]:
        ix=np.flatnonzero(np.array(q["private_mother_sex_codes"])==code)
        if len(ix)>0:
            sub=counts[ix]
            npairs=int(sum(sub.sum(axis=1)*(sub.sum(axis=1)-1)//2))
            groups[code]={"mothers":len(ix),"offspring_with_strict_H_fathers":int(sub.sum()),
                "informative_sibling_pairs":npairs,
                "observed_same_father_pair_fraction":float(np.sum(sub*(sub-1)//2)/npairs) if npairs else None}
    return {
        "year":year,
        "eligible_strict_H_candidate_fathers":k,
        "sampled_mothers_with_eligible_H_offspring":m,
        "source_assigned_offspring":int(counts.sum()),
        "donor_genets_with_at_least_one_assigned":int((donor_sizes>0).sum()),
        "mothers_with_at_least_two_genotyped":int((mother_sizes>=2).sum()),
        "min_max_offspring_in_one_maternal_family":[int(mother_sizes.min()),int(mother_sizes.max())],
        "maternal_family_size_mean":float(mother_sizes.mean()),
        "within_family_sibling_pairs":pairs_total,
        "same_father_sibling_pairs_observed":pairs_same,
        "same_father_sibling_pair_fraction_observed":obs,
        "mother_equal_same_father_pair_fraction_observed":equal_obs,
        "fitted_flower_coefficient_at_lambda10m":beta,
        "inference":"Only children with SSR-assigned strictly H-coded fathers; conditional on available candidate father set and sampled maternal families",
        "independent_mother_genet_null":rate_samples(nullA[:,0],obs),
        "independent_mother_genet_null_mother_equal":rate_samples(nullA[:,1],equal_obs),
        "exact_global_father_mother_totals_exchangeable_null":rate_samples(nullB[:,0],obs),
        "exact_global_father_mother_totals_exchangeable_null_mother_equal":rate_samples(nullB[:,1],equal_obs),
        "spatially_weighted_exact_margins_nulls":lam_results,
        "source_mother_sex_strata_descriptive":groups,
        "limitations":["No independently marked individual pollinator visitation events",
            "Paternity assigned at >=0.95, father candidates incomplete",
            "Published Daphne fruits are single-seeded; offspring can come from multiple separate fruits of a mother",
            "Source is published paternity diversity experiment, not newly collected field data",
            "Paternity concentrates through maternal mating compatibility or temporal phenology as well as insect carryover",
            "Fixed margin null conditions on realized globally successful sires",
            "Metropolis swaps approximate but do not prove exact mixed equilibrium",
            "Source parentage may pool natural/day-open/night-open treatments on the same maternal genet",
            "Original MasterBayes sire assignment includes a distance prior; the source father labels are not distance-blind",
            "Mother-equal pair metric is restricted to mothers with at least two assigned H-sired offspring",
            "Maternal genotype and paternal source identity not output"]}

def main():
    fs=data()
    y=[analyse_year(build(year,fs),year) for year in (2018,2020)]
    assert y[0]["same_father_sibling_pairs_observed"]<=y[0]["within_family_sibling_pairs"]
    assert y[1]["same_father_sibling_pairs_observed"]<=y[1]["within_family_sibling_pairs"]
    return {
      "status":"original_Daphne_within_mother_father_concentration_fixed_margin_v54",
      "paper_doi":"10.1002/ajb2.16239",
      "source_zenodo":"10.5281/zenodo.8253921",
      "source_zip_sha256":ZIP_SHA,
      "rng_seed":RNG_SEED,
      "yearly":y,
      "inferential_limits":{
        "same_original_2023_publication_as_v28_to_v53":True,
        "parentage_source_probability_095_threshold":True,
        "genet_is_observation_unit_of_father_mother_not_each_child":True,
        "null_conditioned_on_original_sampled_mothers_and_year_specific_sires":True,
        "concentrated_fathers_does_not_establish_repeat_individual_pollinator":True,
        "sex_specific_brood_comparisons_descriptive_no_causal_assignment":True,
        "Daphne_pseudomezereum_not_observed":True,
        "no_original_genet_ID_location_or_parentage_rows_exported":True}}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT/"daphne_parentage_sire_clustering_v54.json")
    args=p.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result=main()
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("DAPHNE_WITHIN_MOTHER_V54",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
