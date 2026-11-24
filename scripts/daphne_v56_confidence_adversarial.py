#!/usr/bin/env python3
"""v0.56: exact-source confidence filters and adversarial paternal label sensitivity.

The user-facing claim is about already inferred fathers, NOT a re-estimated
genetic paternity model. The published Zenodo 8253921 archive holds adult
16-SSR alleles but NOT offspring allele genotypes or alternative parentage
posteriors. Thus distance-independent genetic reassignment is unavailable.

1. Keep natural pollination ONLY, original strict H flowering candidate sires,
   >=0.95, >=0.975, >=0.99, >=0.995 MasterBayes posterior calls.
2. Recompute mother-equal within-mother same-father pairing, all exact
   assigned father totals, maternal sample sizes, and the 5m-distance
   Metropolis-swap null (most generous distance alternative of v55).
3. Compute a maximally hostile, genetically UNCONSTRAINED number of
   father-label changes required to reduce observed mother-equal pair rate
   below the fixed-margin spatial null's 97.5% quantile. This is a
   deterministic counterfactual fragility bound, NOT evidence that any
   particular original assignment was wrong.
4. Calculate sum(1-p) only as expected number wrong under the SAME
   source posterior model. It cannot bound errors from mis-specified priors.
All original IDs, alleles, parental edges, coordinates stay private in memory.
Outputs consist of hashes, aggregate counts, scenario estimates.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,stringify,ZIP_SHA
from daphne_jezoensis_treatment_stratified_paternity_v55 import sample_file,setup
from daphne_jezoensis_parentage_clustering_v54 import (
    pair_fraction,equal_mother_pair_fraction,mcmc_fixed_margins,rate_samples)

OUT=Path("analysis_results_daphne_confidence_prior_v56")
SEED=20261124
THRESHOLDS=[.95,.975,.99,.995]
SPATIAL_LAMBDA=5.
N_CHAINS=2
N_SAVED=500

def original_natural_cohort(year,files,meta,setup_baseline,cutoff):
    raw=files[f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv"].copy()
    source=raw.merge(meta,on="DNAID",how="left",validate="one_to_one")
    assert len(source)==(257 if year==2018 else 327)
    assert source.Info.notna().all()
    prob=pd.to_numeric(source.Probability,errors="coerce")
    k=source.loc[source.Info.eq("Natural")&source.Father.notna()&prob.ge(cutoff)].copy()
    k["Father_key"]=stringify(k.Father)
    k["Mother_key"]=stringify(k.Mother)
    k=k.loc[k.Father_key.isin(set(setup_baseline["sireindex"]))].copy()
    assert len(k)>40 and k.Mother_key.nunique()>=12
    # When cutoff is the original .95, EXACTLY recover baseline v55
    if cutoff==.95:
        reference=setup_baseline["good"]
        ref=reference.loc[reference.Info=="Natural"]
        assert len(k)==len(ref)==(110 if year==2018 else 131)
        assert set(k.DNAID)==set(ref.DNAID)
    moms=sorted(k.Mother_key.unique())
    fatherix=setup_baseline["sireindex"]
    C=np.zeros((len(moms),len(fatherix)),np.int32)
    adultXY=setup_baseline["adultXY"]
    for j,genet in enumerate(moms):
        for father in k.loc[k.Mother_key==genet,"Father_key"]:
            C[j,fatherix[father]]+=1
    xy=np.array([[adultXY.loc[id0,"x"],adultXY.loc[id0,"y"]] for id0 in moms],float)
    donor=setup_baseline["sire_xy"]
    D=np.hypot(xy[:,None,0]-donor[None,:,0],xy[:,None,1]-donor[None,:,1])
    assert int(C.sum())==len(k)
    sprob=pd.to_numeric(k.Probability).to_numpy(float)
    assert sprob.min()>=cutoff and (sprob<=1.0).all()
    return C,D,sprob

def adversarial_label_errors(C,upper_null):
    """Minimum number of arbitrary father-label changes to hit upper_null.
    Each replacement creates a uniquely new donor for that mother; does NOT
    preserve father totals nor require genotype Mendelian compatibility.
    Because k-1,k-2,... gains are decreasing in each sibship, greedy
    selection of highest per-step gain is optimal for this relaxed bound.
    """
    current=C.copy()
    n=current.sum(axis=1)
    pairs=n*(n-1)//2
    informative=pairs>0
    M=int(informative.sum())
    assert M>0
    original=equal_mother_pair_fraction(current)
    count=0
    while equal_mother_pair_fraction(current)>upper_null:
        slopes=np.zeros_like(current,dtype=float)
        slopes[informative]=np.maximum(0,current[informative]-1)/pairs[informative,None]/M
        if slopes.max()<=0:
            return {"original_mother_equal":original,"upper_null":upper_null,
                    "minimum_reassignment_for_upper_null":None,
                    "max_possible_observed_decrease":original-equal_mother_pair_fraction(current)}
        idx=np.unravel_index(np.argmax(slopes),slopes.shape)
        j,f=map(int,idx)
        vacant=np.flatnonzero(current[j]==0)
        assert len(vacant)>0, "relaxed reassignment should have unused candidate fathers"
        current[j,f]-=1
        current[j,int(vacant[0])]+=1
        count+=1
        assert count<=int(C.sum())
    assert np.array_equal(current.sum(axis=1),n)
    frac=equal_mother_pair_fraction(current)
    return {"original_mother_equal":original,"upper_null":float(upper_null),
            "minimum_reassignment_for_upper_null":int(count),
            "minimum_fraction_of_assigned_natural_seeds":float(count/C.sum()),
            "new_mother_equal_after_relabel":float(frac),
            "note":"Any reassigned donor genotype/location is unconstrained; source father total NOT retained; this is a deliberately liberal adversarial bound, not a corrected pedigree"}

def estimate(C,D,year,cutoff_ix):
    observed=equal_mother_pair_fraction(C)
    vals=[];chains=[]
    for c in range(N_CHAINS):
        # Original and randomized father-marginal-preserving starts.
        seed=(SEED+year*10+cutoff_ix*100)//10*20+c*10
        x,diagnostic=mcmc_fixed_margins(C,D,SPATIAL_LAMBDA,seed,draws=N_SAVED)
        vals.append(x[:,1]);chains.append(diagnostic)
    samples=np.concatenate(vals)
    return {"mother_equal_sibling_observed":float(observed),
            "null_fixed_sire_totals_spatial_lambda5":rate_samples(samples,observed),
            "chain_mean_difference":float(abs(vals[0].mean()-vals[1].mean())),
            "chain_acceptance":[float(x["conditional_swap_acceptance"]) for x in chains],
            "chain_starts":[x["chain_start"] for x in chains]}

def main():
    fs=data()
    meta,meta_n=sample_file()
    assert meta_n==646
    out=[]
    for year in [2018,2020]:
        baseline=setup(year,fs,meta)
        group=[]
        last_n=10**8
        for i,t in enumerate(THRESHOLDS):
            C,D,pr=original_natural_cohort(year,fs,meta,baseline,t)
            assert int(C.sum())<=last_n
            last_n=int(C.sum())
            q=estimate(C,D,year,i)
            upper=q["null_fixed_sire_totals_spatial_lambda5"]["null_95"][1]
            grp={"threshold":float(t),"n_confident_strictH_natural_offspring":int(C.sum()),
                "n_mothers":len(C),"mothers_with_ge2":int((C.sum(axis=1)>=2).sum()),
                "confidence_min":float(pr.min()),
                "confidence_mean":float(pr.mean()),
                "sum_source_posterior_expected_wrong_sire_calls":float(np.sum(1-pr)),
                "posterior_model_based_probability":True,
                "same_father_pairs_pooled":float(pair_fraction(C)[0]),
                **q,
                "adversarial_paternity_relabelling_bound":adversarial_label_errors(C,upper)}
            group.append(grp)
        out.append({"year":year,"strict_H_eligible_father_genets":len(baseline["sirekeys"]),
                    "natural_only_high_confidence_thresholds":group})
    return {
        "status":"v56_Daphne_source_paternity_confidence_adversarial_bounds",
        "source_zenodo":"10.5281/zenodo.8253921","source_zip_sha256":ZIP_SHA,
        "source_paper":"10.1002/ajb2.16239",
        "source_sample_646_has_genotypes":False,
        "source_adult_16_SSR_alleles_present":True,
        "source_offspring_16_SSR_alleles_present":False,
        "has_independent_original_parentage_inference_without_spatial_prior":False,
        "method":"Natural-only, strict H candidate sires; source MasterBayes confidence threshold sweep with exact maternal sizes/sire totals and spatial 5m margin-preserving swaps; relaxed adversarial labels",
        "source_archive_alternative_paternity_sets_available":False,
        "draws_per_chain":N_SAVED,"chains":N_CHAINS,"rng_seed":SEED,
        "by_year":out,
        "limitations":{
          "paternity_posterior_from_same_source_distance_prior":True,
          "source_parentage_genetics_cannot_be_recomputed_with_no_offspring_allele_data":True,
          "MCMC_draws_correlated_not_independent_genet_samples":True,
          "original_mate_genotype_compatibility_unavailable_for_adversarial_changes":True,
          "source_posterior_expected_wrong_under_model_not_bound_under_misspecification":True,
          "natural_pollination_only_not_independent_site_replication":True,
          "source_distance_kernel_in_null_not_distance_blind_genetic_assignment":True,
          "no_source_individual_identifiers_or_coordinates_in_output":True,
          "no_focal_D_pseudomezereum_reproductive_measurement":True}}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT/"daphne_source_confidence_sensitivity_v56.json")
    a=p.parse_args()
    a.output.parent.mkdir(exist_ok=True,parents=True)
    result=main()
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("DAPHNE_CONFIDENCE_V56",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
