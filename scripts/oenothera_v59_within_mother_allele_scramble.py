#!/usr/bin/env python3
"""v0.59: Is leave-one-SSR-locus prediction child-specific, or just maternal
and within-treatment sibship genetic structure?

The v58 holdout model finds positive mean log predictions for a withheld locus
in Rhodes et al. (2017) Oenothera harringtonii (Dryad p24q3, Zenodo 4942351).

Here we construct a STRONGER NEGATIVE CONTROL: same mother's offspring,
same treatment, same heldout-locus alleles, missing-data patterns, and the same
6-locus genotype-only potential-father weights are retained. The target
child's withheld 7th-locus genotype is reassigned uniformly among its
siblings in the SAME maternal genet x pollinator treatment family. A pure
shared-mother or family-structured association survives this permutation;
only child-specific multilocus genetic information can exceed it.
A second diagnostic permutes within mother, across treatments, to see
whether treatment-specific gene flow matters. Neither is a true sire label.

No locations or parentage priors; parent, offspring IDs and individual
SSR allele sizes never printed, written, or committed. All public outputs
are statistical aggregates. This is reanalysis of the same Rhodes dataset.
"""
from __future__ import annotations
import argparse,collections,json
from pathlib import Path
import numpy as np
from scipy.special import softmax
from oenothera_v57_genotype_only_paternity import source
from oenothera_v58_locus_robustness import probabilities, SCORE_EPS

OUT=Path("analysis_results_oenothera_child_specific_genetic_signal_v59")
SEED=20261124
PERM_DRAWS=2999
MOTHER_BOOT=9999
LOC=7
EXPECTED_V58={
  .002:0.6518456957530187,
  .02:0.5860018809991241,
  .1:0.1617373926266094,
}

def grouping(d,mode):
    mi=d["mother_index"]
    if mode=="same_mother_same_treatment":
        keys=[(int(mi[i]),d["child_treatment"][i]) for i in range(len(mi))]
    elif mode=="same_mother_any_treatment":
        keys=[(int(mi[i]),) for i in range(len(mi))]
    else:
        raise ValueError(mode)
    ix={}
    for i,key in enumerate(keys):ix.setdefault(key,[]).append(i)
    groups=[np.asarray(rows,dtype=int) for rows in ix.values()]
    assert sum(len(r) for r in groups)==521
    mothers=np.unique(mi)
    assert len(mothers)==26
    mids={int(m):i for i,m in enumerate(mothers)}
    groups=[(mids[int(mi[idx[0]])],idx) for idx in groups]
    return groups, {
        "grouping":mode,
        "groups":len(groups),
        "groups_with_two_or_more_children":sum(len(idx)>=2 for m,idx in groups),
        "singletons":sum(len(idx)==1 for m,idx in groups),
        "mothers":len(mothers),
        "group_size_quantiles":np.quantile([len(idx) for m,idx in groups],[0,.25,.5,.75,1]).tolist()
    }

def one_scenario(d,trans,mc,eps,mode,groups,rng):
    n,f,l=trans.shape
    assert (n,f,l)==(521,323,7)
    mi=d["mother_index"]
    mask=np.arange(f)[None,:]!=mi[:,None]
    assert np.all(mask.sum(axis=1)==322)
    logp=np.log(trans+eps)
    full=logp.sum(axis=2)
    nmom=26
    orig=np.zeros(nmom,float)
    expnull=np.zeros(nmom,float)
    nvalid=np.zeros(nmom,dtype=int)
    null=np.zeros((PERM_DRAWS,nmom),float)
    by_locus_orig=np.zeros(LOC,float)
    by_locus_null=np.zeros(LOC,float)
    by_locus_count=np.zeros(LOC,dtype=int)
    for loc in range(LOC):
        # Six discovery loci only: heldout locus is never used for weights.
        weights=softmax(np.where(mask,full-logp[:,:,loc],-1e8),axis=1)
        P=trans[:,:,loc]+eps
        ref=(P*mask).sum(axis=1)/322
        assert np.all(ref>0)
        for mom_id,ids in groups:
            w=weights[ids]
            target=P[ids]  # genes from the SAME mother (within same treatment primary)
            valid=mc[ids,loc].astype(int)
            m=len(ids)
            assert np.all(mi[ids]==mi[ids[0]])
            # score matrix (discovery child i, seventh-locus target child j).
            # Each column's reference uses the same target child j.
            num=w@target.T
            score=np.log(np.maximum(num,1e-300)/ref[ids][None,:])
            diag=np.diag(score)
            orig[mom_id]+=float(np.sum(diag*valid))
            expnull[mom_id]+=float(np.sum(score*valid[None,:])/m)
            nvalid[mom_id]+=int(valid.sum())
            by_locus_orig[loc]+=float(np.sum(diag*valid))
            by_locus_null[loc]+=float(np.sum(score*valid[None,:])/m)
            by_locus_count[loc]+=int(valid.sum())
            if m<2:
                null[:,mom_id]+=float(np.sum(diag*valid))
                continue
            # An exact uniform random permutation in each maternal sibship.
            permutation=np.argsort(rng.random((PERM_DRAWS,m)),axis=1)
            rows=np.arange(m)[None,:]
            draw_score=score[rows,permutation]
            null[:,mom_id]+=(draw_score*valid[permutation]).sum(axis=1)
    assert (nvalid>0).all() and nvalid.sum()==3445,("unexpected typed child-locus counts",nvalid.sum())
    observed=orig/nvalid
    expected=expnull/nvalid
    randomized=null/nvalid[None,:]
    obs=float(observed.mean())
    ex=float(expected.mean())
    diff=observed-expected
    permdiff=float(obs-ex)
    assert abs(randomized.mean()-ex)<.015
    sim=randomized.mean(axis=1)
    # Mother-genet bootstrap for effect size: sibling pairs / seven loci are
    # *not* independent population-level replicates.
    draw_ix=rng.integers(nmom,size=(MOTHER_BOOT,nmom))
    boot=diff[draw_ix].mean(axis=1)
    tail=float((1+np.count_nonzero(sim>=obs))/(len(sim)+1))
    return {
       "epsilon":float(eps),
       "grouping":mode,
       "number_independent_maternal_genets":nmom,
       "child_locus_pairs_evaluable":int(nvalid.sum()),
       "original_heldout_log_score_nats_mother_equal":obs,
       "within_group_shuffled_expected_log_score_nats_mother_equal":ex,
       "child_specific_increment_over_mother_group_null_nats":permdiff,
       "mother_bootstrap_95_increment":np.quantile(boot,[.025,.975]).tolist(),
       "mothers_with_positive_increment":int((diff>0).sum()),
       "mothers_with_negative_increment":int((diff<0).sum()),
       "mothers_with_zero_increment":int((diff==0).sum()),
       "null_exact_within_group_uniform_permutations":PERM_DRAWS,
       "permutation_null_score_mean":float(sim.mean()),
       "permutation_null_score95":np.quantile(sim,[.025,.975]).tolist(),
       "conditional_permutation_upper_tail_p":tail,
       "permutation_null_mother_equal_score_sd":float(sim.std(ddof=1)),
       "per_locus_observed_log_gain_vs_uniform_fathers":[float(by_locus_orig[j]/by_locus_count[j]) for j in range(LOC)],
       "per_locus_sib_permutation_expected_log_gain_vs_uniform_fathers":[float(by_locus_null[j]/by_locus_count[j]) for j in range(LOC)],
       "source_mendelian_father_scores_not_calibrated_posterior":True,
       "marginal_maternal_genotype_and_treatment_unchanged_in_primary":mode=="same_mother_same_treatment",
    }

def main():
    d=source()
    trans,inform,mc,mi_mismatch=probabilities(d)
    assert int(mi_mismatch.sum())==1
    allresults={}
    rng=np.random.default_rng(SEED)
    for mode in ["same_mother_same_treatment","same_mother_any_treatment"]:
        groups,gmeta=grouping(d,mode)
        out=[]
        for eps in SCORE_EPS:
            r=one_scenario(d,trans,mc,eps,mode,groups,rng)
            # reproduce v58 primary without relying on a generated report.
            assert np.isclose(r["original_heldout_log_score_nats_mother_equal"],
                EXPECTED_V58[eps],atol=1e-10,rtol=0)
            out.append(r)
        allresults[mode]={"design":gmeta,"epsilon_sensitivity":out}
    return {
       "status":"original_Oenothera_v59_within_family_child_allele_scrambling_negative_control",
       "source_paper":"10.1111/mec.14115",
       "source_Dryad":"10.5061/dryad.p24q3",
       "source_fixed_Zenodo":"4942351",
       "source_genotypes_sha256":d["original_genotype_sha256"],
       "source_correlated_paternity_sha256":d["original_paternity_sha256"],
       "parent_genets":323,"offspring_genotyped":521,"independent_mothers":26,
       "loci":7,"informative_mother_offspring_locus_pairs":int(mc.sum()),
       "seed":SEED,
       "results":allresults,
       "limitations":{
          "original_2017_source_reused_not_new_field_experiment":True,
          "within_group_permutation_conditions_on_maternal_genotype_and_sibling_treatment":True,
          "excess_six_to_seventh_locus_prediction_not_proof_of_identified_true_father":True,
          "all_candidate_sires_sampled_adults_only_unobserved_fathers_remain":True,
          "SSR_loci_potential_population_structure_and_relatedness_not_eliminated":True,
          "random_within_family_swap_conservative_if_siblings_share_single_father":True,
          "nonzero_allele_code_valid_zero_is_missing":True,
          "no_geographic_distance_or_spatial_prior_input":True,
          "mother_bootstrap_not_seven_loci_independent":True,
          "no_original_alleles_coordinates_parentage_rows_or_ids_exported":True,
          "does_not_reidentify_Daphne_jezoensis_or_D_pseudomezereum_parentage":True}}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT/"oenothera_v59_within_family_allele_scramble.json")
    args=p.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    res=main()
    args.output.write_text(json.dumps(res,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("OENOTHERA_V59_CHILD_PERMUTATION",json.dumps(res,ensure_ascii=False,allow_nan=False),flush=True)
