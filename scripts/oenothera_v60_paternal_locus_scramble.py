#!/usr/bin/env python3
"""v0.60: Does heldout-SSR signal require coherent seven-locus paternal
genotypes, rather than merely per-locus allele frequencies?

Source: Rhodes, Fant & Skogen (2017) *Oenothera harringtonii*, Molecular
Ecology doi 10.1111/mec.14115, original Dryad 10.5061/dryad.p24q3,
SHA-pinned Zenodo 4942351 (reused v57-59 original sample).

Original father candidates: 323 diploid seven-SSR adult genotypes; 521
offspring typed from 26 mothers in 68 mother x pollination-treatment families.
Positive control: shuffle only NON-MATERNAL adult father labels with exactly
the SAME permutation across all seven loci: all posterior-free sire-candidate
scores must be invariant (mere renaming of candidate fathers).
Negative control: independently permute each diploid paternal locus's entire
323-candidate label mapping, keeping all 26 *actual observed maternal genets*
fixed as possible sires and their known maternal genotypes untouched. This
preserves parental allele-frequency / within-locus diplotype distributions,
the sampled maternal genotypes, original mother-offspring pedigrees, offspring
treatment and all missingness per locus. It destroys identity coherence across
paternal loci for 297 non-maternal adult candidates.

Also shuffle independently WITHIN source parent alphabetic-ID prefix strata
for adult paternal donors, a coarse observational structure control. Prefix
groups are not verified habitat/population strata; no biological label used.

For each permutation, re-evaluate 6-to-7 heldout likelihood ratio relative to
uniform candidate fathers without any coordinates or spatial prior. Uncalibrated
Mendelian smoothing epsilon .002/.02/.1. The 26 mother genets are the
statistical resampling units, not 521 progeny x 7 markers or permutations.
This is NOT independently known sire identity, a true sire posterior, a
multisite/phylogenetic replication or pollinator behavior causal inference.

All original individual IDs, SSR allele sizes, parentage records and
locations stay in memory; logs and JSON are aggregate-only.
"""
from __future__ import annotations
import argparse,collections,json,re
from pathlib import Path
import numpy as np
from scipy.special import softmax
from oenothera_v57_genotype_only_paternity import source
from oenothera_v57_source_schema import download,FILES
from oenothera_v57_genotype_key_audit import records
from oenothera_v58_locus_robustness import probabilities

OUT=Path("analysis_results_oenothera_paternal_multilocus_v60")
SEED=20261124
PERMUTATIONS=399
EPSILONS=(.002,.02,.1)
N_MOTHERS=26
KNOWN_V58={.002:.6518456957530187,.02:.5860018809991241,.1:.1617373926266094}
MOTHER_BOOT=9999

def source_sets(d):
    mi=np.asarray(d["mother_index"],dtype=int)
    mothers=np.unique(mi)
    assert len(mothers)==N_MOTHERS
    all_adults=np.arange(len(d["parent_genotype"]),dtype=int)
    free=all_adults[~np.isin(all_adults,mothers)]
    assert len(free)==297
    parentrows=records(download("genotypes.xlsx",FILES["genotypes.xlsx"]))["parents"]
    names=[str(row[0]) for row in parentrows]
    assert len(names)==323
    prefix=[re.match(r"^[A-Za-z]+",name).group(0) for name in names]
    assert all(prefix)
    grouped=collections.defaultdict(list)
    for idx in free:grouped[prefix[int(idx)]].append(int(idx))
    strata=[np.array(ids,dtype=int) for ids in grouped.values()]
    assert np.array_equal(np.sort(np.concatenate(strata)),free)
    return mothers,free,strata,{
      "total_adults":len(all_adults),
      "maternal_known_adults_held_fixed":len(mothers),
      "nonmaternal_adult_candidates_permuted":len(free),
      "nonmaternal_prefix_groups":len(strata),
      "nonmaternal_prefix_group_size_quantiles":np.quantile([len(x) for x in strata],[0,.25,.5,.75,1]).tolist(),
      "unshufflable_nonmaternal_prefix_groups":int(sum(len(x)<2 for x in strata)),
      "fraction_nonmaternal_candidate_rows_in_prefix_groups_ge2":
          float(sum(len(x) for x in strata if len(x)>=2)/len(free)),
      "no_sensitive_source_prefix_labels_emitted":True
    }

def assign_perm(free,strata,rng,mode):
    perm=np.arange(323,dtype=int)
    if mode=="global":
        perm[free]=rng.permutation(free)
    elif mode=="within_prefix":
        for group in strata:
            perm[group]=rng.permutation(group)
    else:
        raise ValueError(mode)
    assert set(perm.tolist())==set(range(323))
    return perm

def heldout_score(t,mc,mi,eps,momrows):
    # Score unlabeled father candidates exactly as in v58.
    n,k,l=t.shape
    assert (n,k,l)==(521,323,7)
    mask=np.arange(k)[None,:]!=mi[:,None]
    s=np.log(t+eps)
    full=s.sum(axis=2)
    per_child=np.zeros(n,dtype=float)
    for locus in range(l):
        w=softmax(np.where(mask,full-s[:,:,locus],-1e8),axis=1)
        test=t[:,:,locus]+eps
        pred=(w*test).sum(axis=1)
        ref=(test*mask).sum(axis=1)/322
        g=np.log(pred/ref)
        assert np.isfinite(g).all()
        per_child+=np.where(mc[:,locus],g,0)
    gm=np.array([per_child[ids].sum()/max(1,mc[ids].sum()) for ids in momrows])
    assert len(gm)==26
    return gm

def main():
    d=source()
    m=np.asarray(d["mother_index"],int)
    t,typed,mc,wrong=probabilities(d)
    assert int(wrong.sum())==1
    moms,free,strata,meta=source_sets(d)
    mother_ix=[np.flatnonzero(m==x) for x in moms]
    for ix in mother_ix:
        assert len(ix)>0 and mc[ix].sum()>0
    observed=np.vstack([heldout_score(t,mc,m,eps,mother_ix) for eps in EPSILONS])
    assert np.allclose(observed.mean(axis=1),
          np.array([KNOWN_V58[e] for e in EPSILONS]),rtol=0,atol=1e-12)
    rng=np.random.default_rng(SEED)
    # Metamorphic identity check: renaming nonmaternal father labels with ONE
    # shared permutation across all seven loci is mathematically invariant.
    consistency_perm=assign_perm(free,strata,rng,"global")
    exact=t[:,consistency_perm,:]
    meta_errors=[float(np.max(np.abs(
      heldout_score(exact,mc,m,eps,mother_ix)-observed[i])))
        for i,eps in enumerate(EPSILONS)]
    assert max(meta_errors)<1e-12,("same coherent father permutation changed predictions",meta_errors)
    results={}
    for mode in ("global","within_prefix"):
        simulated=np.empty((PERMUTATIONS,len(EPSILONS),N_MOTHERS),dtype=float)
        n_changed=np.empty(PERMUTATIONS,dtype=float)
        for draw in range(PERMUTATIONS):
            # independent donor identities per SSR locus, preserving an
            # intact pair of paternal alleles *within* each locus.
            per_locus=[assign_perm(free,strata,rng,mode) for _ in range(7)]
            n_changed[draw]=np.mean([float(np.mean(p[free]!=free)) for p in per_locus])
            shuf=np.stack([t[:,p,locus] for locus,p in enumerate(per_locus)],axis=2)
            assert np.array_equal(shuf[:,moms,:],t[:,moms,:]),"original mother genotypes must remain unchanged in father candidate set"
            for i,eps in enumerate(EPSILONS):
                simulated[draw,i,:]=heldout_score(shuf,mc,m,eps,mother_ix)
        cas=[]
        for i,eps in enumerate(EPSILONS):
            mom=observed[i]
            neg=simulated[:,i,:]
            per_draw=neg.mean(axis=1)
            avg_per_mother=neg.mean(axis=0)
            increments=mom-avg_per_mother
            obs=float(mom.mean())
            null=float(per_draw.mean())
            B=rng.integers(N_MOTHERS,size=(MOTHER_BOOT,N_MOTHERS))
            boot=increments[B].mean(axis=1)
            p=float((1+(per_draw>=obs).sum())/(PERMUTATIONS+1))
            cas.append({
              "epsilon":float(eps),
              "source_original_evaluation_mean_nats":obs,
              "independently_shuffled_sire_locus_null_mean_nats":null,
              "child_specific_paternal_identity_coherence_increment_nats":obs-null,
              "null_95percent_across_locus_permutes":np.quantile(per_draw,[.025,.975]).tolist(),
              "null_max":float(per_draw.max()),
              "conditional_null_upper_p":p,
              "number_null_locus_permutations":PERMUTATIONS,
              "equal_mother_resample_95increment":np.quantile(boot,[.025,.975]).tolist(),
              "mothers_observed_minus_null_positive":int((increments>0).sum()),
              "mothers_observed_minus_null_negative":int((increments<0).sum()),
              "mothers_observed_minus_null_zero":int((increments==0).sum()),
              "no_paternity_truth_labels_available":True,
              "prior_free_genotype_only_smoothing_uncalibrated":True
            })
        results[mode]={
          "null_type":"independent whole-diploid paternal locus genotype permutations among nonmaternal adults; original maternal donors frozen",
          "average_fraction_nonmaternal_sire_label_moved_per_locus":float(n_changed.mean()),
          "epsilon_sensitivity":cas}
    return {
      "status":"source_Oenothera_SSR_paternal_multilocus_identity_coherence_v60",
      "source_study":"10.1111/mec.14115",
      "source_Dryad":"10.5061/dryad.p24q3",
      "source_fixed_Zenodo_record":4942351,
      "source_genotypes_sha256":d["original_genotype_sha256"],
      "source_paternity_reference_sha256":d["original_paternity_sha256"],
      "parents":323,"offspring":521,"maternal_genets":26,"source_SSR_loci":7,
      "maternal_offspring_informative_locus_pairs":int(mc.sum()),
      "permutation_design":meta,
      "coherent_permutation_metmorphic_max_absolute_score_difference":max(meta_errors),
      "coherent_permutation_test_passed":True,
      "results":results,
      "limits":{
        "one_2017_published_Oenothera_dataset":True,
        "no_Daphne_offspring_genotypes_used_or_new_sire_reassignment":True,
        "paternal_diploid_allele_pair_per_marker_kept_intact":True,
        "all_26_maternal_source_donor_genotypes_held_fixed":True,
        "paternal_marker_exchange_removes_individual_cross_locus_identity_only_for_297_nother_genets":True,
        "original_maternal_child_genetic_profiles_and_treatments_unchanged":True,
        "prefix_groups_unverified_as_true_genetic_or_geographic_population":True,
        "marker_identity_predictive_signal_could_include_genetic_relatedness_and_population_structure":True,
        "no_true_paternity_posterior_or_survival_recruitment":True,
        "mother_genet_resampling_not_offspring_independent":True,
        "no_original_identifiers_alleles_genotypes_or_coords_output":True,
        "focal_D_pseudomezereum_not_measured":True}}
if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--output",type=Path,default=OUT/"oenothera_paternal_multilocus_coherence_v60.json")
    arg=cli.parse_args()
    arg.output.parent.mkdir(parents=True,exist_ok=True)
    v=main()
    arg.output.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("OENOTHERA_V60_PATERNAL_LOCUS_PERMUTE",json.dumps(v,ensure_ascii=False,allow_nan=False),flush=True)
