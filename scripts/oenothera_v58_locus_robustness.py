#!/usr/bin/env python3
"""v0.58: 7-SSR leave-one-locus-out and negative-control audit on the same
independent published Oenothera harringtonii family data (Rhodes et al. 2017,
Dryad p24q3; exact pinned Zenodo 4942351; v57 verified reader).

Questions (do NOT misrepresent as calibrated sire posteriors):
1. Do strictly unique Mendelian compatible sampled sires persist after any
   one of seven loci is withheld?
2. Does genotype-only identity evidence from six loci predict the held-out
   seventh offspring allele, over an all-sampled-father uniform baseline?
3. Would arbitrary wrong mother labels give many more Mendelian conflicts?
The independent unit for score uncertainty is the maternal genet (26), not
521 seedlings x 7 loci. Preserve no per-offspring/parent/mother identifiers,
SSR alleles, exact locations or raw data on GitHub or in CI artifacts.

No locations/distance prior, no behavioral inference, no survival data.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.special import softmax
from oenothera_v57_genotype_only_paternity import source, LOC

OUT=Path("analysis_results_oenothera_locus_robustness_v58")
SEED=20261124
SCORE_EPS=[.002,.02,.1]
MOTHER_BOOT=7999
MOTHER_NEGATIVE_DRAWS=1999

def probabilities(dat):
    """Allele-specific Mendelian transmission P(child genotype | mom, dad),
    computed only when both alleles of all three are actually typed.
    Missing alleles are 0-coded. Uninformative father is assigned mean
    sample likelihood and not counted in the strict compatibility gate.
    """
    pp=dat["parent_genotype"]
    cc=dat["offspring_genotype"]
    mi=dat["mother_index"].astype(int)
    n,k=len(cc),len(pp)
    assert pp.shape==(323,7,2) and cc.shape==(521,7,2)
    assert len(mi)==n
    mom=pp[mi]
    fully_mc=(mom!=0).all(axis=2)&(cc!=0).all(axis=2)
    trans=np.ones((n,k,LOC),dtype=np.float64)
    inform=np.zeros((n,k,LOC),dtype=bool)
    incompat_mom=np.zeros((n,LOC),dtype=bool)
    for j in range(LOC):
        m=mom[:,j,:]
        c=cc[:,j,:]
        f=pp[:,j,:]
        has=fully_mc[:,j]
        incompat_mom[:,j]=has & (~((m[:,0]==c[:,0])|(m[:,0]==c[:,1])|
                                    (m[:,1]==c[:,0])|(m[:,1]==c[:,1])))
        p=np.zeros((n,k),dtype=float)
        for a in range(2):
            for b in range(2):
                p+=.25*((((m[:,a,None]==c[:,0,None]) & (f[None,:,b]==c[:,1,None]))|
                           ((m[:,a,None]==c[:,1,None]) & (f[None,:,b]==c[:,0,None]))))
        father_typed=(f!=0).all(axis=1)
        ok=has[:,None]&father_typed[None,:]
        inform[:,:,j]=ok
        # mean over actually typed candidates in the population, not arbitrary
        # assumption that an untyped potential sire matches the offspring.
        # For MC-missing, locus contributes identically to all candidates.
        n_known=np.sum(ok,axis=1)
        bkg=np.divide(np.sum(p*ok,axis=1),np.maximum(1,n_known))
        trans[:,:,j]=np.where(ok,p,np.where(has[:,None],bkg[:,None],1.))
    assert int(incompat_mom.sum())==1, "v57 maternal source sanity"
    return trans,inform,fully_mc,incompat_mom

def eligible_candidates(trans,inform,mother_ix,min_loci=5):
    assert trans.shape==inform.shape
    valid=inform.sum(axis=2)>=min_loci
    bad=(inform & (trans==0)).sum(axis=2)
    ok=valid&(bad==0)
    ok[np.arange(len(ok)),mother_ix]=False
    return ok

def leaveout(trans,inform,mi):
    full=eligible_candidates(trans,inform,mi,5)
    n=full.sum(axis=1)
    unique=full.argmax(axis=1)
    truly=n==1
    by_locus=[]
    same=np.zeros(len(mi),dtype=int)
    robust_unique=np.ones(len(mi),dtype=bool)&truly
    c_alt=np.zeros(len(mi),dtype=int)
    # With 1/7 masked, require four remaining informative loci.
    # Alternative stricter 5/6 support is separately reported.
    for j in range(LOC):
        pm=trans[:,:,np.arange(LOC)!=j]
        im=inform[:,:,np.arange(LOC)!=j]
        compat4=eligible_candidates(pm,im,mi,4)
        compat5=eligible_candidates(pm,im,mi,5)
        ct4=compat4.sum(axis=1)
        ct5=compat5.sum(axis=1)
        same_candidate4=compat4[np.arange(len(mi)),unique]
        stable4=truly&(ct4==1)&same_candidate4
        stable5=truly&(ct5==1)&compat5[np.arange(len(mi)),unique]
        same+=stable4
        robust_unique&=stable4
        c_alt+=(ct4>1)&truly
        by_locus.append({
          "locus_held_out":j+1,
          "offspring_unique_sire_with_four_complete_loci":int((ct4==1).sum()),
          "offspring_zero_candidate_with_four_complete_loci":int((ct4==0).sum()),
          "original_unique_still_same_unique_with_four_complete_loci":int(stable4.sum()),
          "original_unique_still_same_unique_with_five_complete_loci":int(stable5.sum()),
          "original_unique_ambiguous_after_drop":int((truly&(ct4>1)).sum()),
          "original_unique_loses_candidate_after_drop":int((truly&(ct4==0)).sum())
        })
    return {
      "n_offspring":len(mi),
      "full_7_locus_unique_count_5_typed_min":int(truly.sum()),
      "full_7_locus_zero_candidate":int((n==0).sum()),
      "full_7_locus_multiple_candidate":int((n>1).sum()),
      "full_unique_same_sire_under_all_seven_single_locus_omissions_4_typed_min":int(robust_unique.sum()),
      "full_unique_single_locus_omission_stability_count_distribution":
          {str(i):int((same[truly]==i).sum()) for i in range(LOC+1)},
      "full_unique_ever_ambiguous_after_one_locus_drop":int((c_alt>0).sum()),
      "unique_sire_stability_by_omitted_locus":by_locus,
      "not_genetic_paternity_certainty_or_parentage_posterior":True,
    }

def holdout(trans,inform,mc,mi,mother_boot=MOTHER_BOOT):
    """Predict an *unseen* seventh SSR genotype with six-locus weights.
    Scores are log ratios vs exchangeable sampled adult fathers, no geography.
    Maternal overlap and completely untyped locus excluded. A score can be
    nonzero despite no perfect sampled father (soft epsilon scenario).
    """
    n,k,l=trans.shape
    assert (n,k,l)==(521,323,7)
    parents=np.arange(k)[None,:]!=mi[:,None]
    n_parent=parents.sum(axis=1)
    assert np.all(n_parent==322)
    held_ok=mc.astype(bool)
    arr=np.zeros((len(SCORE_EPS),n,LOC),dtype=float)
    for e,eps in enumerate(SCORE_EPS):
        logp=np.log(trans+eps)
        full=logp.sum(axis=2)
        for loc in range(LOC):
            lw=full-logp[:,:,loc]
            lw=np.where(parents,lw,-1e8)
            w=softmax(lw,axis=1)
            test=trans[:,:,loc]+eps
            pred=np.sum(w*test,axis=1)
            ref=np.sum(test*parents,axis=1)/n_parent
            delta=np.log(pred/ref)
            assert np.isfinite(delta).all()
            arr[e,:,loc]=np.where(held_ok[:,loc],delta,0)
    moms=np.unique(mi)
    assert len(moms)==26
    momrows=[np.flatnonzero(mi==m) for m in moms]
    rng=np.random.default_rng(SEED)
    B=rng.integers(len(moms),size=(mother_boot,len(moms)))
    results=[]
    for e,eps in enumerate(SCORE_EPS):
        valid=held_ok
        # mother-mean is the first-order unit. one child with all 7 loci
        # should not count as seven independent mothers.
        gm=np.array([arr[e,idx].sum()/max(1,valid[idx].sum()) for idx in momrows])
        boots=gm[B].mean(axis=1)
        locs=[]
        for locus in range(LOC):
            ix=valid[:,locus]
            vals=arr[e,ix,locus]
            group=np.array([arr[e,ix_m,locus][valid[ix_m,locus]].mean()
                if valid[ix_m,locus].any() else np.nan for ix_m in momrows])
            group=group[np.isfinite(group)]
            locs.append({
              "locus_held_out":locus+1,"evaluable_offspring":int(ix.sum()),
              "mean_log_heldout_likelihood_gain_over_uniform_sire":float(vals.mean()),
              "mother_mean_gain":float(group.mean()),
              "n_informative_mothers":len(group),
            })
        results.append({
         "epsilon":float(eps),
         "n_maternal_genets":len(moms),
         "informative_offspring_locus_pairs":int(valid.sum()),
         "child_locus_mean_log_likelihood_ratio":float(arr[e].sum()/valid.sum()),
         "equal_mother_mean_log_likelihood_ratio":float(gm.mean()),
         "mother_bootstrap95":np.quantile(boots,[.025,.975]).tolist(),
         "maternal_genets_with_positive_gain":int((gm>0).sum()),
         "maternal_genets_with_negative_gain":int((gm<0).sum()),
         "by_held_out_locus":locs})
    return results

def mother_negative_control(dat,fully_mc,incompat,draws=MOTHER_NEGATIVE_DRAWS):
    p=dat["parent_genotype"];c=dat["offspring_genotype"];mi=dat["mother_index"]
    true_bad=int(incompat.any(axis=1).sum())
    assert true_bad==1
    eligible_mom=np.unique(mi)
    rng=np.random.default_rng(SEED+188)
    bads=np.empty(draws,dtype=int);mismatches=np.empty(draws,dtype=int)
    true_loci=int(fully_mc.sum())
    for k in range(draws):
        # Swap real mothers among observed maternal-genet pool, forcing each
        # child to be paired with a different original mother genotype.
        wrong=eligible_mom[rng.integers(len(eligible_mom)-1,size=len(mi))]
        wrong+=(wrong>=mi) if False else 0 # no index arithmetic on parent ID
        # Correct exclusion by indexing the *position in the sorted mother
        # pool*, not raw 323-parent indices.
        eligible_pos=np.searchsorted(eligible_mom,mi)
        draw_pos=rng.integers(len(eligible_mom)-1,size=len(mi))
        draw_pos+=draw_pos>=eligible_pos
        wrong=eligible_mom[draw_pos]
        mm=p[wrong]
        ok=(mm!=0).all(axis=2)&(c!=0).all(axis=2)
        compat=(mm[:,:,0]==c[:,:,0])|(mm[:,:,0]==c[:,:,1])|(
                 mm[:,:,1]==c[:,:,0])|(mm[:,:,1]==c[:,:,1])
        mismatch=ok&(~compat)
        bads[k]=np.count_nonzero(mismatch.any(axis=1))
        mismatches[k]=int(mismatch.sum())
    return {
       "actual_mother_offspring_with_one_or_more_mendelian_exclusions":true_bad,
       "wrong_mother_randomized_offspring_with_mendelian_exclusions_mean":float(bads.mean()),
       "wrong_mother_randomized_offspring_with_mendelian_exclusions95":np.quantile(bads,[.025,.975]).tolist(),
       "randomized_min":int(bads.min()),
       "randomized_draws":draws,
       "wrong_mother_randomized_incompatible_loci_total_mean":float(mismatches.mean()),
       "original_mother_child_completely_typed_locus_pairs":true_loci,
       "randomized_p_lower_one_sided":float((1+(bads<=true_bad).sum())/(draws+1)),
       "no_raw_mother_identity_emitted":True}

def main():
    d=source()
    p,inf,mc,wrong=probabilities(d)
    return {
       "status":"independent_Oenothera_leave_one_locus_out_and_maternal_negative_control_v58",
       "source_doi":"10.5061/dryad.p24q3",
       "source_published_article":"10.1111/mec.14115",
       "pinned_genotype_sha256":d["original_genotype_sha256"],
       "pinned_paternity_sha256":d["original_paternity_sha256"],
       "genotype_loci":LOC,"source_mothers":int(len(np.unique(d["mother_index"]))),
       "strict_mendelian_uniqueness_stability":leaveout(p,inf,d["mother_index"]),
       "heldout_locus_predictive_evidence":holdout(p,inf,mc,d["mother_index"]),
       "wrong_mother_negative_control":mother_negative_control(d,mc,wrong),
       "limits":{
          "off_locus_prediction_is_out_of_marker_not_out_of_population":True,
          "no_spatial_or_distance_prior_in_model":True,
          "missing_zero_alleles_ignored_or_background_locus_frequency":True,
          "soft_eps_weight_not_calibrated_paternity_posterior":True,
          "genotyped_sampled_adults_only_no_external_fathers":True,
          "marker_dropout_stability_not_survival_paternity_validation":True,
          "random_wrong_mother_genotypes_may_capture_real_relatives":True,
          "one_original_source_paper_not_new_field_replicate":True,
          "no_original_genotypes_ID_or_coordinates_written":True,
          "focal_Daphne_pseudomezereum_not_measured":True}}
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=OUT/"oenothera_locus_robustness_v58.json")
    o=ap.parse_args()
    o.output.parent.mkdir(parents=True,exist_ok=True)
    res=main()
    o.output.write_text(json.dumps(res,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("OENOTHERA_V58_LOLO_AND_HOLDOUT",json.dumps(res,ensure_ascii=False,allow_nan=False),flush=True)
