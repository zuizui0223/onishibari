#!/usr/bin/env python3
"""v0.57 independently typed Oenothera harringtonii parent/offspring (7 SSR).
Distance-free Mendelian parental compatibility and pollinator-group comparison.

Same original Rhodes et al. 2017 Mol Ecol dataset, Dryad DOI 10.5061/dryad.p24q3,
official preservation Zenodo 4942351 source exact archive MD5/SHA checked.

IMPORTANT: Not a father posteriors model calibrated to CERVUS or MLTR.
Complete-genotype exclusion is deterministic Mendelian compatibility. A
"soft" compatibility weight below is only a sensitivity rank, NOT validated
paternity confidence. Candidate fathers are all sampled adults minus known
mother; unknown/unsampled fathers, genotyping errors, and population allele
frequencies have not been modeled. Never report individual parental or
offspring IDs, allele genotypes or sample-level parentage.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,pearsonr
from oenothera_v57_source_schema import download,FILES
from oenothera_v57_genotype_key_audit import records

OUT=Path("analysis_results_independent_oenothera_v57")
LOC=7
EPS_LIST=(0.002,0.02,0.10)
TRT_MAP={"C":"c","D":"de","N":"ne"}

def transform(row):
    assert len(row)==15
    assert all(str(x).strip() for x in row[1:15])
    return np.array([int(float(x)) for x in row[1:15]],dtype=np.int32).reshape(LOC,2)

def source():
    geno=download("genotypes.xlsx",FILES["genotypes.xlsx"])
    dry=download("multiplePaternity.csv",FILES["multiplePaternity.csv"])
    tbl=records(geno)
    parents=tbl["parents"]; children=tbl["offspring"]
    parent_raw=[str(r[0]) for r in parents]
    parent4=[id[-4:] for id in parent_raw]
    assert len(set(parent4))==len(parent4)==323, "Four-digit parent suffix must be globally unique"
    map4={id:i for i,id in enumerate(parent4)}
    genop=np.stack([transform(r) for r in parents])
    genoc=np.stack([transform(r) for r in children])
    m=[];trt=[];mother_key=[];offsids=[]
    for r in children:
        sid=str(r[0])
        i=sid.split("_")
        assert len(i)==2 and len(i[0])==4 and len(i[1])==2
        assert i[0] in map4 and i[1][0].upper() in TRT_MAP
        m.append(map4[i[0]])
        mother_key.append(i[0])
        trt.append(TRT_MAP[i[1][0].upper()])
        offsids.append(sid)
    assert len(offsids)==len(set(offsids))==521
    label=dry.decode("utf-8-sig").replace(chr(13)+chr(10),chr(10)).replace(chr(13),chr(10))
    meta=list(csv.DictReader(io.StringIO(label,newline="")))
    assert len(meta)==60
    family_ref={(r["plantID"],r["treatment"]):float(r["correlatedPaternity"]) for r in meta}
    assert len(family_ref)==60
    return {
      "parent_genotype":genop,"offspring_genotype":genoc,"mother_index":np.array(m),
      "family":list(zip(mother_key,trt)),"source_ref":family_ref,
      "child_treatment":np.array(trt),
      "original_genotype_sha256":hashlib.sha256(geno).hexdigest(),
      "original_paternity_sha256":hashlib.sha256(dry).hexdigest(),
    }

def meiosis_probability(m,f,c):
    # Exactly four equiprobable gamete pairings (two maternal x two paternal).
    # If a genotype is homozygous, duplicated gametes count twice correctly.
    p=0
    for ma in m:
        for fa in f:
            if (ma==c[0] and fa==c[1]) or (ma==c[1] and fa==c[0]):
                p+=.25
    return p

def primary(d):
    par=d["parent_genotype"]
    child=d["offspring_genotype"]
    mother=d["mother_index"]
    n_par=len(par);n_off=len(child)
    # 0 is an untyped SSR allele sentinel in this published workbook, NOT
    # a biological microsatellite size. Both mother/offspring alleles must
    # be nonzero before exact Mendelian exclusion is meaningful.
    exact=np.ones((n_off,n_par,LOC),dtype=float)
    father_locus_informative=np.zeros((n_off,n_par,LOC),dtype=bool)
    mother_child_locus_informative=np.zeros((n_off,LOC),dtype=bool)
    mother_inconsistent=np.zeros((n_off,LOC),dtype=bool)
    for ix in range(n_off):
        mom=par[mother[ix]]
        for locus in range(LOC):
            c=child[ix,locus]
            if 0 in mom[locus] or 0 in c:
                continue
            mother_child_locus_informative[ix,locus]=True
            mother_inconsistent[ix,locus]=not np.any(c[:,None]==mom[locus][None,:])
            typed_fathers=[]
            for f in range(n_par):
                if 0 in par[f,locus]:
                    continue
                father_locus_informative[ix,f,locus]=True
                prob=meiosis_probability(mom[locus],par[f,locus],c)
                exact[ix,f,locus]=prob
                if f!=mother[ix]:
                    typed_fathers.append(prob)
            # Missing paternal allele = uninformative locus, not a perfect
            # paternal genetic match. Replace by random sampled-father mean
            # likelihood (source population) so missing genotypes cannot gain
            # a favorable sire-score merely by being untyped.
            background=float(np.mean(typed_fathers)) if typed_fathers else 0.
            for f in range(n_par):
                if not father_locus_informative[ix,f,locus]:
                    exact[ix,f,locus]=background
        exact[ix,mother[ix],:]=0
    # Strict Mendelian exclusion requires at least five informative loci
    # in both child/mother and candidate sire; partial founder genotypes
    # are not mistaken for 0-mismatch certainty.
    zeros=((exact==0) & father_locus_informative).sum(axis=2)
    usable=father_locus_informative.sum(axis=2)
    zeros[np.arange(n_off),mother]=LOC+1
    usable[np.arange(n_off),mother]=0
    valid=(usable>=5)
    compat=(zeros==0)&valid
    near=(zeros<=1)&valid
    compatible=compat.sum(axis=1)
    near1=near.sum(axis=1)
    minm=np.where(valid,zeros,LOC+1).min(axis=1)
    per_locus=[]
    for locus in range(LOC):
        pp=par[:,locus,:]
        cc=child[:,locus,:]
        per_locus.append({
            "locus_index_one_based":locus+1,
            "mother_offspring_mendelian_mismatch":int(mother_inconsistent[:,locus].sum()),
            "adult_zero_coded_allele_cells":int((pp==0).sum()),
            "offspring_zero_coded_allele_cells":int((cc==0).sum()),
            "adult_observed_allele_diversity":int(len(np.unique(pp))),
            "offspring_observed_allele_diversity":int(len(np.unique(cc))),
            "offspring_mother_locus_fully_called":int(mother_child_locus_informative[:,locus].sum()),
            "offspring_with_at_least_one_zero_mismatch_father_after_omitting_this_locus":
                int(((((exact[:,:,np.arange(LOC)!=locus]==0)&father_locus_informative[:,:,np.arange(LOC)!=locus]).sum(axis=2)==0)&
                     (father_locus_informative[:,:,np.arange(LOC)!=locus].sum(axis=2)>=4)).any(axis=1).sum())})
    by_treat={}
    for label in ("c","de","ne"):
        v=np.array(d["child_treatment"])==label
        q=compatible[v]
        by_treat[label]={
         "offspring":int(v.sum()),
         "fully_mendel_compatible_father_zero":int((q==0).sum()),
         "fully_mendel_compatible_father_one":int((q==1).sum()),
         "fully_mendel_compatible_father_more_than_one":int((q>1).sum()),
         "median_compatible_candidate_count":float(np.median(q)),
         "mean_compatible_candidate_count":float(np.mean(q)),
         "median_le1_mismatch_candidate_count":float(np.median(near1[v])),
        }
    no_dist={"parents_sampled":len(par),"offspring_genotyped":n_off,
       "loci":LOC,"mother_genotype_inconsistent_loci_total":int(mother_inconsistent.sum()),
       "offspring_with_any_mother_genotype_exclusion":int(mother_inconsistent.any(axis=1).sum()),
       "offspring_with_zero_compatible_father":int((compatible==0).sum()),
       "offspring_with_unique_zero_mismatch_father":int((compatible==1).sum()),
       "offspring_with_multiple_zero_mismatch_fathers":int((compatible>1).sum()),
       "all_offspring_median_compatible_fathers":float(np.median(compatible)),
       "all_offspring_max_compatible_fathers":int(compatible.max()),
       "mother_offspring_mismatch_count_per_child":{str(int(k)):int(v) for k,v in collections.Counter(mother_inconsistent.sum(axis=1)).items()},
       "source_zero_coded_alleles_are_treated_as_missing":True,
       "candidate_father_minimum_complete_informative_loci":5,
       "mother_child_observed_loci_median":float(np.median(mother_child_locus_informative.sum(axis=1))),
       "mother_child_observed_loci_min":int(mother_child_locus_informative.sum(axis=1).min()),
       "maternal_compatibility_by_locus":per_locus,
       "mean_minimal_paternal_locus_mismatches":float(minm.mean()),
       "offspring_with_min_2_plus_father_mismatch":int((minm>=2).sum()),
       "by_original_pollination_treatment":by_treat,
       "original_source_parentage_identification_is_genetic_only":True,
       "parent_genotype_is_original_census_no_distance_used":True}
    return exact,no_dist

def family_weights(d,exact,eps):
    # Non-calibrated likelihood, NO spatial kernels or distances.
    # Genotyping error floor eps allows unsampled fathers. In the absence of
    # precise allele errors or allele frequency likelihood, normalized values
    # are *weights*, NOT validated paternity posteriors.
    ls=np.log(exact+eps).sum(axis=2)
    ls[np.arange(len(ls)),d["mother_index"]]=-1e8
    mx=ls.max(axis=1,keepdims=True)
    ww=np.exp(ls-mx);ww/=ww.sum(axis=1,keepdims=True)
    assert np.allclose(ww.sum(axis=1),1)
    fam=collections.defaultdict(list)
    for i,key in enumerate(d["family"]):fam[key].append(i)
    per_group=[]
    for key,ids in fam.items():
        if len(ids)<2:continue
        w=ww[ids]
        kin=w@w.T
        upper=kin[np.triu_indices(len(ids),1)]
        per_group.append({
          "key":key,
          "n_offspring":len(ids),
          "p_same":float(upper.mean()),
          "effective_donors_average":float(1/np.mean(np.sum(w*w,axis=1))),
          "source_mltr_correlatedPaternity":float(d["source_ref"][key]) if key in d["source_ref"] else None,
          "treat":key[1]})
    matched=[z for z in per_group if z["source_mltr_correlatedPaternity"] is not None]
    corr=np.array([z["source_mltr_correlatedPaternity"] for z in matched])
    est=np.array([z["p_same"] for z in matched])
    all_group=set(d["family"])
    by={}
    for t in ("c","de","ne"):
        g=[z for z in per_group if z["treat"]==t]
        gref=[z for z in g if z["source_mltr_correlatedPaternity"] is not None]
        by[t]={
          "families_with_two_or_more_children":len(g),
          "mean_genotype_only_pair_cofather_weight":float(np.mean([v["p_same"] for v in g])),
          "families_in_published_MLTR_subset":len(gref),
          "mean_source_MLTR_correlated_paternity":float(np.mean([v["source_mltr_correlatedPaternity"] for v in gref])) if gref else None,
          "median_children_per_family":float(np.median([v["n_offspring"] for v in g]))
        }
    # Matched, repeated-mother comparison: both pollinator exclusion
    # treatments were applied to the same maternal genet, so maternal
    # isolation and SSR genotype are held constant by design.
    bymom=collections.defaultdict(dict)
    for z in per_group:bymom[z["key"][0]][z["treat"]]=z
    def pair_test(metric):
        dif=[]
        for group in bymom.values():
            if "de" not in group or "ne" not in group:continue
            a,b=group["de"][metric],group["ne"][metric]
            if a is None or b is None:continue
            dif.append(float(b-a))  # night excluded (bee) minus day excluded (moth)
        vals=np.array(dif,float)
        assert len(vals)>=8
        rng=np.random.default_rng(14115+int(round(eps*100000))+len(vals))
        bi=rng.integers(len(vals),size=(9999,len(vals)))
        boot=vals[bi].mean(axis=1)
        swap=rng.choice([-1.,1.],size=(19999,len(vals)))
        null=(swap*vals).mean(axis=1)
        return {
          "paired_maternal_genets":len(vals),
          "night_excluded_bee_minus_day_excluded_moth_mean":float(vals.mean()),
          "paired_mother_boot_95":np.quantile(boot,[.025,.975]).tolist(),
          "paired_signflip_two_sided_p":float((1+(np.abs(null)>=abs(vals.mean())).sum())/(1+len(null))),
          "positive_mother_differences":int((vals>0).sum()),
          "negative_mother_differences":int((vals<0).sum()),
          "not_independent_family_units":True}
    pair={"distance_free_genotype_only":pair_test("p_same"),
          "published_same_family_MLTR":pair_test("source_mltr_correlatedPaternity")}
    return {
       "within_mother_paired_exclusion_test":pair,
       "genotyping_error_floor":eps,
       "method":"un-calibrated prior-free per-parent transmission product with epsilon floor, exclude mother, normalize across all sampled adults",
       "all_source_families":len(all_group),
       "all_source_MLTR_families":len(d["source_ref"]),
       "families_with_pair_information":len(per_group),
       "genotyped_families_with_original_MLTR_match":len(matched),
       "genotyped_families_not_in_original_MLTR_table":len(per_group)-len(matched),
       "families_posterior_correlated_paternity_spearman":float(spearmanr(corr,est).statistic),
       "families_posterior_correlated_paternity_pearson":float(pearsonr(corr,est).statistic),
       "estimated_mean_pair_cofather_weight":float(est.mean()),
       "source_mean_MLTR_correlated_paternity":float(corr.mean()),
       "by_treatment":by
    }

def main():
    d=source()
    exact,quality=primary(d)
    w=[family_weights(d,exact,eps) for eps in EPS_LIST]
    assert d["original_genotype_sha256"]=="8348dd28cbc4532cea22f80676d241cea322527bec35e8ec1995edd812be30e7"
    assert d["original_paternity_sha256"]=="c703d56287b62f7550931e46f54591fa7a872438420bf401814eea3186320eb6"
    return {"status":"independent_Oenothera_source_genotype_only_no_spatial_prior_v57",
      "published_study":"10.1111/mec.14115","Dryad":"10.5061/dryad.p24q3",
      "fixed_Zenodo_record":4942351,
      "genotype_workbook_sha256":d["original_genotype_sha256"],
      "fruit_MLTR_source_sha256":d["original_paternity_sha256"],
      "Mendelian_identity":quality,
      "uncalibrated_genotype_only_weight_sensitivity":w,
      "limitations":{
        "source_has_independent_7_locus_offspring_and_parent_genotypes":True,
        "no_spatial_prior_or_locations_used":True,
        "paternity_sample_source_published_2017_not_new_2026_fielddata":True,
        "candidate_fathers_only_sampled_adults_no_external_pollen":True,
        "zero_encoded_SSR_alleles_are_missing_not_true_size_zero":True,
        "strict_sire_compatibility_requires_at_least_5_informative_loci":True,
        "Mendelian_disagreements_may_mean_genotyping_error_mutation_or_nonmatches":True,
        "soft_error_floors_not_empirical_paternity_posteriors":True,
        "same_source_MLTR_correlatedPaternity_not_independent_same_family_replicate":True,
        "not_focal_Daphne_pseudomezereum_nor_reestimation_of_D_jezoensis":True,
        "no_original_individual_identifiers_or_alleles_exported":True}}
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=OUT/"oenothera_genotype_distance_free_v57.json")
    p=a.parse_args()
    p.output.parent.mkdir(exist_ok=True,parents=True)
    o=main()
    p.output.write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("OENOTHERA_GENOTYPE_ONLY_V57",json.dumps(o,ensure_ascii=False,allow_nan=False),flush=True)
