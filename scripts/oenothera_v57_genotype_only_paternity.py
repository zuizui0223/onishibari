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
import collections,csv,hashlib,io,json,math
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
    exact=np.zeros((n_off,n_par,LOC),dtype=float)
    mother_inconsistent=np.zeros((n_off,LOC),dtype=bool)
    for ix in range(n_off):
        mom=par[mother[ix]]
        for locus in range(LOC):
            c=child[ix,locus]
            mother_inconsistent[ix,locus]=not np.any(c[:,None]==mom[None,:])
            for f in range(n_par):
                exact[ix,f,locus]=meiosis_probability(mom[locus],par[f,locus],c)
        exact[ix,mother[ix],:]=0
    # Fixed father's true Mendelian match could be absent if parent unsampled
    # or allele calling error; do not impute a father for that offspring.
    zeros=(exact==0).sum(axis=2)
    zeros[np.arange(n_off),mother]=LOC+1
    compatible=(zeros==0).sum(axis=1)
    near1=(zeros<=1).sum(axis=1)
    minm=np.min(zeros,axis=1)
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
          "source_mltr_correlatedPaternity":float(d["source_ref"][key]),
          "treat":key[1]})
    corr=np.array([z["source_mltr_correlatedPaternity"] for z in per_group])
    est=np.array([z["p_same"] for z in per_group])
    all_group=set(d["family"])
    by={}
    for t in ("c","de","ne"):
        g=[z for z in per_group if z["treat"]==t]
        by[t]={
          "families_with_two_or_more_children":len(g),
          "mean_genotype_only_pair_cofather_weight":float(np.mean([v["p_same"] for v in g])),
          "mean_source_MLTR_correlated_paternity":float(np.mean([v["source_mltr_correlatedPaternity"] for v in g])),
          "median_children_per_family":float(np.median([v["n_offspring"] for v in g]))
        }
    return {
       "genotyping_error_floor":eps,
       "method":"un-calibrated prior-free per-parent transmission product with epsilon floor, exclude mother, normalize across all sampled adults",
       "all_source_families":len(all_group),
       "all_source_MLTR_families":len(d["source_ref"]),
       "families_with_pair_information":len(per_group),
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
