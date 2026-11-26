#!/usr/bin/env python3
"""v0.55 treatment-stratified original SSR paternal co-clustering in D. jezoensis.

Corrects v0.54 aggregation across natural, day-open and night-open treatment
in a subset of female mother genets. Uses author's public R codebook:
Figures_Analyses.R joins Sample_DATA.201820_seed.csv DNAID -> Info.
Tests natural-only as the primary independent cohort, plus the three
treatments separately and a pooled mother x treatment statistic in which
all father counts are fixed WITHIN treatment.

Source and new code are all publicly reproducible; DO NOT print original
genet IDs, DNAIDs, genotypes, locations or pairwise parents.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import softmax
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,get,API,ZIP_SHA,stringify
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import build,estimate
from daphne_jezoensis_parentage_clustering_v54 import (
    pair_fraction,equal_mother_pair_fraction,
    independent_null,random_margins_null,mcmc_fixed_margins,rate_samples)

OUT=Path("analysis_results_daphne_treatment_specific_sire_clustering_v55")
SEED=20261124
TREATMENTS=["Natural","Day open","Night open"]
SAMPLE_SHA="739dafbbee07de46ee3c6a10086b754f1e33b697d111a9ec8ee22d5a8a35d9b1"
SAMPLES_NAME="Sample_DATA.201820_seed.csv"
DRAWS_PER_CHAIN=800

def sample_file():
    meta=json.loads(get(API))
    fs=meta["files"]
    entries=fs.get("entries",{}) if isinstance(fs,dict) else {i["key"]:i for i in fs}
    blob=get(entries["Data for Figures_Analyses .zip"]["links"].get("content") or
             entries["Data for Figures_Analyses .zip"]["links"]["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        files=[p for p in z.namelist() if Path(p).name==SAMPLES_NAME]
        assert len(files)==1
        raw=z.read(files[0])
    assert hashlib.sha256(raw).hexdigest()==SAMPLE_SHA
    df=pd.read_csv(io.BytesIO(raw))
    assert {"DNAID","Info"}<=set(df) and not df.DNAID.duplicated().any()
    assert set(df.Info.dropna())>=set(TREATMENTS)
    return df[["DNAID","Info"]].copy(),len(df)

def setup(year,files,sample):
    source=files[f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv"]
    assert len(source)==(257 if year==2018 else 327)
    all_seed=source.merge(sample,on="DNAID",how="left",validate="one_to_one")
    assert len(all_seed)==len(source) and all_seed.Info.notna().all()
    counts=all_seed.Info.value_counts()
    expected={2018:{"Natural":158,"Day open":37,"Night open":62},
              2020:{"Natural":199,"Day open":67,"Night open":61}}[year]
    assert {x:int(counts[x]) for x in TREATMENTS}==expected
    assigned=all_seed.loc[all_seed.Father.notna()&
                          pd.to_numeric(all_seed.Probability,errors="coerce").ge(.95)].copy()
    expected_p095={2018:{"Natural":132,"Day open":32,"Night open":45},
                   2020:{"Natural":163,"Day open":52,"Night open":46}}[year]
    asscounts=assigned.Info.value_counts()
    assert {x:int(asscounts[x]) for x in TREATMENTS}==expected_p095

    a=files["GenetID_geno_phenotype_1820Nopadult_211109.csv"].copy()
    assert len(a)==492
    a["key"]=stringify(a.GenetID2)
    assert a.key.is_unique
    a["x"]=pd.to_numeric(a["East.g"],errors="coerce")
    a["y"]=pd.to_numeric(a["North.g"],errors="coerce")
    a["natural_flowers"]=pd.to_numeric(a[f"NaturalFlowerNo_{year}.g"],errors="coerce")
    a["sex"]=a[f"Sex_{year}.g"].astype(str)
    fathers=a[(a.sex=="H")&a.natural_flowers.gt(0)&a.x.notna()&a.y.notna()].copy().reset_index(drop=True)
    assert len(fathers)==(105 if year==2018 else 71)
    sire_keys=fathers.key.tolist(); sireindex={key:i for i,key in enumerate(sire_keys)}
    base=build(year,files,include_HF=False)
    assert np.array_equal(np.log(fathers.natural_flowers.to_numpy()),base["logflowers"])
    assigned["Mother_key"]=stringify(assigned.Mother)
    assigned["Father_key"]=stringify(assigned.Father)
    good=assigned.loc[assigned.Father_key.isin(sireindex)].copy()
    assert len(good)==(180 if year==2018 else 206)
    # Despite source adult coding such as F,F or F,dead, all sampled mother
    # records have a female component (see original paper design).
    mother_sex=a.set_index("key").sex
    tags=good.Mother_key.map(mother_sex)
    assert tags.notna().all() and tags.str.startswith("F").all()
    adultXY=a.set_index("key")[["x","y"]]
    xy=fathers[["x","y"]].to_numpy(dtype=float)
    assert good.Mother_key.map(adultXY.x).notna().all()
    assert good.Mother_key.map(adultXY.y).notna().all()
    # Critical gate: aggregate source good records by maternal genet, ignoring
    # treatment; exactly match previous v30 parentage arrays, in original
    # mother and father candidate ordering.
    mm=sorted(good.Mother_key.unique())
    mat=np.zeros((len(mm),len(fathers)),dtype=np.int32)
    idx={name:i for i,name in enumerate(mm)}
    for mother,father in zip(good.Mother_key,good.Father_key):
        mat[idx[mother],sireindex[father]]+=1
    assert np.array_equal(mat,base["counts"]), "Treatment join changed original v30 mother/father counts"
    d=np.array([[adultXY.loc[z,"x"],adultXY.loc[z,"y"]] for z in mm])
    ds=np.hypot(d[:,None,0]-xy[None,:,0],d[:,None,1]-xy[None,:,1])
    assert np.allclose(ds,base["distance"],atol=1e-8)
    return {
      "year":year,"source_assigned_full":len(assigned),
      "strict_H_assigned":len(good),
      "source_full_genotyped_by_treatment":expected,
      "source_all_p095_by_treatment":expected_p095,
      "source_strict_H_by_treatment":{t:int((good.Info==t).sum()) for t in TREATMENTS},
      "good":good,
      "adultXY":adultXY,"sire_xy":xy,
      "sirekeys":sire_keys,"sireindex":sireindex,
      "logflowers":base["logflowers"],
      "same_mother_2018_2020_shared_status":"No new independent site; strict H father set year-specific",
      "reference_beta_lambda10":float(estimate(base,10.)["beta_hat"]),
    }

def subset(q,source,treatment):
    rr=source["good"]
    d=rr.loc[rr.Info==treatment].copy()
    groups=[(str(k),g) for k,g in d.groupby("Mother_key",sort=True)]
    assert len(groups)>=4
    mom=np.array([[source["adultXY"].loc[name,"x"],source["adultXY"].loc[name,"y"]]
                  for name,g in groups])
    ix=source["sireindex"]
    c=np.zeros((len(groups),len(source["sirekeys"])),dtype=np.int32)
    for j,(key,g) in enumerate(groups):
        for father in g.Father_key:c[j,ix[father]]+=1
    assert int(c.sum())==len(d)
    distance=np.hypot(mom[:,None,0]-source["sire_xy"][None,:,0],
                      mom[:,None,1]-source["sire_xy"][None,:,1])
    return {"counts":c,"distance":distance,"logflowers":source["logflowers"],
            "mothers":[name for name,g in groups],"treatment":treatment}

def pairStats(counts):
    obs,pairs,same=pair_fraction(counts)
    eq=equal_mother_pair_fraction(counts)
    return {"mothers":len(counts),"offspring":int(counts.sum()),
            "mothers_with_2_or_more":int((counts.sum(axis=1)>=2).sum()),
            "pairs_all":int(pairs),"same_father_pairs":int(same),
            "pooled_pair_observed":obs,"mother_equal_observed":eq}

def observed_and_null(q,year,treatment,lambdas):
    counts=q["counts"]
    ob=pairStats(counts)
    rng=np.random.default_rng(SEED+year*13+TREATMENTS.index(treatment)*133)
    beta=float(estimate(q,10.)["beta_hat"])
    probs=softmax(-q["distance"]/10+beta*q["logflowers"][None,:],axis=1)
    nulA=independent_null(counts,probs,rng,ob["pairs_all"])
    nulB=random_margins_null(counts,10000,rng,ob["pairs_all"])
    result={"treatment":treatment,
        **ob,"model_beta_10m":beta,
        "iid_null_mother_equal":rate_samples(nulA[:,1],ob["mother_equal_observed"]),
        "fixed_margin_uniform_null_mother_equal":rate_samples(nulB[:,1],ob["mother_equal_observed"])}
    sp={}
    for lambda_i,lam in enumerate(lambdas):
        vals=[];qc=[]
        for chain in (0,1):
            seed0=(SEED+year*10+TREATMENTS.index(treatment)*100+lambda_i*10)
            seed=(seed0//10)*20+chain*10
            x,m=mcmc_fixed_margins(counts,q["distance"],float(lam),seed,draws=DRAWS_PER_CHAIN)
            vals.append(x);qc.append(m)
        comb=np.concatenate(vals)
        sp[str(lam)]={
          "pooled_pair_null":rate_samples(comb[:,0],ob["pooled_pair_observed"]),
          "mother_equal_null":rate_samples(comb[:,1],ob["mother_equal_observed"]),
          "chain_mother_equal_means":[float(x[:,1].mean()) for x in vals],
          "chain_mother_equal_mean_diff":float(abs(vals[0][:,1].mean()-vals[1][:,1].mean())),
          "chain_acceptance":[float(m["conditional_swap_acceptance"]) for m in qc],
          "n_draws":len(comb)}
    result["spatial_fixed_margin_nulls"]=sp
    return result

def main():
    files=data()
    samples,nsamples=sample_file()
    result=[]
    for year in (2018,2020):
        source=setup(year,files,samples)
        per=[]
        for name in TREATMENTS:
            q=subset(None,source,name)
            # all treatments use same explicit 10m spatial prior; natural-only
            # sensitivity uses 5/10/20m. Pooled one-third more dependent because
            # day/night from same source female plants.
            per.append(observed_and_null(q,year,name,([5,10,20] if name=="Natural" else [10])))
        pooled_groups=sum(t["mothers_with_2_or_more"] for t in per)
        # Conditional on exact sire totals WITHIN each treatment: combine
        # expected fractions after dividing every maternal-treatment group
        # into its own independent pair study. This is a simple mean of
        # per-treatment null expectations weighted by informative groups.
        combined_obs=sum(t["mother_equal_observed"]*t["mothers_with_2_or_more"] for t in per)/pooled_groups
        combined_exp=sum(t["spatial_fixed_margin_nulls"]["10"]["mother_equal_null"]["null_mean"]*
                          t["mothers_with_2_or_more"] for t in per)/pooled_groups
        result.append({
          "year":year,
          "source_all_genotyped_treatment_counts":source["source_full_genotyped_by_treatment"],
          "source_p095_treatment_counts":source["source_all_p095_by_treatment"],
          "strict_H_father_p095_by_treatment":source["source_strict_H_by_treatment"],
          "strict_H_father_total":source["strict_H_assigned"],
          "full_data_reference_beta":source["reference_beta_lambda10"],
          "by_treatment":per,
          "mother_treatment_stratified_observed_equal":float(combined_obs),
          "mother_treatment_stratified_null_mean_equal":float(combined_exp),
          "n_informative_mother_treatment_groups":pooled_groups,
          "note":"Null father totals constrained within each original pollination treatment; natural and day/night arms differ in source experimental design"
          })
    return {
      "status":"source_original_Daphne_treatment_specific_full_sibling_structure_v55",
      "doi_paper":"10.1002/ajb2.16239","doi_original":"10.5281/zenodo.8253921",
      "zip_sha256":ZIP_SHA,"sample_treatment_csv_sha256":SAMPLE_SHA,
      "sample_metadata_table_rows":nsamples,
      "years":result,
      "limits":{
         "original_parentage_source_uses_spatial_prior":True,
         "source_one_seed_per_fruit":True,
         "day_night_treatments_share_mother_genets":True,
         "natural_only_is_primary_unmixed_exposure_analysis":True,
         "natural_and_exposure_treatments_cannot_be_merged_for_behavior_causation":True,
         "maternal_genet_sample_sizes_unequal_and_some_very_small":True,
         "father_total_is_conditioned_on_observed_success_within_each_treatment":True,
         "individual_pollen_transfer_or_visitor_identity_not_observed":True,
         "no_maternal_ID_sire_ID_genotype_or_site_coordinates_exported":True,
         "focal_Daphne_pseudomezereum_not_measured":True}}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=OUT/"daphne_treatment_specific_maternal_sire_clustering_v55.json")
    a=p.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    res=main()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(res,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("DAPHNE_TREATMENT_STRATIFIED_V55",json.dumps(res,ensure_ascii=False,allow_nan=False),flush=True)
