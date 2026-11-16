#!/usr/bin/env python3
"""v0.29 original Daphne 2018/2020: HF father distances and 5m mate sensitivity.

Paternal evidence verifies some HF-coded genets can sire offspring; selective
inclusion of *observed* fathers is outcome-conditioned and NOT a causal
predictor for maternal fruit set. Source 2018/2020 natural female fruit and
genetic parentage are derived from Shibata & Kudo 2023 Zenodo 8253921.
Only aggregate counts/intervals; no source individuals or positions disclosed.
"""
from __future__ import annotations
import json
import numpy as np,pandas as pd
from pathlib import Path
from daphne_jezoensis_2018_2020_local_mate_fruit_v27 import source,year_cohort,spearman_safe,RNG as PREV_RNG
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,stringify,cats
OUT=Path("analysis_results_daphne_HF_true_nearby_v29")
RNG=np.random.default_rng(20261116)
BOOT=8000
PFILES={2018:"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
        2020:"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"}
def safe(x):
    return None if not np.isfinite(x) else round(float(x),7)
def ci(arr):
    return [safe(x) for x in np.quantile(np.asarray(arr,dtype=float),[.025,.975])]
def parse_parentage(parentage,adult,year):
    sub=parentage[PFILES[year]]
    q=sub[sub.Father.notna()&pd.to_numeric(sub.Probability,errors="coerce").ge(.95)].copy()
    assert len(q)==(209 if year==2018 else 261)
    base=adult[["GenetID2","East.g","North.g",f"Sex_{year}.g"]].copy()
    base["key"]=stringify(base.GenetID2)
    assert not base.key.duplicated().any()
    mum=base.rename(columns={"key":"mum","East.g":"mx","North.g":"my"})
    sire=base.rename(columns={"key":"sire","East.g":"fx","North.g":"fy",
                              f"Sex_{year}.g":"sire_sex"})
    q["mum"]=stringify(q.Mother)
    q["sire"]=stringify(q.Father)
    v=q.merge(mum[["mum","mx","my"]],on="mum",how="left",validate="m:1")
    v=v.merge(sire[["sire","fx","fy","sire_sex"]],on="sire",how="left",validate="m:1")
    assert v[["mx","my","fx","fy"]].notna().all().all()
    v["dist"]=np.hypot(
        pd.to_numeric(v.mx)-pd.to_numeric(v.fx),
        pd.to_numeric(v.my)-pd.to_numeric(v.fy))
    v["sex_category"]=v.sire_sex.map(cats)
    assert len(v)==len(q)
    out={}
    for cat in ("H_exact","HF_exact","male_mixed_code","other_unknown"):
        z=v[v.sex_category==cat]
        if len(z)==0:
            out[cat]={"assigned_seeds":0,"n_fathers":0}
            continue
        out[cat]={
          "assigned_seeds":len(z),
          "distinct_father_genets":int(z.sire.nunique()),
          "distinct_mother_genets":int(z.mum.nunique()),
          "n_dist_lt5":int(z.dist.lt(5).sum()),
          "n_dist_lt10":int(z.dist.lt(10).sum()),
          "fraction_lt5":safe(z.dist.lt(5).mean()),
          "fraction_lt10":safe(z.dist.lt(10).mean()),
          "mean_distance_m":safe(z.dist.mean()),
          "median_distance_m":safe(z.dist.median()),
          "max_distance_m":safe(z.dist.max())}
    verified=set(v[v.sex_category=="HF_exact"].sire.astype(str))
    assert len(verified)==(1 if year==2018 else 3)
    return out,verified

def per_year(yr,control,adult,parentage):
    fam,qc=year_cohort(control,adult,yr)
    fathers,known_HF=parse_parentage(parentage,adult,yr)
    colsex=f"Sex_{yr}.g"
    colflo=f"NaturalFlowerNo_{yr}.g"
    q=adult.copy()
    q["key"]=stringify(q.GenetID2)
    q["fl"]=pd.to_numeric(q[colflo],errors="coerce")
    q["x"]=pd.to_numeric(q["East.g"],errors="coerce")
    q["y"]=pd.to_numeric(q["North.g"],errors="coerce")
    confirmed=q[(q[colsex]=="HF")&q.key.isin(known_HF)&q.fl.gt(0)].copy()
    assert len(confirmed)==len(known_HF)
    dist=np.hypot(
      fam.plotX.to_numpy(float)[:,None]-confirmed.x.to_numpy(float)[None,:],
      fam.plotY.to_numpy(float)[:,None]-confirmed.y.to_numpy(float)[None,:])
    srcH=fam.H_only_flowerswithin5.to_numpy(float)
    srcHall=fam.H_plus_HF_flowerswithin5.to_numpy(float)
    srcHconfirmed=srcH+(np.where(dist<5,confirmed.fl.to_numpy(float)[None,:],0)).sum(axis=1)
    target=fam.fruit_rate.to_numpy(float)
    rs=[spearman_safe(x,target) for x in (srcH,srcHconfirmed,srcHall)]
    samples=RNG.integers(0,len(fam),size=(BOOT,len(fam)))
    delta=[]
    for ix in samples:
        x0=spearman_safe(srcH[ix],target[ix])
        x1=spearman_safe(srcHconfirmed[ix],target[ix])
        if x0 is not None and x1 is not None:delta.append(x1-x0)
    # Grid block mean analysis at 20m/40m, two fixed origins each.
    blocks=[]
    for size in (20,40):
        for shift in (0,size/2):
            grid=np.floor((fam[["plotX","plotY"]].to_numpy(float)+shift)/size).astype(int)
            t=pd.DataFrame({"gx":grid[:,0],"gy":grid[:,1],
              "local_H":srcH,"local_H_HF_all":srcHall,"fruit":target})
            means=t.groupby(["gx","gy"],as_index=False).mean(numeric_only=True)
            n_blocks=len(means)
            blocks.append({"grid_m":size,"shift_m":shift,
               "n_blocks_with_females":n_blocks,
               "rho_block_mean_H_only":spearman_safe(means.local_H,means.fruit),
               "rho_block_mean_H_plus_all_HF":spearman_safe(means.local_H_HF_all,means.fruit)})
    return {
       "year":yr,
       "n_source_F_matched_plant":len(fam),
       "n_true_HF_sire_genets":len(known_HF),
       "father_distances_by_code":fathers,
       "spearman_F_fruit_rate_vs_5m_flowers":{
          "H_only":rs[0],
          "H_plus_only_HF_with_confirmed_paternity":rs[1],
          "H_plus_all_HF":rs[2],
          "confirmed_HF_increment_delta_rho":safe(rs[1]-rs[0]),
          "confirmed_HF_increment_plant_boot_95%":ci(delta)},
       "grid_mean_sensitivity":blocks,
       "limitations":{
          "observed_sire_filter_is_outcome_conditioned":True,
          "source_mother_fruit_rows_not_same_as_all_genotyped_offspring":True,
          "mother_level_5m_fecundity_test_is_observational":True,
          "paternity_father_distances_have_nonnindependent_seed_siblings":True,
          "no_D_pseudomezereum_data":True
       }}
def main():
    OUT.mkdir(exist_ok=True)
    f=source()
    p=data()
    control=f["Nopporo_control.csv"]
    a=f["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    j={"status":"Daphne_2018_2020_HF_true_father_distance_sensitivity_v29",
      "source_paper":"10.1002/ajb2.16239",
      "source_zenodo":"10.5281/zenodo.8253921",
      "by_year":[per_year(y,control,a,p) for y in [2018,2020]],
      "not_a_same_flower_mediation_estimate":True,
      "not_independent_research_replication":True,
      "no_individual_ID_or_location_released":True}
    (OUT/"daphne_2018_2020_HF_true_paternity_5m_sensitivity_v29.json").write_text(json.dumps(j,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_CONFIRMED_HF_NEARBY_V29",json.dumps(j,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
