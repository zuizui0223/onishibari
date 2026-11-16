#!/usr/bin/env python3
"""v0.27: exploratory natural female fruit set versus local potential pollen donors.

Original Shibata & Kudo 2023, DOI 10.1002/ajb2.16239, Zenodo 8253921.
Uses validated within-source RametID -> Ramet match (same sex + same flowers +
same fruits) before assigning adult genotype coordinates. Parentage study
establishes a local pollen-dispersal context but this does NOT link 2015
phenology to 2018/2020 individual fruit observations.

Main descriptive estimand is within-year, plant-equal Spearman correlation
of fruit set with H-only 10 m flowering abundance. HF and wider radii as
sensitivity. NO individual coordinates, source IDs, genotypes or rows output.
Not a causal pollination-limitation test; genetic mate compatibility,
flower-level synchrony, and female resource limits remain unobserved.
"""
from __future__ import annotations
import hashlib,io,json,urllib.request,zipfile,math
from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import spearmanr
OUT=Path("analysis_results_daphne_2018_2020_local_mates_v27")
RECORD="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
HASHES={
 "Nopporo_control.csv":"8b36a936ac7880d6fd7abb9928cd1b0677ffc75e80f1abe6a967c0e67118b07b",
 "GenetID_geno_phenotype_1820Nopadult_211109.csv":"97790d8b3467d2f9e860f97607aa60525a17ddf1b97b0645fdb1b179ce99135e"}
RNG=np.random.default_rng(20261116)
BS=5500

def req(url,cap=2000000):
    r=urllib.request.Request(url,headers={"User-Agent":"original-Daphne-local-potential-mate-v27","Accept":"*/*"})
    with urllib.request.urlopen(r,timeout=80) as f:b=f.read(cap+1)
    if len(b)>cap:raise ValueError("file exceeds cap")
    return b
def source():
    rec=json.loads(req(RECORD))
    fs=rec["files"];es=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    link=es["Data for Figures_Analyses .zip"]["links"]
    blob=req(link.get("content") or link["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    data={}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for name in z.namelist():
            base=Path(name).name
            if base in HASHES:
                b=z.read(name)
                assert hashlib.sha256(b).hexdigest()==HASHES[base]
                data[base]=pd.read_csv(io.BytesIO(b))
    assert set(data)==set(HASHES)
    return data

def normalized(x):
    return x.astype(str).str.strip().str.replace(r"\.0$","",regex=True)
def safe(v):
    return None if not np.isfinite(v) else round(float(v),7)
def spearman_safe(x,y):
    if len(x)<5 or len(np.unique(x))<2 or len(np.unique(y))<2:return None
    return safe(spearmanr(x,y).statistic)

def year_cohort(c,adult,year):
    sexcol=f"Sex_{year}.g"
    flcol=f"NaturalFlowerNo_{year}.g"
    frcol=f"NaturalFruitNo_{year}.g"
    for col in [sexcol,flcol,frcol,"East.g","North.g","Ramet","GenetID2"]:
        assert col in adult
    table=adult.copy()
    assert table.Ramet.notna().all() and not normalized(table.Ramet).duplicated().any()
    assert not table.GenetID2.duplicated().any()
    table["srcID"]=normalized(table.Ramet)
    c=c[(c.Year==year)&c.Treatment.eq("Control")&c.Sex.eq("F")].copy()
    c["srcID"]=normalized(c.RametID)
    c["source_FlowerNo"]=pd.to_numeric(c.FlowerNo,errors="coerce")
    c["source_FruitNo"]=pd.to_numeric(c.FruitNo,errors="coerce")
    merge=c.merge(table,on="srcID",how="left",validate="m:1",suffixes=("","_adult"))
    matched=merge[merge.GenetID2.notna()].copy()
    originalsex=matched.Sex.astype(str)
    geno_sex=matched[sexcol].astype(str)
    adult_fl=pd.to_numeric(matched[flcol],errors="coerce")
    adult_fr=pd.to_numeric(matched[frcol],errors="coerce")
    validkey=(originalsex.eq("F")&geno_sex.eq("F")&
              matched.source_FlowerNo.eq(adult_fl)&matched.source_FruitNo.eq(adult_fr))
    # Some source control rows have missing count values (known 2020 issue).
    # Keep track of those mismatches and EXCLUDE; never force equality or zero-fill.
    n_id_but_no_count_concordance=int((~validkey).sum())
    if year==2018: assert n_id_but_no_count_concordance==0
    else: assert n_id_but_no_count_concordance<=2
    assert originalsex.eq(geno_sex).all(), "sex label mismatch on a shared RametID"
    matched=matched[validkey].copy()
    matched["mxf"]=pd.to_numeric(matched["East.g"],errors="coerce")
    matched["myf"]=pd.to_numeric(matched["North.g"],errors="coerce")
    validout=(matched.source_FlowerNo.gt(0)&
              matched.source_FruitNo.ge(0)&
              matched.source_FruitNo.le(matched.source_FlowerNo)&
              matched.mxf.notna()&matched.myf.notna())
    good=matched[validout].copy()
    good["fruit_rate"]=good.source_FruitNo/good.source_FlowerNo
    good["fruit_any"]=(good.source_FruitNo>0).astype(int)
    good["plotX"]=good.mxf.astype(float)
    good["plotY"]=good.myf.astype(float)
    # Adult published sex code allows ambiguous composite labels (e.g. H,H);
    # only one-letter H is considered a verified candidate in primary analysis.
    # This is a functional pollen source PROXY, not actual verified donor.
    d=table.copy()
    d["x"]=pd.to_numeric(d["East.g"],errors="coerce")
    d["y"]=pd.to_numeric(d["North.g"],errors="coerce")
    d["flowers"]=pd.to_numeric(d[flcol],errors="coerce")
    d["sex"]=d[sexcol].astype(str)
    valid=d.x.notna()&d.y.notna()&d.flowers.gt(0)
    sets={
      "H_only":d[valid & d.sex.eq("H")],
      "H_plus_HF":d[valid & d.sex.isin(["H","HF"])],
      "all_male_labeled_composite":d[valid & d.sex.str.contains(r"(^|,)H(?:,|$)|HF",regex=True)]
    }
    candidate_counts={}
    for k,pop in sets.items():
        candidate_counts[k]=int(len(pop))
        assert not pop.GenetID2.duplicated().any()
        loc=pop[["x","y"]].to_numpy(float)
        fem=good[["plotX","plotY"]].to_numpy(float)
        if len(loc)==0 or len(fem)==0:raise ValueError("No donor candidates or female recipients")
        dist=np.sqrt(np.square(fem[:,None,0]-loc[None,:,0])+
                     np.square(fem[:,None,1]-loc[None,:,1]))
        for radius in [5,10,20]:
            good[f"{k}_nwithin{radius}"]=(dist<radius).sum(axis=1)
            good[f"{k}_flowerswithin{radius}"]=(np.where(dist<radius,pop.flowers.to_numpy(float)[None,:],0)).sum(axis=1)
        good[f"{k}_exp10_flower"]=np.sum(pop.flowers.to_numpy(float)[None,:]*np.exp(-dist/10),axis=1)
        good[f"{k}_nearest_m"]=dist.min(axis=1)
    result={
      "year":year,
      "n_control_F_source":len(c),
      "n_matched_F_to_source_adult_with_validated_counts":len(matched),
      "n_id_only_matches_excluded_due_to_nonconcordance":n_id_but_no_count_concordance,
      "n_valid_matched_F":len(good),
      "n_nonmatching_F":len(c)-len(matched),
      "n_matched_F_excluded_bad_outcome":len(matched)-len(good),
      "n_candidate_flowering_adult_genets":candidate_counts,
      "female_fruit_any_count":int(good.fruit_any.sum()),
      "female_fruit_rate_pooled":safe(good.source_FruitNo.sum()/good.source_FlowerNo.sum()),
      "female_fruit_rate_plant_mean":safe(good.fruit_rate.mean()),
      "female_fruit_rate_plant_median":safe(good.fruit_rate.median()),
      "source_female_flower_sum":int(good.source_FlowerNo.sum()),
      "source_female_fruit_sum":int(good.source_FruitNo.sum()),
      "local_neighbor_summary":{},
      "association":{}}
    for k in sets:
        for radius in (5,10,20):
            x=good[f"{k}_flowerswithin{radius}"].to_numpy(float)
            key=f"{k}_flowerswithin{radius}"
            result["local_neighbor_summary"][key]={
                "n_zero_candidate_flowers":int((x==0).sum()),
                "median_candidate_flowers":safe(np.median(x)),
                "q75_candidate_flowers":safe(np.quantile(x,.75)),
                "max_candidate_flowers":safe(np.max(x))}
            result["association"][key]={
                "spearman_fruit_rate":spearman_safe(x,good.fruit_rate.to_numpy()),
                "spearman_any_fruit":spearman_safe(x,good.fruit_any.to_numpy())}
    primary="H_only_flowerswithin10"
    x=good[primary].to_numpy(float)
    if len(np.unique(x))>1:
        # Quantiles based on ranks avoid empty bins where many 0 donor counts.
        rank=pd.Series(x).rank(method="average",pct=True).to_numpy()
        low=good[rank<=.5];high=good[rank>.5]
        result["primary_median_rank_strata"]={
         "n_low":len(low),"n_high":len(high),
         "low_any_fruit_pct":safe(100*low.fruit_any.mean()) if len(low) else None,
         "high_any_fruit_pct":safe(100*high.fruit_any.mean()) if len(high) else None,
         "low_mean_fruit_rate":safe(low.fruit_rate.mean()) if len(low) else None,
         "high_mean_fruit_rate":safe(high.fruit_rate.mean()) if len(high) else None}
    # Return internal arrays for cluster bootstrap, no raw coordinates or IDs output.
    return good,result

def cluster_results(a,b):
    combined=pd.concat([a,b],ignore_index=True)
    source_ids=combined.GenetID2.astype(str)
    unique=source_ids.drop_duplicates().to_numpy()
    by_year={2018:a,2020:b}
    estimates=[]
    for yy,g in by_year.items():
        estimates.append(spearman_safe(g.H_only_flowerswithin10.to_numpy(float),g.fruit_rate.to_numpy(float)))
    # Weighted mean of within-year Spearman, descriptive, no model-based variance.
    point=sum((len(g)-3)*v for g,v in zip((a,b),estimates))/(len(a)+len(b)-6) if all(v is not None for v in estimates) else None
    boots=[]
    for _ in range(BS):
        chosen=unique[RNG.integers(0,len(unique),len(unique))]
        count=pd.Series(chosen).value_counts()
        stats=[]
        for yy,g in by_year.items():
            freq=g.GenetID2.astype(str).map(count).fillna(0).astype(int)
            gg=g.loc[g.index.repeat(freq)].copy()
            rho=spearman_safe(gg.H_only_flowerswithin10.to_numpy(float),gg.fruit_rate.to_numpy(float))
            if rho is None:break
            stats.append((len(gg)-3,rho))
        if len(stats)==2 and sum(w for w,_ in stats)>0:
            boots.append(sum(w*r for w,r in stats)/sum(w for w,_ in stats))
    result={"primary":"Spearman r of individual female fruit rate versus strict-H 10m flower opportunity, stratified by year",
       "2018_spearman":estimates[0],"2020_spearman":estimates[1],
       "n_unique_identified_genets":int(len(unique)),
       "n_genets_observed_both_years":int((source_ids.value_counts()>1).sum()),
       "year_stratified_plant_equal_weighted_r":safe(point) if point is not None else None,
       "genet_cluster_boot_95%":[safe(v) for v in np.quantile(boots,[.025,.975])] if boots else None,
       "n_boot_valid":len(boots),
       "not_spatial_block_adjusted":True,
       "exploratory_not_causal":True}
    return result

def main():
    OUT.mkdir(exist_ok=True)
    ff=source()
    ctrl=ff["Nopporo_control.csv"];ad=ff["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    a,aa=year_cohort(ctrl,ad,2018)
    b,bb=year_cohort(ctrl,ad,2020)
    primary=cluster_results(a,b)
    output={"status":"original_source_local_mate_fecundity_descriptive_v27",
       "original_article":"10.1002/ajb2.16239",
       "author_source":"10.5281/zenodo.8253921",
       "source_zip_sha256":ZIP_SHA,
       "source_csv_sha256":HASHES,
       "year_profiles":[aa,bb],
       "primary_two_year_sensitivity":primary,
       "limitations":{
         "sex_code_F_H_function_not_verified_here":True,
         "HF_and_composite_donor_codes_are_hypothetical_sensitivities":True,
         "local_abundance_is_not_compatible_pollen_arrival":True,
         "2015_temporal_phenology_not_linked_to_2018_2020_female_outcomes":True,
         "maternal_fruit_set_can_be_resource_limited":True,
         "spatial_autocorrelation_can_inflate_nominal_precision":True,
         "no_unobserved_recipient_genets_imputed":True,
         "not_D_pseudomezereum_field_measurement":True}}
    file=OUT/"daphne_2018_2020_mate_10m_female_fruit_v27.json"
    file.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_LOCAL_MATES_FRUIT_V27",json.dumps(output,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
