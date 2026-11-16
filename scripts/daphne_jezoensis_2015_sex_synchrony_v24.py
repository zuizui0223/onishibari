#!/usr/bin/env python3
"""v0.24: original-source 2015 Daphne jezoensis H-F flowering synchrony audit.

Goal: distinguish female temporal display from temporal availability of pollen donors.
Author study: Shibata & Kudo (2023), DOI 10.1002/ajb2.16239.
Source: Zenodo 8253921, FlowerLongevity_2015.csv (5 censuses, 181 ramets).
H pollen production and F stigma receptivity are only PROXIES from open-flower counts.
No local neighbors, compatibility, insects, paternal outcomes, or fitness inferred.
Reuses original author sample; not an independent replication.
"""
from __future__ import annotations
import hashlib,io,json,urllib.request,zipfile,datetime
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

OUT=Path("analysis_results_daphne_2015_sex_synchrony_v24")
RECORD="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
CSV_SHA="1c3daed5ad7b61af2175ae02194c948e91fb0e3cbfc0ed0ab1c49b5d5e608968"
DAYS=np.array([0.,9.,18.,27.,35.])
LABELS=["Apr 9","Apr 18","Apr 27","May 6","May 14"]
BOOT_N=12000
RNG=np.random.default_rng(20261116)

def get(url,cap=2_000_000):
    req=urllib.request.Request(url,headers={"User-Agent":"onishibari-source-audited-sex-synchrony-v24","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as h:
        b=h.read(cap+1)
    if len(b)>cap:raise RuntimeError("source file exceeded cap")
    return b

def load():
    j=json.loads(get(RECORD))
    fs=j["files"];es=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    zurl=es["Data for Figures_Analyses .zip"]["links"]
    blob=get(zurl.get("content") or zurl["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        paths=[p for p in z.namelist() if Path(p).name=="FlowerLongevity_2015.csv"]
        assert len(paths)==1
        b=z.read(paths[0])
    assert hashlib.sha256(b).hexdigest()==CSV_SHA
    d=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig")
    assert len(d)==905 and d.RametID.nunique()==181
    assert d.groupby("Sex").RametID.nunique().to_dict()=={"F":76,"H":90,"HF":15}
    assert set(d.Year)=={2015}
    assert not d.duplicated(["RametID","Date"]).any()
    d["date"]=pd.to_datetime(d["Date"],format="%Y/%m/%d",errors="raise")
    assert sorted(d.date.dt.strftime("%m-%d").unique())==["04-09","04-18","04-27","05-06","05-14"]
    d["day"]=(d.date-d.date.min()).dt.days.astype(float)
    assert np.array_equal(np.sort(d.day.unique()),DAYS)
    d["open"]=pd.to_numeric(d.FloweringNo,errors="coerce")
    d["total"]=pd.to_numeric(d.FlowerNo,errors="coerce")
    assert d[["open","total"]].notna().all().all()
    assert d.open.ge(0).all() and d.total.gt(0).all() and d.open.le(d.total).all()
    assert d.groupby("RametID").total.nunique().eq(1).all()
    return d

def time_stats(f,h):
    # f, h: temporal sums of observed open flower numbers, each sex group.
    nf=np.trapezoid(f,DAYS,axis=-1)
    nh=np.trapezoid(h,DAYS,axis=-1)
    cenf=np.trapezoid(f*DAYS,DAYS,axis=-1)/nf
    cenh=np.trapezoid(h*DAYS,DAYS,axis=-1)/nh
    latef=0.5*(f[...,3]+f[...,4])*8/nf
    lateh=0.5*(h[...,3]+h[...,4])*8/nh
    # Temporal overlap coefficient of day-normalized display densities:
    # polygonal interpolation onto 0.25-day intervals; only measured 35-day window.
    dense=np.arange(0,35.00001,.25)
    # interpolation vectorization for potentially batched data
    pf=np.array([np.interp(dense,DAYS,vec)/den for vec,den in zip(np.atleast_2d(f),np.atleast_1d(nf))])
    ph=np.array([np.interp(dense,DAYS,vec)/den for vec,den in zip(np.atleast_2d(h),np.atleast_1d(nh))])
    overlap=np.trapezoid(np.minimum(pf,ph),dense,axis=1)
    if np.ndim(f)==1:
        return {"F_centroid":float(cenf),"H_centroid":float(cenh),
                "F_minus_H_centroid":float(cenf-cenh),
                "F_late_share":float(latef),"H_late_share":float(lateh),
                "time_overlap_coefficient":float(overlap[0])}
    return {"F_centroid":cenf,"H_centroid":cenh,
            "F_minus_H_centroid":cenf-cenh,"F_late_share":latef,
            "H_late_share":lateh,"time_overlap_coefficient":overlap}

def summary(d):
    groups={}
    for sex in ("F","H"):
        a=d[d.Sex==sex]
        wide=a.pivot(index="RametID",columns="day",values="open")[list(DAYS)]
        assert not wide.isna().any().any()
        groups[sex]=wide.to_numpy(float)
    f=groups["F"];h=groups["H"]
    # Temporal sums for the originally sampled populations: not all flowering plants in area.
    fs=f.sum(axis=0);hs=h.sum(axis=0)
    ffrac=f.mean(axis=0); hfrac=h.mean(axis=0)
    temporal=time_stats(fs,hs)
    details=[]
    for i,day in enumerate(DAYS):
        details.append({"label":LABELS[i],"day_since_Apr9":int(day),
          "female_active_plants":int((f[:,i]>0).sum()),
          "hermaphrodite_active_plants":int((h[:,i]>0).sum()),
          "female_open_flowers_in_sample":int(fs[i]),
          "hermaphrodite_open_flowers_in_sample":int(hs[i]),
          "female_open_flowers_per_sampled_plant":float(ffrac[i]),
          "hermaphrodite_open_flowers_per_sampled_plant":float(hfrac[i]),
          "source_sample_H_to_F_open_flower_ratio":float(hs[i]/fs[i]) if fs[i]>0 else None,
          "sex_balanced_H_to_F_open_flowers_per_plant_ratio":float(hfrac[i]/ffrac[i]) if ffrac[i]>0 else None})
    # Source sampled H abundance and F abundance might not be spatially representative;
    # no conclusion about genetically compatible, nearby donors is drawn.
    source_fisher=fisher_exact([[int((f[:,-1]>0).sum()),int((f[:,-1]==0).sum())],
                               [int((h[:,-1]>0).sum()),int((h[:,-1]==0).sum())]])
    # independent-plant stratified non-parametric bootstrap
    fi=RNG.integers(0,len(f),size=(BOOT_N,len(f)))
    hi=RNG.integers(0,len(h),size=(BOOT_N,len(h)))
    fboot=f[fi].sum(axis=1)
    hboot=h[hi].sum(axis=1)
    stat_boot=time_stats(fboot,hboot)
    ci={}
    for k,v in stat_boot.items():
        ci[k]=[float(x) for x in np.quantile(v,[.025,.975])]
    out={"status":"source_rechecked_2015_FH_display_temporal_overlap_v24",
       "source_article":"10.1002/ajb2.16239",
       "source_data":"10.5281/zenodo.8253921",
       "ZIP_SHA256":ZIP_SHA,"CSV_SHA256":CSV_SHA,
       "n_plant_F":len(f),"n_plant_H":len(h),"n_unclassified_HF_excluded":15,
       "date_observation_window":LABELS,"census_day_offsets":[int(x) for x in DAYS],
       "census":details,
       "temporal_metrics":temporal,
       "plant_bootstrap_CI95":ci,
       "last_day_active_rate_F":float((f[:,-1]>0).mean()),
       "last_day_active_rate_H":float((h[:,-1]>0).mean()),
       "last_day_F_vs_H_active_fisher_odds_ratio":float(source_fisher.statistic),
       "last_day_F_vs_H_active_fisher_p":float(source_fisher.pvalue),
       "limits":{
         "same_original_2023_study_data":True,
         "pollen_donor_proxy_from_open_H_flowers_not_verified_pollen_release":True,
         "female_receipt_proxy_from_open_F_flowers_not_verified_stigma_receptivity":True,
         "no_local_male_density_or_compatible_mate_distance":True,
         "no_time_of_day_visits_or_mature_fruit_in_2015_census":True,
         "no_causal_claim_that_late_female_flowers_are_unfertilized":True,
         "no_focal_D_pseudomezereum_data":True,
         "five_census_trapezoid_estimation_only":True,
         "2015_cohort_not_independent_research_replication":True
       }}
    return out

def main():
    OUT.mkdir(exist_ok=True)
    result=summary(load())
    assert result["n_plant_F"]==76 and result["n_plant_H"]==90
    assert result["census"][-1]["female_active_plants"]==61
    assert result["census"][-1]["hermaphrodite_active_plants"]==16
    assert result["temporal_metrics"]["time_overlap_coefficient"]>=0
    assert result["temporal_metrics"]["time_overlap_coefficient"]<=1
    path=OUT/"daphne_2015_sex_open_bloom_timing_overlap_v24.json"
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_FH_SYNCHRONY_V24",json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
