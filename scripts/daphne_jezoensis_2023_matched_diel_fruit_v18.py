#!/usr/bin/env python3
"""Original 2023 Daphne jezoensis day- versus night-open fruit set, v0.18.

Independent *secondary analysis* of Shibata & Kudo (2023), AJB 110:e16239.
Reuses authors' publicly deposited source from https://zenodo.org/records/8253921.
Confirm source SHA, 208 valid inflorescences, 23 year×female pairs (22 ramet IDs).
Aggregate per plant and treatment; do not mistake 208 inflorescences as independent plants.
Primary descriptive contrast is mean across female plant-year pairs of night minus day
proportion of flowers setting a fruit. Bootstrap unit is the paired plant. Paired
sign flips are only a diagnostic sensitivity, NOT a design-randomization p-value.
No direct flow from visitation to genetic siring, and NO Daphne pseudomezereum data.
"""
from __future__ import annotations
import hashlib, io, json, urllib.request, zipfile, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest

OUT=Path("analysis_results_daphne_jezoensis_paired_diel_v18")
SOURCE="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
SOURCE_CSV_SHA="28f7c9d54cdf92e95815ba70c0308f6a6e9fd116ddf6396efa71f88ceb1cc7de"
RNG=np.random.default_rng(20261116)
NBOOT=12000
GROUPS=("Day open","Night open")

def fetch(url,cap=1000000):
    request=urllib.request.Request(url,headers={"User-Agent":"published-original-Daphne-diel-paired-secondary-v18","Accept":"*/*"})
    with urllib.request.urlopen(request,timeout=80) as res:raw=res.read(cap+1)
    if len(raw)>cap:raise RuntimeError("source bytes beyond cap")
    return raw

def load():
    manifest=json.loads(fetch(SOURCE,2500000))
    file_inventory=manifest["files"]
    entries=file_inventory.get("entries",{}) if isinstance(file_inventory,dict) else {x.get("key"):x for x in file_inventory}
    url=entries["Data for Figures_Analyses .zip"]["links"].get("content") or entries["Data for Figures_Analyses .zip"]["links"].get("self")
    b=fetch(url,1500000)
    assert hashlib.sha256(b).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        matches=[x for x in z.namelist() if x.split("/")[-1]=="Nopporo_daynoon.csv"]
        assert len(matches)==1
        bb=z.read(matches[0])
    assert hashlib.sha256(bb).hexdigest()==SOURCE_CSV_SHA
    d=pd.read_csv(io.BytesIO(bb),encoding="utf-8-sig")
    return d

def clean(d):
    assert len(d)==216,"expected original 216 rows, including 8 missing outcomes"
    assert {"Year","Sex","Treatment","RametID","InfloID","FlowerNo","FruitNo"} <=set(d)
    assert set(d["Treatment"])==set(GROUPS)
    assert set(d["Year"])=={2018,2020}
    assert d["Sex"].eq("F").all()
    assert d.FruitNo.isna().sum()==8
    a=d.dropna(subset=["FlowerNo","FruitNo"]).copy()
    a["FlowerNo"]=pd.to_numeric(a.FlowerNo,errors="raise").astype(int)
    a["FruitNo"]=pd.to_numeric(a.FruitNo,errors="raise").astype(int)
    assert len(a)==208 and a.FlowerNo.ge(1).all()
    assert a.FruitNo.ge(0).all() and a.FruitNo.le(a.FlowerNo).all()
    assert a.groupby(["Year","Treatment"]).size().to_dict()=={
       (2018,"Day open"):40,(2018,"Night open"):40,
       (2020,"Day open"):64,(2020,"Night open"):64}
    assert not a.duplicated(["Year","RametID","InfloID","Treatment"]).any(), "repeated flower head identifier"
    n=a.groupby(["Year","RametID"])["Treatment"].nunique()
    assert len(n)==23 and n.eq(2).all()
    assert a.RametID.nunique()==22
    return a

def aggregate(a):
    p=a.groupby(["Year","RametID","Treatment"],as_index=False).agg(
        fruits=("FruitNo","sum"),flowers=("FlowerNo","sum"),
        inflorescences=("InfloID","size"))
    assert len(p)==46
    p["rate"]=p.fruits/p.flowers
    wide=p.pivot(index=["Year","RametID"],columns="Treatment",values=["fruits","flowers","inflorescences","rate"])
    assert wide.shape[0]==23 and not wide.isna().any().any()
    wide.columns=["_".join(x) for x in wide.columns]
    wide=wide.reset_index()
    wide["delta_night_minus_day"]=wide["rate_Night open"]-wide["rate_Day open"]
    return p,wide

def signflip_exact(values):
    x=np.asarray(values,dtype=float)
    n=len(x)
    assert 0<n<=17
    stat=abs(x.mean())
    ex=0
    for bits in itertools.product((-1,1),repeat=n):
        y=np.asarray(bits)*x
        if abs(y.mean())>=stat-1e-12:ex+=1
    return ex/(2**n)

def bstrat(a):
    # Bootstrap paired maternal individuals rather than inflorescences.
    n=len(a)
    pairs=np.arange(n)
    ids=RNG.integers(0,n,(NBOOT,n))
    d=a["delta_night_minus_day"].to_numpy()
    x=a["fruits_Day open"].to_numpy()
    y=a["fruits_Night open"].to_numpy()
    nx=a["flowers_Day open"].to_numpy()
    ny=a["flowers_Night open"].to_numpy()
    vals=d[ids].mean(axis=1)
    pool=y[ids].sum(axis=1)/ny[ids].sum(axis=1)-x[ids].sum(axis=1)/nx[ids].sum(axis=1)
    return [float(v) for v in np.quantile(vals,[.025,.975])],[float(v) for v in np.quantile(pool,[.025,.975])]

def row_stats(g):
    x=g["delta_night_minus_day"].to_numpy()
    b,bc=bstrat(g)
    day_f=int(g["fruits_Day open"].sum())
    night_f=int(g["fruits_Night open"].sum())
    day_n=int(g["flowers_Day open"].sum())
    night_n=int(g["flowers_Night open"].sum())
    p_day=day_f/day_n
    p_night=night_f/night_n
    npos=int((x>0).sum())
    nneg=int((x<0).sum())
    nzero=int((x==0).sum())
    sig=binomtest(npos,npos+nneg,.5,alternative="two-sided").pvalue if (npos+nneg)>0 else 1
    return {
      "paired_plant_years":len(g),
      "distinct_ramets":int(g.RametID.nunique()),
      "inflorescence_day":int(g["inflorescences_Day open"].sum()),
      "inflorescence_night":int(g["inflorescences_Night open"].sum()),
      "day_fruits":day_f,"day_flowers":day_n,
      "night_fruits":night_f,"night_flowers":night_n,
      "pooled_day_rate":round(p_day,7),"pooled_night_rate":round(p_night,7),
      "pooled_night_minus_day":round(p_night-p_day,7),
      "pooled_diff_paired_plant_boot_CI95":[round(v,7) for v in bc],
      "mean_plant_day_rate":round(float(g["rate_Day open"].mean()),7),
      "mean_plant_night_rate":round(float(g["rate_Night open"].mean()),7),
      "mean_matched_plant_rate_diff":round(float(x.mean()),7),
      "mean_matched_plant_boot_CI95":[round(v,7) for v in b],
      "median_matched_plant_rate_diff":round(float(np.median(x)),7),
      "n_night_higher":npos,"n_day_higher":nneg,"n_tied":nzero,
      "sign_test_two_sided":round(float(sig),7),
      "paired_signflip_two_sided_sensitivity_NOT_randomization_p":round(float(signflip_exact(x)),7) if len(x)<=17 else None,
      "pooled_night_day_risk_ratio_descriptive":round(float(p_night/p_day),7) if p_day else None
    }

def two_years(w):
    stats={str(y):row_stats(g) for y,g in w.groupby("Year")}
    # Genet-resampling yields uncertainty on treatment-by-year contrast.
    rng=np.random.default_rng(20261117)
    all_ids=w.RametID.drop_duplicates().to_numpy()
    # same sample set for both years in contrast: independent sequence
    difference=[]
    for _ in range(12000):
        draw=all_ids[rng.integers(0,len(all_ids),len(all_ids))]
        ct=pd.Series(draw).value_counts()
        x=w.assign(_n=w.RametID.map(ct).fillna(0))
        a=x[(x.Year==2018)&(x._n>0)]
        b=x[(x.Year==2020)&(x._n>0)]
        if len(a) and len(b):
            difference.append(float(np.average(a.delta_night_minus_day,weights=a._n)-
                                    np.average(b.delta_night_minus_day,weights=b._n)))
    stats["year_interaction_diagnostic"]={
      "mean_2018_minus_2020_night_day_diff":round(
         stats["2018"]["mean_matched_plant_rate_diff"]-
         stats["2020"]["mean_matched_plant_rate_diff"],7),
      "unique_ramet_cluster_boot_CI95":[round(float(v),7) for v in np.quantile(difference,[.025,.975])],
      "n_unique_ramets_in_two_years":int(len(all_ids)),
      "shared_ramets_across_both_years":int(w.RametID.value_counts().gt(1).sum()),
      "not_a_design_randomization_p":True
    }
    return stats

def main():
    OUT.mkdir(exist_ok=True)
    raw=load()
    d=clean(raw)
    plant,wide=aggregate(d)
    output={
     "status":"success_original_2023_daphne_jezoensis_matched_maternal_fruit_secondary",
     "source_original_2023_article_doi":"10.1002/ajb2.16239",
     "source_original_2023_zenodo":"10.5281/zenodo.8253921",
     "zip_sha256":ZIP_SHA,"csv_sha256":SOURCE_CSV_SHA,
     "n_source_rows":len(raw),"n_rows_with_valid_mature_fruit":len(d),
     "n_missing_fruit_excluded":int(raw.FruitNo.isna().sum()),
     "n_distinct_female_ramets":int(d.RametID.nunique()),
     "n_plant_year_pairs":int(wide.shape[0]),
     "design":"same female plant receives day- and night-open on DIFFERENT inflorescences; paired maternal effect, not pollen sire effect",
     "source_published_hour_exposure":{"2018":{"day":37.0,"night":32.5},"2020":{"day":95.5,"night":91.0}},
     "source_published_exposure_periods_days":{"2018":{"day":3,"night":4},"2020":{"day":8,"night":8}},
     "by_year":two_years(wide),
     "no_self_incompatibility_or_paternity_inference":True,
     "not_a_direct_flower_biomass_measure":True,
     "not_new_onishibari_data":True,
     "not_a_novel_experiment_source_authors_previously_published_effect":True
    }
    # Note: 216 raw vs 208 eligible precisely equals paper-design original 80+128 inflorescences.
    assert output["n_plant_year_pairs"]==23
    assert output["by_year"]["2018"]["paired_plant_years"]==12
    assert output["by_year"]["2020"]["paired_plant_years"]==11
    f=OUT/"shibata_kudo2023_diel_paired_fruits_v18.json"
    f.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_DIEL_PAIRED_V18",json.dumps(output,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
