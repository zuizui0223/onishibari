#!/usr/bin/env python3
"""v0.23: 2015 natural flowering source aggregation, independent year of author study.

Quantify normalized and absolute 35-day open-flower opportunity, by sex class.
This is NOT individual flower lifespan, reproductive seed output, resource
allocation, sex randomized intervention, or independent publication replication.
"""
import hashlib,io,json,urllib.request,zipfile,datetime
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
OUT=Path("analysis_results_daphne_2015_flower_days_v23")
ZIPSHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
CSVSHA="1c3daed5ad7b61af2175ae02194c948e91fb0e3cbfc0ed0ab1c49b5d5e608968"
NBOOT=15000
RNG=np.random.default_rng(20261116)
def get(url,limit=2000000):
    req=urllib.request.Request(url,headers={"User-Agent":"Daphne-published-2015-sex-duration-v23","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as f:r=f.read(limit+1)
    if len(r)>limit:raise RuntimeError("source too large")
    return r
def source():
    j=json.loads(get("https://zenodo.org/api/records/8253921"))
    files=j["files"];es=files.get("entries",{}) if isinstance(files,dict) else {x["key"]:x for x in files}
    zurl=es["Data for Figures_Analyses .zip"]["links"]
    archive=get(zurl.get("content") or zurl["self"])
    assert hashlib.sha256(archive).hexdigest()==ZIPSHA
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        names=[n for n in z.namelist() if Path(n).name=="FlowerLongevity_2015.csv"]
        assert len(names)==1
        content=z.read(names[0])
    assert hashlib.sha256(content).hexdigest()==CSVSHA
    d=pd.read_csv(io.BytesIO(content))
    return d
def safe(x):
    return round(float(x),7)
def month_day(s):
    return datetime.datetime.strptime(str(s).strip(),"%Y/%m/%d")
def analyse(d):
    assert len(d)==905
    assert set(d.Sex)=={"F","H","HF"}
    assert d.RametID.nunique()==181
    assert d.groupby("Sex").RametID.nunique().to_dict()=={"F":76,"H":90,"HF":15}
    assert set(d.Year)=={2015}
    assert not d.duplicated(["RametID","Date"]).any()
    d=d.copy()
    d["dt"]=d.Date.map(month_day)
    start=d.dt.min()
    d["day"]=(d.dt-start).dt.days
    days=sorted(d.day.unique())
    assert days==[0,9,18,27,35],days
    count=pd.to_numeric(d.FloweringNo,errors="coerce")
    total=pd.to_numeric(d.FlowerNo,errors="coerce")
    assert count.notna().all() and total.notna().all()
    assert total.gt(0).all() and count.ge(0).all() and count.le(total).all()
    d["open"]=count
    d["total"]=total
    assert d.groupby("RametID").total.nunique().eq(1).all()
    units=d.groupby(["RametID","Sex"]).size()
    assert len(units)==181 and units.eq(5).all()
    rows=[]
    for (plant,sex),v in d.groupby(["RametID","Sex"]):
        v=v.sort_values("day")
        assert list(v.day)==days
        n=float(v.total.iloc[0])
        opens=v.open.to_numpy(dtype=float)
        frac=opens/n
        rows.append({"id_private":plant,"sex":sex,
          "total_flowers":n,
          "absolute_open_flower_days":float(np.trapezoid(opens,days)),
          "normalized_open_days":float(np.trapezoid(frac,days)),
          "mean_open_fraction":float(np.trapezoid(frac,days))/35,
          "last_observed_positive":int(max(v.day.loc[v.open.gt(0)])) if v.open.gt(0).any() else None})
    u=pd.DataFrame(rows)
    return u,days
def meanboot(a,b):
    x=a.to_numpy(dtype=float);y=b.to_numpy(dtype=float)
    ix=RNG.integers(0,len(x),size=(NBOOT,len(x)))
    iy=RNG.integers(0,len(y),size=(NBOOT,len(y)))
    diff=x[ix].mean(axis=1)-y[iy].mean(axis=1)
    ratio=x[ix].mean(axis=1)/y[iy].mean(axis=1)
    return {"diff_ci95":[safe(v) for v in np.quantile(diff,[.025,.975])],
      "ratio_ci95":[safe(v) for v in np.quantile(ratio,[.025,.975])]}
def group(u):
    out={}
    for sex,v in u.groupby("sex"):
        out[sex]={"n_plant":len(v),
          "mean_total_flower_per_plant":safe(v.total_flowers.mean()),
          "median_total_flower_per_plant":safe(v.total_flowers.median()),
          "mean_normalized_open_days":safe(v.normalized_open_days.mean()),
          "median_normalized_open_days":safe(v.normalized_open_days.median()),
          "mean_absolute_flower_days":safe(v.absolute_open_flower_days.mean()),
          "median_absolute_flower_days":safe(v.absolute_open_flower_days.median()),
          "mean_open_fraction_35_day_window":safe(v.mean_open_fraction.mean()),
          "n_with_last_census_positive":int(v.last_observed_positive.eq(35).sum()),
          "n_with_first_census_positive":None}
    return out
def main():
    OUT.mkdir(exist_ok=True)
    d=source()
    u,days=analyse(d)
    a=u[u.sex=="F"];b=u[u.sex=="H"]
    contrast={}
    for metric in ["normalized_open_days","absolute_open_flower_days","total_flowers"]:
        x=a[metric];y=b[metric]
        ci=meanboot(x,y)
        contrast[metric]={"F_minus_H_mean":safe(x.mean()-y.mean()),
          "F_div_H_mean_ratio":safe(x.mean()/y.mean()),
          **ci}
    # Census source has n=15 "HF" source codes, but the authors excluded them
    # from the original F-versus-H Figure S3. Do not assume their functional sex.
    out={"status":"PASS_original_2015_Daphne_window_exposure_v23",
      "source_doi":"10.1002/ajb2.16239",
      "zenodo_doi":"10.5281/zenodo.8253921",
      "zip_sha256":ZIPSHA,"csv_sha256":CSVSHA,
      "n_source_rows":len(d),"n_individual_ramets":len(u),
      "census_35day_offsets":days,
      "sex_counts":{"F":len(a),"H":len(b),"HF_undefined_excluded_from_FH":int(u.sex.eq("HF").sum())},
      "sex_group_summary":group(u),
      "F_vs_H":contrast,
      "unmeasured_true_flower_lifespan":True,
      "original_author_reported_same_pattern_in_2015":True,
      "separate_census_from_2018_5F_9H_experiment":True,
      "abs_flower_days_are_display_not_paternal_fitness":True,
      "no_Daphne_pseudomezereum_functional_sex_inferred":True}
    file=OUT/"daphne_2015_166_FH_seasonal_open_days_v23.json"
    file.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_2015_NATURAL_FLOWER_DAYS_V23",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
