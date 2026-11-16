#!/usr/bin/env python3
"""Source-first structure audit of Daphne jezoensis floral longevity (Shibata & Kudo 2023).
Only aggregated source schema and experimental units; no genotype/plant identifiers output.
Author data Zenodo 10.5281/zenodo.8253921.
"""
from __future__ import annotations
import hashlib,io,json,urllib.request,zipfile
from pathlib import Path
import pandas as pd,numpy as np

OUT=Path("analysis_results_daphne_longevity_v20")
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
CSV_SHA="140aeae0064ff9e13ff4199794abe7dc533d5ac4f3fa0e929fc815eb0d4df017"
def get(url, cap):
    req=urllib.request.Request(url,headers={"User-Agent":"published-Daphne-floral-longevity-audit-v20","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=75) as res:r=res.read(cap+1)
    if len(r)>cap:raise RuntimeError("Exceeded retrieval cap")
    return r
def main():
    OUT.mkdir(exist_ok=True)
    rec=json.loads(get("https://zenodo.org/api/records/8253921",2_000_000))
    fs=rec.get("files")
    ff=fs.get("entries",{}) if isinstance(fs,dict) else {x.get("key"):x for x in fs}
    url=ff["Data for Figures_Analyses .zip"]["links"].get("content") or ff["Data for Figures_Analyses .zip"]["links"]["self"]
    archive=get(url,1_000_000)
    assert hashlib.sha256(archive).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        matches=[x for x in z.namelist() if Path(x).name=="FlowerLongevity.csv"]
        assert len(matches)==1
        b=z.read(matches[0])
    assert hashlib.sha256(b).hexdigest()==CSV_SHA
    d=pd.read_csv(io.BytesIO(b))
    assert len(d)==408
    assert {"RametID","InfloID","Sex","Treatment","Date","FloweringNo","FlowerNo"}<=set(d.columns)
    result={"status":"source_design_inventory_only",
        "source_article":"10.1002/ajb2.16239",
        "source_data":"10.5281/zenodo.8253921",
        "source_archive_sha":ZIP_SHA,"source_flower_table_sha":CSV_SHA,
        "n_records":len(d),"missing":{str(k):int(v) for k,v in d.isna().sum().items()},
        "n_ramets":int(d.RametID.nunique()),"n_inflo_id":int(d.InfloID.nunique()),
        "full_duplicates":int(d.duplicated().sum()),
        "key_duplicates":int(d.duplicated(["RametID","InfloID","Treatment","Date"]).sum()),
        "columns":list(d.columns)}
    for var in ("FloweringNo","FlowerNo"):
        a=pd.to_numeric(d[var],errors="coerce")
        result[var]={"range":[float(a.min()),float(a.max())],"n_non_numeric":int(a.isna().sum()),
            "n_negative":int((a<0).sum()),"n_positive":int((a>0).sum()),
            "n_zero":int((a==0).sum()),
            "n_above_FlowerNo":int((a>pd.to_numeric(d.FlowerNo,errors="coerce")).sum())}
    dt=pd.to_datetime(d.Date,errors="coerce")
    result["date"]={"n_parsed":int(dt.notna().sum()),
       "n_unique_dates":int(dt.nunique()),
       "minimum_date":str(dt.min())[:10],
       "maximum_date":str(dt.max())[:10],
       "source_example_date_value":str(d.Date.iloc[0])[:30]}
    # Group levels and unit counts. Do not output any RametID or InfloID.
    units=d.groupby(["RametID","InfloID","Sex","Treatment"],dropna=False).agg(
        rows=("Date","size"),n_dates=("Date","nunique"),
        first=("Date","min"),last=("Date","max")).reset_index()
    result["unit_integrity"]={
       "n_unique_ramet_inflo_sex_treatment":len(units),
       "min_census_per_unit":int(units.rows.min()),
       "max_census_per_unit":int(units.rows.max()),
       "census_frequency":{str(k):int(v) for k,v in units.rows.value_counts().items()},
       "n_ramets_with_many_sexes":int((d.groupby("RametID").Sex.nunique()>1).sum())}
    result["groups"]=[]
    for (sex,treat),a in d.groupby(["Sex","Treatment"]):
        unit=a[["RametID","InfloID"]].drop_duplicates()
        result["groups"].append({"sex":str(sex),"treatment":str(treat),
           "n_rows":len(a),"n_ramets":int(a.RametID.nunique()),
           "n_unit_pairs":len(unit),
           "unique_dates":int(a.Date.nunique()),
           "n_obs_FloweringNo_gt0":int((pd.to_numeric(a.FloweringNo)>0).sum()),
           "mean_FloweringNo":round(float(pd.to_numeric(a.FloweringNo).mean()),5),
           "mean_FlowerNo":round(float(pd.to_numeric(a.FlowerNo).mean()),5),
           "median_FloweringNo":round(float(pd.to_numeric(a.FloweringNo).median()),5)})
    # Which pairs have both F bag & F outcross / H emasculate_bag & emasculate_outcross?
    combos=d.groupby(["RametID","Sex"]).Treatment.agg(lambda v:set(v)).reset_index()
    for sex,a,b in [("F","bag","outcross"),("H","emasculate_bag","emasculate_outcross"),
                    ("H","emasculate_bag","self"),("H","emasculate_bag","natural")]:
        x=combos.loc[combos.Sex.eq(sex),"Treatment"]
        result[f"n_{sex}_same_ramet_with_{a}_and_{b}"]=int(sum(a in y and b in y for y in x))
    # Check which inflorescences appear under multiple treatments: potential bad key
    byinf=d.groupby(["RametID","InfloID"]).Treatment.nunique()
    result["n_ramet_inflo_with_multiple_treatments"]=int((byinf>1).sum())
    p=OUT/"daphne_2023_original_flower_longevity_schema_v20.json"
    p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_LONGEVITY_SOURCE_V20",json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
