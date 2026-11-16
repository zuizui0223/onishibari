#!/usr/bin/env python3
"""v0.22: source audit of 2015 natural flowering censuses in Daphne jezoensis.

A reproducibility check on original authors' prior observation, NOT an
independent population or causal longevity/fitness experiment.
"""
import hashlib,io,json,urllib.request,zipfile,datetime
from pathlib import Path
import pandas as pd,numpy as np
OUT=Path("analysis_results_daphne_2015_phenology_v22")
ZIPSHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
CSVSHA="1c3daed5ad7b61af2175ae02194c948e91fb0e3cbfc0ed0ab1c49b5d5e608968"
def request(url,limit=1500000):
    req=urllib.request.Request(url,headers={"User-Agent":"daphne-2015-original-phenology-source-audit-v22","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as res:b=res.read(limit+1)
    if len(b)>limit:raise RuntimeError("limit exceeded")
    return b
def main():
    OUT.mkdir(exist_ok=True)
    j=json.loads(request("https://zenodo.org/api/records/8253921",2000000))
    v=j.get("files");es=v.get("entries",{}) if isinstance(v,dict) else {x["key"]:x for x in v}
    links=es["Data for Figures_Analyses .zip"]["links"]
    b=request(links.get("content") or links["self"])
    assert hashlib.sha256(b).hexdigest()==ZIPSHA
    with zipfile.ZipFile(io.BytesIO(b)) as zipf:
        match=[x for x in zipf.namelist() if Path(x).name=="FlowerLongevity_2015.csv"]
        assert len(match)==1
        c=zipf.read(match[0])
    assert hashlib.sha256(c).hexdigest()==CSVSHA
    d=pd.read_csv(io.BytesIO(c))
    assert len(d)==905
    assert {"RametID","Sex","Year","Date","FloweringNo","FlowerNo"}<=set(d)
    out={"status":"original_2015_source_structure_audit",
         "source_article_doi":"10.1002/ajb2.16239","author_data_doi":"10.5281/zenodo.8253921",
         "archive_SHA256":ZIPSHA,"table_SHA256":CSVSHA,
         "n_records":len(d),"n_ramet":int(d.RametID.nunique()),
         "n_sex":{str(k):int(v) for k,v in d.groupby("Sex").RametID.nunique().items()},
         "year":{str(k):int(v) for k,v in d.Year.value_counts().items()},
         "source_date_example":str(d.Date.iloc[0]),
         "n_source_dates":int(d.Date.nunique()),
         "n_duplicates_by_ramet_year_date":int(d.duplicated(["RametID","Year","Date"]).sum()),
         "n_complete_duplicate_rows":int(d.duplicated().sum()),
         "n_ramet_with_both_sex_codes":int((d.groupby("RametID").Sex.nunique()>1).sum()),
         "n_dates_per_ramet":{
             str(k):int(v) for k,v in d.groupby("RametID").Date.nunique().value_counts().items()},
         "original_date_level_counts":{str(k):int(v) for k,v in d.Date.value_counts().sort_index().items()},
         "fields":{}}
    for col in ["FloweringNo","FlowerNo"]:
        n=pd.to_numeric(d[col],errors="coerce")
        out["fields"][col]={"n_missing_or_non_numeric":int(n.isna().sum()),
            "n_zero":int((n==0).sum()),"min":float(n.min()),"max":float(n.max()),
            "mean":float(n.mean()),"n_negative":int((n<0).sum())}
    groups=d.groupby("RametID")
    out["flower_total_invariant_by_ramet"]={
        "n_with_different_totals":int((groups.FlowerNo.nunique()>1).sum()),
        "n_missing_totals":int(d.FlowerNo.isna().sum())}
    out["source_year_is_2015"]=(d.Year.nunique()==1 and int(d.Year.iloc[0])==2015)
    dt=[]
    for val in d.Date:
        raw=str(val).strip()
        try:
            # explicitly year from source column; M/D dates only
            date=datetime.datetime.strptime(f"2015/{raw}","%Y/%m/%d")
            dt.append(date.timetuple().tm_yday)
        except Exception:dt.append(None)
    out["date_parse_M_D_2015"]={"n_ok":int(sum(x is not None for x in dt)),
        "n_missing":int(sum(x is None for x in dt)),
        "min_doy":int(min(x for x in dt if x is not None)),
        "max_doy":int(max(x for x in dt if x is not None))}
    # If FlowerNo is not invariant, avoid automatic per-total normalization.
    out["can_normalize_per_ramet_total"]=(
        out["flower_total_invariant_by_ramet"]["n_with_different_totals"]==0 and
        out["fields"]["FlowerNo"]["n_missing_or_non_numeric"]==0 and
        out["fields"]["FlowerNo"]["min"]>0)
    (OUT/"daphne_2015_natural_flower_census_schema_v22.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_2015_SEX_SOURCE_V22",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
