#!/usr/bin/env python3
"""v0.17 source diagnostic for 2023 Daphne jezoensis diurnal/nocturnal treatment data.

Author data: Zenodo 8253921; published Shibata & Kudo (2023), AJB 16239.
Produces aggregated coding and unique plant/cohort diagnostics only.
NO biological effects, raw IDs, genotype rows, coordinates, or individual results.
"""
from __future__ import annotations
import io, json, hashlib, urllib.request, zipfile, re
from pathlib import Path
from collections import Counter
import pandas as pd
import numpy as np

API="https://zenodo.org/api/records/8253921"
ZIP_EXPECT="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
OUT=Path("analysis_results_daphne_diel_design_v17")
TARGETS=[
 "Nopporo_daynoon.csv","Nopporo_control.csv","Nopporo_outcross.csv",
 "FlowerLongevity.csv","TotalVolatiles.csv",
 "GenetID_geno_phenotype_1820Nopadult_211109.csv",
 "NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
 "NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"
]
def get(url, cap=8000000):
    req=urllib.request.Request(url,headers={"User-Agent":"Daphne-2023-published-source-aggregate-v17","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as h:b=h.read(cap+1)
    if len(b)>cap:raise RuntimeError("source cap exceeded")
    return b
def read():
    j=json.loads(get(API,2500000))
    fs=j.get("files",{})
    f=fs.get("entries",{}) if isinstance(fs,dict) else {a.get("key"):a for a in fs}
    def url(key):
        return f[key]["links"].get("content") or f[key]["links"].get("self")
    zbyte=get(url("Data for Figures_Analyses .zip"))
    assert hashlib.sha256(zbyte).hexdigest()==ZIP_EXPECT, "Source zip checksum mismatch"
    files={}
    with zipfile.ZipFile(io.BytesIO(zbyte)) as z:
        for n in z.namelist():
            short=Path(n).name
            if short in TARGETS:
                b=z.read(n)
                files[short]=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig")
    rcode={}
    for key in ["Figures_Analyses.R","Paternity analyses.R"]:
        if key in f:
            b=get(url(key))
            rcode[key]=b.decode("utf8",errors="replace")
    return files,rcode

def normal_value(v):
    if pd.isna(v):return "NA"
    return str(v)[:45]
def summed(x,col):
    if col not in x:return None
    a=pd.to_numeric(x[col],errors="coerce")
    return {"n_numeric":int(a.notna().sum()),"n_negative":int((a<0).sum()),
            "n_fractional":int((a.notna()&((a%1)!=0)).sum()),
            "min":float(a.min()) if a.notna().any() else None,
            "max":float(a.max()) if a.notna().any() else None,
            "sum":round(float(a.sum()),6)}
def counts(x,col):
    if col not in x:return None
    return {normal_value(k):int(v) for k,v in x[col].value_counts(dropna=False).items()}
def profile(d,name):
    v={"name":name,"n_rows":len(d),"columns":list(d.columns),
       "key_levels":{c:counts(d,c) for c in ["Year","Sex","Treatment","Use"] if c in d},
       "n_unique_ramet":int(d.RametID.nunique()) if "RametID" in d else None,
       "n_unique_inflo":int(d.InfloID.nunique()) if "InfloID" in d else None,
       "n_unique_year_ramet":int(d[["Year","RametID"]].drop_duplicates().shape[0]) if {"Year","RametID"}<=set(d) else None,
       "n_unique_year_ramet_treatment":int(d[["Year","RametID","Treatment"]].drop_duplicates().shape[0]) if {"Year","RametID","Treatment"}<=set(d) else None,
       "flower":summed(d,"FlowerNo"),"fruit":summed(d,"FruitNo"),
       "duplicate_complete_rows":int(d.duplicated().sum())}
    if "FlowerNo" in d and "FruitNo" in d:
        f=pd.to_numeric(d.FlowerNo,errors="coerce")
        p=pd.to_numeric(d.FruitNo,errors="coerce")
        v["n_fruit_gt_flower"]=int(((p>f)&p.notna()&f.notna()).sum())
        v["n_flower_zero"]=int((f==0).sum())
        v["n_flower_missing"]=int(f.isna().sum())
        v["n_fruit_missing"]=int(p.isna().sum())
        valid=d[f.notna()&p.notna()&(f>0)&(p>=0)&(p<=f)].copy()
        v["analysis_eligible_rows"]=len(valid)
        v["groups"]=[]
        keys=[c for c in ["Year","Sex","Treatment"] if c in valid]
        if keys:
            for k,g in valid.groupby(keys,dropna=False):
                if not isinstance(k,tuple):k=(k,)
                v["groups"].append({"levels":dict(zip(keys,map(normal_value,k))),
                  "rows":len(g),
                  "unique_ramets":int(g.RametID.nunique()) if "RametID" in g else None,
                  "flowers":int(pd.to_numeric(g.FlowerNo).sum()),
                  "fruits":int(pd.to_numeric(g.FruitNo).sum()),
                  "rows_with_fruits":int((pd.to_numeric(g.FruitNo)>0).sum())})
        if {"Year","RametID","Treatment"}<=set(valid.columns):
            by=valid.groupby(["Year","RametID"])["Treatment"].nunique()
            v["n_year_ramet_multiple_treatments"]=int((by>1).sum())
            v["n_year_ramet_only_one_treatment"]=int((by==1).sum())
            q=valid.groupby(["Year","RametID","Treatment"]).size()
            v["n_year_ramet_treatment_repeats"]=int((q>1).sum())
            v["max_rows_per_year_ramet_treatment"]=int(q.max())
        elif {"RametID","Treatment"}<=set(valid.columns):
            by=valid.groupby(["RametID"])["Treatment"].nunique()
            v["n_ramet_multiple_treatments"]=int((by>1).sum())
            q=valid.groupby(["RametID","Treatment"]).size()
            v["n_ramet_treatment_repeats"]=int((q>1).sum())
    return v

def comparison(files):
    v={}
    for key in ["Nopporo_daynoon.csv","Nopporo_control.csv","Nopporo_outcross.csv"]:
        d=files.get(key)
        if d is None:continue
        for other in ["Nopporo_daynoon.csv","Nopporo_control.csv","Nopporo_outcross.csv"]:
            if other<=key:continue
            b=files.get(other)
            if b is None:continue
            common={"RametID","Year"}<=set(d)&set(b)
            if common:
                a=set(zip(d["Year"].astype(str),d["RametID"].astype(str)))
                z=set(zip(b["Year"].astype(str),b["RametID"].astype(str)))
                v[f"{key} vs {other}"]={"shared_year_ramet_n":len(a&z),
                    "only_first":len(a-z),"only_second":len(z-a)}
    return v

def code_excerpt(rcode):
    excerpts={}
    terms=r"(?:Nopporo_daynoon|Nopporo_control|Nopporo_outcross|night|day|noon|fruit|paternity|treatment)"
    for name,content in rcode.items():
        lines=content.splitlines()
        filtered=[]
        for idx,line in enumerate(lines):
            if re.search(terms,line,re.I) and len(line)<180:
                # No source-level data records here, only analysis code / comments
                if line.lstrip().startswith(("#","read.","data.","write.","glm(","glmer(","lmer(","ggplot(","levels(","factor(")) or "<-" in line:
                    filtered.append({"line":idx+1,"text":line.strip()[:160]})
        excerpts[name]={"n_lines":len(lines),"candidate_lines":filtered[:70]}
    return excerpts

def main():
    OUT.mkdir(exist_ok=True)
    files,rcode=read()
    needed={"Nopporo_daynoon.csv","Nopporo_control.csv","Nopporo_outcross.csv"}
    assert needed<=set(files),f"missing source {needed-set(files)}"
    data={"status":"source_verified_design_aggregate_no_new_effect",
          "source":"https://zenodo.org/records/8253921",
          "sha256_original_zip":ZIP_EXPECT,
          "profiles":[profile(d,n) for n,d in files.items()],
          "file_overlap_year_ramet":comparison(files),
          "original_R_script_context":code_excerpt(rcode)}
    (OUT/"design_diagnostic_v17.json").write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n")
    printable={k:v for k,v in data.items() if k!="original_R_script_context"}
    for p in printable["profiles"]:
        if p["name"] not in needed:
            p.pop("groups",None)
    print("DAPHNE_DIEL_DESIGN_V17",json.dumps(printable,ensure_ascii=False),flush=True)
    print("ORIGINAL_R_CODE_CONTEXT_V17",json.dumps(data["original_R_script_context"],ensure_ascii=False)[:24000],flush=True)

if __name__=="__main__":main()
