#!/usr/bin/env python3
"""v0.56 source-archive audit of the identification priors and alternative outputs.

Reads exact pinned author ZIP (Zenodo 8253921) via v28 and summarizes only
metadata and aggregate paternity-confidence diagnostics. No original IDs,
genotypes, exact coordinates, or parentage rows are printed or saved.
"""
from __future__ import annotations
import hashlib, io, json, re, zipfile, collections
from pathlib import Path
import numpy as np
import pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data, get, API, ZIP_SHA, stringify

def sourcezip():
    meta=json.loads(get(API))
    fs=meta["files"]
    e=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    u=e["Data for Figures_Analyses .zip"]["links"]
    blob=get(u.get("content") or u["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    return zipfile.ZipFile(io.BytesIO(blob))

def source_index(z):
    entries=[]
    for obj in z.infolist():
        name=Path(obj.filename).name
        if not name or name.startswith("."):continue
        ext=Path(name).suffix.lower()
        if ext not in (".csv",".r",".rmd",".txt",".pdf",".xlsx",".xls",".png",".jpg",".dat",".rds",".rdata",".zip",".tsv"):continue
        # Redact any source path components; no source row values.
        category=(
          "paternity_csv" if re.search(r"^NF2018|^NF2020",name,re.I) else
          "adult_genotype" if "geno_phenotype" in name.lower() else
          "genotype_other" if "geno" in name.lower() or "allele" in name.lower() else
          "seed_codebook" if "sample_data" in name.lower() else
          "analysis_script" if ext in (".r",".rmd") else
          "other")
        entries.append(dict(filename=name,bytes=obj.file_size,category=category,
                            member_ext=ext,crc32=str(obj.CRC)))
    print("PINNED_ARCHIVE_CONTENT_INDEX_V56",json.dumps({
       "archive_sha256":ZIP_SHA,"files":entries,"n_files":len(entries)},ensure_ascii=False))

def source_csv_genotype_inventory(z):
    for member in z.namelist():
        file=Path(member).name
        if not file.lower().endswith(".csv"):continue
        try:
            df=pd.read_csv(io.BytesIO(z.read(member)),encoding="utf-8-sig")
        except Exception as exc:
            print("SOURCE_CSV_READABILITY_V56",json.dumps({"name":file,"reason":type(exc).__name__}))
            continue
        cols=list(map(str,df.columns))
        alleles=[k for k in cols if re.match(r"^DpJe.*\.(x|y)$",k)]
        possible=[k for k in cols if re.search(r"offspring|mother|father|dnaid|sample|genet|seed|probability|sex|info",k,re.I)]
        ssr_pairs=len({re.sub(r"\.(x|y)$","",k) for k in alleles})
        print("SOURCE_CSV_GENETIC_SCHEMA_V56",json.dumps({
            "file":file,"rows":len(df),"cols":len(cols),"titles":cols if len(cols)<=75 else cols[:75],
            "ssr_allele_columns":len(alleles),"ssr_loci_count":ssr_pairs,
            "nonmissing_allele_cells":int(df[alleles].notna().sum().sum()) if alleles else 0,
            "candidate_join_fieldnames":possible[:25]},ensure_ascii=False))

def original_codebook_audit(z):
    hits=[]
    pattern=re.compile(r"MasterBayes|masterbayes|posterior|paternity|geneal|lambda|distance|probab|beta|seed|father|genotyp",re.I)
    for obj in z.infolist():
        nm=Path(obj.filename).name
        if not nm.lower().endswith((".r",".rmd",".txt")) or obj.file_size>120000:continue
        txt=z.read(obj.filename).decode("utf-8-sig",errors="replace")
        statements=[s.strip() for s in txt.splitlines()
                    if pattern.search(s) and len(s.strip())<270 and s.strip() and not s.strip().startswith("#")]
        # Share only codebook *hints*, never full original line. Capture
        # counts for code structure, and which model families appear.
        hit={k:len(re.findall(k,txt,re.I)) for k in
             ["MasterBayes","sire","patern","distance","beta","p95","0\\.95",
              "genotyp","MCMC","unsampled","posterior","SeedID","Probability","Natural",
              "Day open","Night open"]}
        hits.append(dict(filename=nm,length=len(txt),code_tokens=hit,
                         count_model_relevant_code_lines=len(statements)))
    print("ORIGINAL_METHOD_CODE_INDEX_V56",json.dumps({"scripts":hits},ensure_ascii=False))

def confidence_and_genotypes(d):
    ad=d["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    cols=[str(c) for c in ad.columns]
    geno_cols=[c for c in cols if re.search(r"allel|ssr|locus|DpJe|primer|[0-9]{3,}",c,re.I)]
    print("ORIGINAL_ADULT_SCHEMA_V56",json.dumps({
      "adults":len(ad),"columns_total":len(cols),
      "column_titles":cols[:110],
      "possible_genotype_columns":geno_cols[:80]},ensure_ascii=False))
    for year in (2018,2020):
        raw=d[f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv"]
        assert len(raw)==(257 if year==2018 else 327)
        prob=pd.to_numeric(raw.Probability,errors="coerce")
        rs={"year":year,"source_rows":len(raw),
            "matched_parent_source_rows":int(raw.Father.notna().sum()),
            "nonmissing_probability_rows":int(prob.notna().sum()),
            "prob_min_max": [float(prob.min()),float(prob.max())],
            "prob_quantiles":{str(z):float(prob.quantile(z)) for z in [0,.01,.05,.1,.25,.5,.75,.9,.95,.99,1]},
            "all_assigned_by_threshold":{},
            "strict_H_father_assigned_by_threshold":{},
            "rows_with_no_assigned_father":int(raw.Father.isna().sum()),
            }
        s=raw.assign(father_key=stringify(raw.Father))
        ad0=ad[["GenetID2",f"Sex_{year}.g",f"NaturalFlowerNo_{year}.g"]].copy()
        ad0["father_key"]=stringify(ad0.GenetID2)
        ad0["flowers"]=pd.to_numeric(ad0[f"NaturalFlowerNo_{year}.g"],errors="coerce")
        ad0["eligible_h"]=ad0[f"Sex_{year}.g"].astype(str).eq("H") & ad0.flowers.gt(0)
        assert not ad0.father_key.duplicated().any()
        joined=s.merge(ad0[["father_key","eligible_h"]],on="father_key",how="left",validate="m:1")
        for threshold in [0.5,.8,.9,.95,.975,.99,.995]:
            selected=prob.ge(threshold)&raw.Father.notna()
            rs["all_assigned_by_threshold"][str(threshold)]=int(selected.sum())
            rs["strict_H_father_assigned_by_threshold"][str(threshold)]=int((selected&joined.eligible_h.eq(True)).sum())
        print("ORIGINAL_POSTERIOR_CONFIDENCE_V56",json.dumps(rs,ensure_ascii=False))

if __name__=="__main__":
    with sourcezip() as z:
        source_index(z)
        source_csv_genotype_inventory(z)
        original_codebook_audit(z)
    confidence_and_genotypes(data())
