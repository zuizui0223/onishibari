#!/usr/bin/env python3
"""Baseline data-coverage audit for reproductive trait decoupling research.

Sources:
  Wang et al. (2021) sexual system Dryad replica on Zenodo, record 4895147.
  Bennett et al. (2018) GloPL original Dryad, doi:10.5061/dryad.dt437.
  DiveRS PROTEUS 2025 authors' analysis GitHub (cleaned_data.csv).
This script does not treat genus-imputed states as observed species data,
or interpret observational joins as evidence of adaptation/speciation.
"""
from __future__ import annotations
from collections import Counter
import csv, hashlib, io, json, re, time, urllib.parse, urllib.request
from pathlib import Path

OUT=Path("analysis_results_broad")
OUT.mkdir(exist_ok=True)
UA={"User-Agent":"angiosperm-reproductive-decoupling/0.1 scientific research"}
SOURCES={
  "wang":"https://zenodo.org/api/records/4895147",
  "glopl":"https://doi.org/10.5061/dryad.dt437",
  "proteus":"https://raw.githubusercontent.com/divers-it/rs-traitspace/main/supplementary_data/cleaned_data.csv",
}
HEADERS={}
NAMES=["GloSexSyst.csv","SplistWithSexSyst.csv","TransitionBetweenSexSyst.csv"]

def get(url,limit=9_000_000,timeout=90):
    errors=[]
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers=UA)
            with urllib.request.urlopen(req,timeout=timeout) as r:
                b=r.read(limit+1)
                if len(b)>limit:raise ValueError(f"response too large {url}")
                return b,r.geturl()
        except Exception as e:
            errors.append(str(e))
            if retry<2:time.sleep(2**retry)
    raise RuntimeError(f"URL inaccessible: {url}: {errors[-1]}")

def csvrows(raw):
    text=raw.decode("utf-8-sig",errors="replace")
    dialect=csv.Sniffer().sniff(text[:4096],delimiters=",\t;")
    r=csv.DictReader(io.StringIO(text), dialect=dialect)
    rows=list(r)
    return r.fieldnames, rows

def clean_name(name):
    if name is None:return ""
    s=str(name).strip().replace("_"," ")
    s=re.sub(r"\s+"," ",s)
    s=re.sub(r"[\u00d7×]","x",s)
    return s.lower().strip()

def audit_table(source,filename,raw,license_note=""):
    cols,rows=csvrows(raw)
    if not cols or not rows:raise ValueError(f"empty data {filename}")
    return {
        "source":source,"file":filename,"bytes":len(raw),
        "sha256":hashlib.sha256(raw).hexdigest(),"rows":len(rows),
        "columns":cols,"examples":rows[:2],"license_note":license_note,
    },rows

def find_glopl():
    """Prefer original Dryad GloPL.csv; author GitHub code repository is fallback."""
    url="https://datadryad.org/api/v2/datasets/"+urllib.parse.quote("doi:10.5061/dryad.dt437",safe="")
    errors=[]
    try:
        blob,_=get(url,limit=2_000_000)
        meta=json.loads(blob)
        (OUT/"dryad_glopl_metadata.json").write_text(json.dumps(meta,indent=2)[:200000])
        version_endpoint=(meta.get("_links",{}).get("stash:version",{}) or {}).get("href")
        if not version_endpoint:raise ValueError("Dataset lacks stash:version link")
        if not version_endpoint.startswith("http"):
            version_endpoint="https://datadryad.org"+version_endpoint
        raw,_=get(version_endpoint,limit=2_000_000)
        vmeta=json.loads(raw)
        (OUT/"dryad_glopl_version.json").write_text(json.dumps(vmeta,indent=2)[:200000])
        vlinks=vmeta.get("_links",{})
        endpoints=[
            item.get("href") for key,item in vlinks.items()
            if "files" in key.lower() and isinstance(item,dict) and item.get("href")]
        if "id" in vmeta:endpoints.append(f"https://datadryad.org/api/v2/versions/{vmeta['id']}/files")
        for p in endpoints:
            try:
                if not p.startswith("http"):p="https://datadryad.org"+p
                bs,_=get(p,limit=2_000_000)
                j=json.loads(bs)
                listing=j.get("_embedded",{}).get("stash:files",[]) or j.get("files",[])
                for x in listing:
                    if x.get("path")!="GloPL.csv" and x.get("filename")!="GloPL.csv":continue
                    fid=x.get("id")
                    href=(x.get("_links",{}).get("stash:download",{}) or {}).get("href")
                    for dl in [href,f"https://datadryad.org/api/v2/files/{fid}/download"]:
                        if not dl:continue
                        if not dl.startswith("http"):dl="https://datadryad.org"+dl
                        try:
                            data,actual=get(dl,limit=4_000_000)
                            if data[:100].lstrip().startswith(b"<html"):continue
                            return data,{"api":url,"download":actual,"release":"Dryad original"}
                        except Exception as e:errors.append("file download: "+str(e))
            except Exception as e:errors.append("file listings: "+str(e))
    except Exception as e:errors.append("Dryad API: "+str(e))
    # Curated working copy from the original GloPL authors' analysis repository.
    fallback=("https://raw.githubusercontent.com/idiv-biodiversity/"
      "pollen-limitation-data-descriptor/master/Data/GloPL_with_id.csv")
    data,actual=get(fallback,limit=4_000_000)
    if not data[:100].lstrip() or data[:100].lstrip().startswith(b"<html"):
        raise ValueError("GloPL author repo did not return valid CSV")
    return data,{"api":url,"download":actual,"release":"Authors' curated GitHub working copy",
                 "original_dryad_errors":errors}

def profile():
    report={"as_of":"2026-11-15","sources":SOURCES,"tables":[],
            "join_analysis":"exact normalized binomial names only; no synonym resolution; no imputed cross-mapping",
            "join_counts":{},"errors":[], "done":False}
    tables={}
    try:
        meta_blob,zenodo_url=get(SOURCES["wang"],limit=3_000_000)
        meta=json.loads(meta_blob)
        report["wang_release"]={"record_id":meta.get("id"),
                                "files":[{"key":f.get("key"),"size":f.get("size")} for f in meta.get("files",[])]}
        wang_files={x.get("key"):x for x in meta.get("files",[])}
        for filename in NAMES:
            if filename not in wang_files:raise ValueError("Missing Wang dataset: "+filename)
            f=wang_files[filename]
            url=(f.get("links") or {}).get("self")
            if not url:raise ValueError("Missing Wang file URL "+filename)
            raw,actual=get(url,limit=6_000_000)
            if len(raw)!=f.get("size"):raise ValueError("Wang file size mismatch "+filename)
            details,rows=audit_table("Wang2021",filename,raw,
                  "2021 sexual system labels have distinct genus/species evidence; do not assume each row independently measured")
            details["download"]=actual
            report["tables"].append(details)
            tables[filename]=rows
        raw,actual=get(SOURCES["proteus"],limit=1_000_000)
        d,rows=audit_table("PROTEUS_DiveRS","cleaned_data.csv",raw,
             "Authors' cleaned research dataset, may exclude species or recode traits")
        d["download"]=actual
        report["tables"].append(d)
        tables["PROTEUS"]=rows
        try:
            raw,gloplmeta=find_glopl()
            d,rows=audit_table("GloPL_2018",("GloPL.csv" if gloplmeta["release"]=="Dryad original" else "GloPL_with_id.csv"),raw,
                    "Pollen supplementation experiments; each row can be site-year-treatment or morph")
            d["download"]=gloplmeta["download"]
            d["release"]=gloplmeta["release"]
            report["tables"].append(d)
            tables["GloPL"]=rows
            report["glopl_access"]="successful"
        except Exception as ex:
            report["glopl_access"]="unavailable"
            report["errors"].append("GloPL download: "+str(ex))
        spec_names=set()
        for row in tables["SplistWithSexSyst.csv"]:
            vals=row.values()
            species_col=next((k for k in row.keys() if "species" in (k or "").lower() or "binomial" in (k or "").lower()),None)
            if species_col: name=row.get(species_col)
            else: name=next((v for v in vals if v and re.match(r"^[A-Z][a-z]+[ _][a-z-]+",v)),None)
            if name: spec_names.add(clean_name(name))
        prot_rows=tables["PROTEUS"]
        p_names=[clean_name(row.get("") or row.get("Species") or row.get("species")) for row in prot_rows]
        w_names=spec_names
        pset=set(filter(None,p_names))
        report["join_counts"]["Wang_full_species_unique"]=len(w_names)
        report["join_counts"]["PROTEUS_unique"]=len(pset)
        report["join_counts"]["Wang_full_to_PROTEUS_exact"]=len(pset&w_names)
        report["join_counts"]["PROTEUS_sex_system_raw_categories"]=dict(Counter(
          str(r.get("SexualSystem","NA")) for r in prot_rows))
        report["join_counts"]["PROTEUS_mating_raw_categories"]=dict(Counter(
          str(r.get("Mating","NA")) for r in prot_rows))
        report["join_counts"]["PROTEUS_Pollination_raw_categories"]=dict(Counter(
          str(r.get("Pollination","NA")) for r in prot_rows))
        report["join_counts"]["PROTEUS_flower_sex_categories"]=dict(Counter(
          str(r.get("FlowerSex","NA")) for r in prot_rows))
        report["join_counts"]["PROTEUS_Lifespan_categories"]=dict(Counter(
          str(r.get("Lifespan","NA")) for r in prot_rows))
        report["join_counts"]["overlap_examples"]=sorted(pset&w_names)[:20]
        # GloPL species column is audited but NOT assumed when unrecognized.
        if "GloPL" in tables:
            gro=tables["GloPL"]
            keys=list(gro[0]) if gro else []
            report["join_counts"]["GloPL_species_column_candidates"]=[
                k for k in keys if any(s in k.lower() for s in ["species","taxon","plantname","binomial"])]
        report["done"]=True
    except Exception as ex:
        report["errors"].append("Unrecoverable source/profile error: "+type(ex).__name__+": "+str(ex))
    (OUT/"coverage_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,default=str)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False,default=str)[:60000])
    if not report["done"]:raise RuntimeError(report["errors"][-1])
if __name__=="__main__":profile()
