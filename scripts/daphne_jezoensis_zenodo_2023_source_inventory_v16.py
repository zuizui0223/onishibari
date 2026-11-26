#!/usr/bin/env python3
"""Inventory author-deposited Daphne jezoensis 2023 day/night source data.

This is a source discovery audit; no inference before parsing data provenance.
Never print public individual IDs, geolocations, or full record rows.
"""
from __future__ import annotations
import csv, hashlib, io, json, re, urllib.request, urllib.error, zipfile
from pathlib import Path
import pandas as pd

RECORD_ID = 8253921
URL = f"https://zenodo.org/api/records/{RECORD_ID}"
OUT = Path("analysis_results_daphne_zenodo2023_v16")
MAX_FILE = 8_000_000

def get(url, cap=MAX_FILE):
    req=urllib.request.Request(url,headers={
        "User-Agent":"ScientificSourceProvenanceAudit/1.0 (published open-data replication)",
        "Accept":"application/json, application/octet-stream, */*"})
    with urllib.request.urlopen(req,timeout=60) as r:
        b=r.read(cap+1)
        if len(b)>cap: raise ValueError("oversize response")
        return b

def parsed(b, name):
    o={"file_name":str(name),"file_bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}
    kind=name.lower()
    if kind.endswith((".csv",".tsv",".txt")):
        for enc in ["utf-8-sig","cp932","latin-1"]:
            try:
                s=b.decode(enc)
                sep= "\t" if kind.endswith(".tsv") else None
                d=pd.read_csv(io.StringIO(s),sep=sep,engine="python")
                o["format"]="table";o["encoding"]=enc
                o["rows"]=len(d)
                o["columns"]=list(map(str,d.columns))
                o["column_non_missing"]={str(k):int(v) for k,v in d.notna().sum().items()}
                return o
            except Exception:continue
        o["parse"]="unsupported_table_encoding_or_delimiter"
    elif kind.endswith((".xls",".xlsx")):
        try:
            x=pd.ExcelFile(io.BytesIO(b))
            o["format"]="excel"
            o["sheets"]={}
            for sheet in x.sheet_names[:25]:
                d=pd.read_excel(x,sheet_name=sheet)
                o["sheets"][sheet]={"rows":len(d),"columns":list(map(str,d.columns)),
                "column_non_missing":{str(k):int(v) for k,v in d.notna().sum().items()}}
        except Exception as e:o["parse_error"]=type(e).__name__+":"+str(e)[:170]
    elif kind.endswith((".r",".rmd",".rscript",".py")):
        o["format"]="analysis_code"
        s=b.decode("utf-8-sig","replace")
        o["code_lines"]=len(s.splitlines())
        o["terms"]={k:len(re.findall(k,s,re.I)) for k in
         ("fruit","flower","pollen","paternity","night","day","glmm","lmer","female","year")}
        o["read_file_calls"]=list(dict.fromkeys(re.findall(r'(?:read\\.csv|read\\.table|read_xlsx|read_excel)\\s*\\(\\s*["\\x27]([^"\\x27]+)',s)))[:40]
    elif kind.endswith((".rds",".rda",".rdata")):
        o["format"]="R_serialized_not_yet_parsed"
    elif zipfile.is_zipfile(io.BytesIO(b)):
        o["format"]="zip"
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            o["members"]=[]
            for info in z.infolist()[:50]:
                if not info.is_dir() and info.file_size<=MAX_FILE:
                    o["members"].append(parsed(z.read(info),Path(info.filename).name))
    else:o["format"]="not_typed_as_table_or_code"
    return o

def main():
    OUT.mkdir(exist_ok=True)
    record=json.loads(get(URL, 3_000_000))
    files=record.get("files",{})
    if isinstance(files,dict): files=files.get("entries",{})
    if isinstance(files,dict): entries=[dict({"key":key},**val) for key,val in files.items()]
    elif isinstance(files,list):entries=files
    else:raise ValueError("Zenodo file entry missing")
    report={
        "source_record":RECORD_ID,"source_doi":"10.5281/zenodo.8253921",
        "source_article_doi":"10.1002/ajb2.16239",
        "published_source_title":record.get("metadata",{}).get("title"),
        "official_api":"zenodo.org/api/records/8253921",
        "n_files":len(entries),
        "files":[],
        "status":"source_inventory_completed_without_individual_rows"
    }
    if not entries:raise ValueError("No official public source files.")
    for entry in entries:
        name=entry.get("key") or entry.get("filename") or "unknown"
        links=entry.get("links") or {}
        link=links.get("content") or links.get("self") or entry.get("links",{}).get("download")
        if not link:link=f"https://zenodo.org/api/records/{RECORD_ID}/files/{name}/content"
        meta={"key":name, "declared_size":entry.get("size"),"declared_checksum":entry.get("checksum"),
              "url_origin":"zenodo.org"}
        size=entry.get("size",0)
        if size and size>MAX_FILE:
            meta["parse"]="SKIPPED_LIMIT"
        else:
            try:
                b=get(link)
                meta.update(parsed(b,name))
                ck=entry.get("checksum")
                if ck and ":" in ck:
                    alg,exp=ck.split(":",1)
                    if alg.lower() in ("md5","sha256"):
                        got=hashlib.new(alg.lower(),b).hexdigest()
                        if got.lower()!=exp.lower():raise ValueError("source checksum mismatch")
                        meta["declared_checksum_verified"]=True
            except Exception as e:
                meta["fetch_or_parse_error"]=type(e).__name__+":"+str(e)[:180]
        report["files"].append(meta)
    if not any(x.get("sha256") for x in report["files"]):
        raise RuntimeError("Could not retrieve any source file.")
    (OUT/"shibata_kudo_2023_zenodo_inventory_v16.json").write_text(
      json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_2023_SOURCE_INVENTORY_V16",json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
