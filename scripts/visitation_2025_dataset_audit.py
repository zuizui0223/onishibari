#!/usr/bin/env python3
"""Inspect an effort-standardized flower visitation/seed-set public dataset.

Chacko, Moretti & Frey 2025, EnviDat 10.16904/envidat.676.
Only headers, counts and non-spatial metadata are exported. Do not republish
precise coordinates of private home gardens as new public GitHub artifacts.
"""
import csv,hashlib,io,json,time,urllib.request,zipfile
from collections import Counter
from pathlib import Path

OUT=Path("analysis_results_pollinator_scarcity")
OUT.mkdir(exist_ok=True)
DATA="https://www.envidat.ch/dataset/d8e61f31-c2e5-4052-8784-b1624dfbfc40/resource/cef4a88c-0ef6-4327-a4a0-fc7d33b468e5/download/rejichacko_etal_2025_envidat.zip"

def load():
    for i in range(3):
        try:
            req=urllib.request.Request(DATA,headers={"User-Agent":"onishibari-pollinator-scarcity-audit/0.1"})
            with urllib.request.urlopen(req,timeout=120) as r:
                raw=r.read(80_000_001)
            if len(raw)>80_000_000:raise ValueError("data ZIP exceeds max expected")
            if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("Not ZIP response")
            return raw
        except Exception:
            if i==2:raise
            time.sleep(2**i)

def meta_for_csv(raw,path):
    text=raw.decode("utf-8-sig",errors="replace")
    reader=csv.DictReader(io.StringIO(text))
    cols=reader.fieldnames or []
    rows=list(reader)
    if len(cols)>300:raise ValueError("Unexpectedly wide table")
    nonmissing=Counter()
    for row in rows:
        for k in cols:
            if str(row.get(k,"")).strip() not in ("","NA","nan","NaN"):
                nonmissing[k]+=1
    return dict(file=path,rows=len(rows),columns=cols,
        nonmissing_by_col=dict(nonmissing),
        column_count=len(cols),
        source_sha256=hashlib.sha256(raw).hexdigest())

def main():
    result={"study":"Reji Chacko et al. 2025 EnviDat DOI 10.16904/envidat.676",
       "source_url":DATA,"site_data_private_gardens":True,
       "export_policy":"No private garden coordinates or raw observer identifier rows exported",
       "status":"failure","files":[],"csv_profiles":[]}
    try:
        raw=load()
        result["download_bytes"]=len(raw)
        result["zip_sha256"]=hashlib.sha256(raw).hexdigest()
        with zipfile.ZipFile(io.BytesIO(raw)) as zipf:
            infos=[i for i in zipf.infolist() if not i.is_dir()]
            result["files"]=[dict(path=i.filename,bytes=i.file_size) for i in infos]
            for info in infos:
                if info.filename.lower().endswith(".csv") and info.file_size<=20_000_000:
                    if "site_coordinates" in info.filename:
                        # Read source schema only, never output location values.
                        text=zipf.read(info).decode("utf-8-sig",errors="replace")
                        reader=csv.DictReader(io.StringIO(text))
                        colnames=reader.fieldnames or []
                        n_rows=sum(1 for _ in reader)
                        result["csv_profiles"].append(dict(file=info.filename,rows=n_rows,
                             columns=colnames,location_values_redacted=True,
                             description="Only header and row count; private garden coordinates not exported"))
                        continue
                    result["csv_profiles"].append(meta_for_csv(zipf.read(info),info.filename))
        result["status"]="success"
    except Exception as e:
        result["error"]=type(e).__name__+": "+str(e)
    (OUT/"public_visitation_dataset_inventory.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"error":result.get("error"),"files":result["files"],"csv_profiles":result["csv_profiles"]},indent=2)[:30000])
    if result["status"]!="success":raise RuntimeError(result["error"])
if __name__=="__main__":main()
