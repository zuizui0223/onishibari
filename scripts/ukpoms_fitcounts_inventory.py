#!/usr/bin/env python3
"""Inventory official UK PoMS 2017-2023 FIT Count download.

Provenance:
  UK Pollinator Monitoring Scheme (2026),
  https://doi.org/10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618
  EIDC / OGL, acknowledges UKCEH, JNCC and governments.
No raw site coordinates or surveyor identifiers are republished.
This does NOT create an NDVI association or infer effective pollen delivery.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path

OUT=Path("analysis_results_ukpoms")
OUT.mkdir(parents=True,exist_ok=True)
IDENT="d706088c-9d48-4be5-b156-6ed6ad5e0618"
CANDIDATES=[
  f"https://data-package.ceh.ac.uk/data/{IDENT}",
  f"https://data-package.ceh.ac.uk/data/{IDENT}.zip",
  f"https://data-package.ceh.ac.uk/data/{IDENT}?download=true",
]
MAX_BYTES=100_000_000
HEAD={"User-Agent":"onishibari-reproductive-systems-sampling/0.1 (academic)"}

def attempt_source(url):
    req=urllib.request.Request(url,headers=HEAD)
    with urllib.request.urlopen(req,timeout=100) as res:
        content_type=res.headers.get("Content-Type","")
        length_header=res.headers.get("Content-Length")
        final_url=res.geturl()
        data=res.read(MAX_BYTES+1)
    if len(data)>MAX_BYTES:raise ValueError("Dataset download exceeds 100MB bound")
    return data,content_type,length_header,final_url

def safe_profile_zip(data):
    result=[]
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for item in z.infolist():
            if item.is_dir():continue
            info=dict(file=item.filename,bytes=item.file_size)
            if item.file_size>75_000_000:
                info["note"]="Not profiled because file exceeds conservative limit"
                result.append(info);continue
            name=item.filename.lower()
            if name.endswith((".csv",".tsv",".txt")):
                raw=z.read(item)
                info["sha256"]=hashlib.sha256(raw).hexdigest()
                if b"\x00" not in raw[:500] and name.endswith((".csv",".tsv")):
                    sample=raw[:2048].decode("utf-8-sig","replace")
                    delimiter="\t" if name.endswith(".tsv") else ","
                    rdr=csv.DictReader(io.StringIO(sample),delimiter=delimiter)
                    info["columns"]=rdr.fieldnames or []
                    # Count lines only where row structure likely to be CSV and no quoting newlines.
                    reader=csv.reader(io.TextIOWrapper(io.BytesIO(raw),encoding="utf-8-sig",errors="replace"),delimiter=delimiter)
                    try:
                        info["rows_excluding_header"]=max(0,sum(1 for _ in reader)-1)
                    except Exception as e:
                        info["row_count_error"]=str(e)[:140]
                    # Expose schema only; not sample-site information.
                    info["column_roles"]={
                       "date":[x for x in info["columns"] if any(y in x.lower() for y in ("date","year","month"))],
                       "pollinator":[x for x in info["columns"] if any(y in x.lower() for y in ("bee","fly","hover","insect","pollinator","butterfly","wasp"))],
                       "habitat":[x for x in info["columns"] if any(y in x.lower() for y in ("flower","habitat","land","vegetation","urban"))],
                       "effort":[x for x in info["columns"] if any(y in x.lower() for y in ("duration","time","effort","count"))],
                       "location":[x for x in info["columns"] if any(y in x.lower() for y in ("longitude","latitude","easting","northing","grid","site","square"))],
                    }
            result.append(info)
    return result

def main():
    report={
      "source_doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "study_year_range":[2017,2023],
      "source":"UK Pollinator Monitoring Scheme, 2026 EIDC",
      "license":"Open Government Licence, attribution/PoMS acknowledgement required",
      "visit_protocol":"10-minute observations within 0.5m x 0.5m quadrat on target flowers",
      "outcomes_not_measured":["pollen grain deposition per insect","seed set","plant reproductive biomass"],
      "coordinate_values_exported":False,"attempts":[],
      "status":"failed"}
    for uri in CANDIDATES:
        try:
            data,typ,length,final=attempt_source(uri)
            summary=dict(request=uri,final_url=final,
                         content_type=typ,advertised_length=length,downloaded_bytes=len(data),
                         sha256=hashlib.sha256(data).hexdigest(),
                         first_bytes_hex=data[:20].hex(),is_zip=zipfile.is_zipfile(io.BytesIO(data)))
            report["attempts"].append(summary)
            if summary["is_zip"]:
                report["archive_files"]=safe_profile_zip(data)
                report["data_url_verified"]=final
                report["status"]="success"
                break
            else:
                # EIDC may return JSON listing/directory with links rather than ZIP
                text=data[:20000].decode("utf-8","replace")
                report["non_zip_preview"]=text[:4000]
                try:
                    obj=json.loads(text)
                    report["response_keys"]=list(obj) if isinstance(obj,dict) else []
                except Exception:pass
        except Exception as exc:
            report["attempts"].append({"request":uri,"error":type(exc).__name__+": "+str(exc)[:600]})
        time.sleep(0.5)
    (OUT/"source_inventory.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:32000])
    if report["status"]!="success":
        raise RuntimeError("Official EIDC data file not yet retrievable through public endpoints; do not claim source rows")
if __name__=="__main__":main()
