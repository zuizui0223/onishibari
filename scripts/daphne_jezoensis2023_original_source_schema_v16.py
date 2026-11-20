#!/usr/bin/env python3
"""Public original source profile for Zenodo 8253921. No effect estimates.

Read original files with bounded downloads; expose only aggregate schema,
never individual IDs, genotypes, locations, pollen counts, or sensitive data.
All source IDs, treatment semantics and experimental units require review.
"""
from __future__ import annotations
import hashlib
import io
import json
import pathlib
import re
import urllib.request
import zipfile
import pandas as pd

RECORD_ID = "8253921"
API = f"https://zenodo.org/api/records/{RECORD_ID}"
OUT = pathlib.Path("source_audit_daphne_jezoensis_diel_v16")
MAX_FILE = 12_000_000
MAX_TOTAL = 35_000_000
EXPECTED = {
    "Nopporo_daynoon.csv": "28f7c9d54cdf92e95815ba70c0308f6a6e9fd116ddf6396efa71f88ceb1cc7de",
    "Nopporo_control.csv": "8b36a936ac7880d6fd7abb9928cd1b0677ffc75e80f1abe6a967c0e67118b07b",
}

def get(url: str, cap=MAX_FILE) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent":"Reproductive-Stage-Source-Audit/0.16","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=90) as handle:
        b=handle.read(cap+1)
    if len(b)>cap: raise RuntimeError("source byte cap exceeded")
    return b

def nested_entries(metadata):
    raw = metadata.get("files", [])
    if isinstance(raw, list):
        return [(x.get("key") or x.get("filename") or "", x.get("links",{}).get("self")
                 or x.get("links",{}).get("content") or "") for x in raw]
    entries=raw.get("entries", {}) if isinstance(raw,dict) else {}
    return [(k, (v.get("links") or {}).get("content") or "")
            for k,v in entries.items()]

def profile(key,b):
    name=pathlib.PurePosixPath(key).name
    result={"name":name,"size":len(b),"sha256":hashlib.sha256(b).hexdigest()}
    if name in EXPECTED:
        result["expected_sha256_matches"]=(result["sha256"]==EXPECTED[name])
        if not result["expected_sha256_matches"]:
            raise RuntimeError("known original source hash changed: "+name)
    if name.lower().endswith(".zip"):
        result["archive_entries"]=[]
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            for x in z.infolist():
                if x.is_dir(): continue
                if x.file_size>MAX_FILE:
                    result["archive_entries"].append({"name":pathlib.PurePosixPath(x.filename).name,
                                                      "skipped":"size cap"})
                else:
                    result["archive_entries"].append(profile(x.filename,z.read(x)))
        return result
    if not name.lower().endswith((".csv",".tsv")): return result
    for encoding in ("utf-8-sig","cp932","latin-1"):
        try:
            separator="\t" if name.lower().endswith(".tsv") else ","
            d=pd.read_csv(io.BytesIO(b),encoding=encoding,sep=separator)
            break
        except (UnicodeDecodeError,UnicodeError):
            continue
    else:
        result["parse_status"]="unsupported_encoding"
        return result
    result["rows"]=int(len(d))
    result["columns"]=[str(c) for c in d.columns]
    result["missing_counts"]={str(k):int(v) for k,v in d.isna().sum().items()}
    result["duplicate_full_rows"]=int(d.duplicated().sum())
    result["field_profiles"]={}
    for col in d.columns:
        c=str(col)
        vals=d[col]
        info={"n_nonmissing":int(vals.notna().sum()),"n_distinct":int(vals.nunique(dropna=True))}
        # Only report level counts for predefined semantics without divulging genet IDs.
        if re.search(r"^(treatment|group|sex|year|experimental|exposure|condition|season)$",c,re.I):
            vc=vals.astype("string").value_counts(dropna=False)
            if len(vc)<=15:
                info["level_counts"]={str(k)[:50]:int(v) for k,v in vc.items()}
        result["field_profiles"][c]=info
    keys=[c for c in d.columns if re.search(r"(plant|ramet|genet|mother|individual|id|year|treat|day|night|sex)",str(c),re.I)]
    result["candidate_repeated_measure_fields"]=[str(c) for c in keys]
    if keys:
        # Row-level distinct counts, never release key values.
        for c in keys[:8]:
            result["field_profiles"][str(c)]["largest_cluster_n"]=int(d[c].value_counts(dropna=True).max()) if d[c].notna().any() else 0
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    metadata=json.loads(get(API,3_000_000))
    files=nested_entries(metadata)
    report={"status":"source_inventory_only_no_reproductive_fitness_estimates",
            "record":RECORD_ID,
            "doi":metadata.get("doi"),
            "title":(metadata.get("metadata") or {}).get("title") or metadata.get("title"),
            "published":(metadata.get("metadata") or {}).get("publication_date") or metadata.get("created"),
            "n_file_entries":len(files),"files":[]}
    total=0
    for name,url in files:
        if not url: continue
        if not (name.lower().endswith((".csv",".tsv",".zip"))): continue
        raw=get(url)
        total+=len(raw)
        if total>MAX_TOTAL:raise RuntimeError("total byte cap exceeded")
        report["files"].append(profile(name,raw))
    if not report["files"]:raise RuntimeError("no retrievable original structured files")
    (OUT/"source_schema_no_row_data_v16.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_DIEL_ORIGINAL_SCHEMA_V16", json.dumps(report,ensure_ascii=False)[:25000],flush=True)
    print("NO_FIT_PERFORMED: dates, units, denominators and treatment semantics require manual identification",flush=True)

if __name__=="__main__":
    main()
