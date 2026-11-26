#!/usr/bin/env python3
"""Original public Dryad Clarkia concinna five-study source quality audit (v0.13).

All source rows are consumed in-memory from official published-version ZIP.
Does not export original plant IDs, locations, or raw data. The source README
is read ONLY as research documentation, not as executable instructions.
"""
from __future__ import annotations
from collections import Counter
import hashlib,io,json,re,urllib.request,urllib.parse,zipfile
from pathlib import Path
import pandas as pd

OUT=Path("analysis_results_clarkia_concinna_design_v13");OUT.mkdir(exist_ok=True)
DOI="10.5061/dryad.8pk0p2nz2"
BASE="https://datadryad.org"
def get(url,cap=2_000_000):
    r=urllib.request.Request(url,headers={"User-Agent":"Reproduction-stage-clarkia-source-audit/0.13","Accept":"application/json, application/zip"})
    with urllib.request.urlopen(r,timeout=70) as res:
        blob=res.read(cap+1)
        if len(blob)>cap:raise ValueError("download cap violated")
        return blob
def main():
    url=BASE+"/api/v2/datasets/"+urllib.parse.quote("doi:"+DOI,safe="")
    meta=json.loads(get(url))
    ver=json.loads(get(url+"/versions"))
    versions=ver.get("_embedded",{}).get("stash:versions",[])
    ids=[]
    for v in versions:
        m=re.search(r"/versions/(\d+)",v.get("_links",{}).get("self",{}).get("href",""))
        if m:ids.append((int(v.get("versionNumber",0)),int(m.group(1))))
    if not ids:raise ValueError("official public Dryad versions absent")
    vid=max(ids)[1]
    raw=get(BASE+f"/api/v2/versions/{vid}/download",2_000_000)
    if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("unexpected published archive")
    out={"status":"failed","doi":DOI,"version_id":vid,"original_data_visibility":meta.get("visibility"),
         "license":meta.get("license"),"n_files":0,"files":[],"no_original_record_ids_exported":True}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for zi in z.infolist():
            if zi.is_dir():continue
            name=Path(zi.filename).name
            if not name:continue
            b=z.read(zi)
            file={"name":name,"source_bytes":len(b),"original_uncompressed_SHA256":hashlib.sha256(b).hexdigest()}
            if name.lower().endswith(".csv"):
                df=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig",on_bad_lines="warn")
                df.columns=[str(c).strip() for c in df.columns]
                file["n_rows"]=len(df);file["columns"]=list(df)
                file["missing_by_column"]={str(k):int(v) for k,v in df.isna().sum().items() if v}
                file["categorical_or_bounded_columns"]={}
                file["numeric_column_profiles"]={}
                for c in df:
                    x=df[c];counter=x.astype("string").str.strip().value_counts(dropna=False)
                    if 1<len(counter)<=14 and "notes" not in c.lower():
                        file["categorical_or_bounded_columns"][c]={str(a):int(b) for a,b in counter.items()}
                    nu=pd.to_numeric(x,errors="coerce")
                    if nu.notna().sum()>max(5,.75*len(df)):
                        file["numeric_column_profiles"][c]={"n_numeric":int(nu.notna().sum()),"mean":float(nu.mean()),
                             "min":float(nu.min()),"max":float(nu.max()),"n_zero":int((nu==0).sum())}
            elif name.lower().endswith(".md"):
                txt=b.decode("utf-8-sig","replace")
                file["readme_first_1800_chars"]=txt[:1800]
                file["readme_headings"]=re.findall(r"(?m)^#{1,4} (.+)$",txt)[:65]
            out["files"].append(file)
            print("CLARKIA_SOURCE_FILE",json.dumps(file,ensure_ascii=False)[:15000],flush=True)
    expected={f"Floral_display_and_geitonogamy_data_study_{i}.csv" for i in range(1,6)}
    got={f["name"] for f in out["files"]}
    if not expected<=got:raise ValueError("required original five study files unavailable "+str(expected-got))
    out["n_files"]=len(out["files"])
    out["status"]="success_original_source_schema_and_sha_checked"
    (OUT/"clarkia_concinna_five_study_source_audit_v13.json").write_text(
       json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("CLARKIA_SOURCE_AUDIT_COMPLETE",out["n_files"],flush=True)
if __name__=="__main__":main()
