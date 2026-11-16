#!/usr/bin/env python3
"""Audit two more published experiments at their original biological grain.

Silene dioica fake visual flowers / male and female reproduction:
  Dryad published v251371, doi:10.5061/dryad.0rxwdbs4d
Clarkia four species / whole-plant pollen supplementation:
  Dryad published v332601, doi:10.5061/dryad.k3j9kd56r

Only original public version ZIP; file SHA256 independently verified using the
previous source inventory. No individual, candidate-father, or plot IDs saved.
Produces source design/denominator profile FIRST; analysis must not assign
biological meanings to obscure columns before README verification.
"""
from __future__ import annotations
import csv,hashlib,io,json,urllib.request,zipfile
from collections import Counter
from pathlib import Path
import numpy as np,pandas as pd
OUT=Path("analysis_results_silene_clarkia_design_v08");OUT.mkdir(exist_ok=True)
STUDIES={
"Silene_dioica_2022":{"version":251371,"doi":"10.5061/dryad.0rxwdbs4d",
 "expected":{"allFruits_FloDis.txt":"9165c2c3c80d234b6395e91147d8f20a85fe210d0a76f28767e11b8790178bb2",
 "dataID_FloDis.txt":"1d245bb2dd6222bdaa05b0e23521eafd286cfa533b3826f8a5d8a65f715d3f9",
 "paternity_FloDis.txt":"2b440e3fad9cfd2801cba56e9de9941d13f4ee936cb9a79938d22610e62050b9",
 "README.md":"2456c5c16aec6bfc398dede02e0d99e7e3940cf607b4efdbcfa0892abfd8bfda"}},
"Clarkia_four_species_2024":{"version":332601,"doi":"10.5061/dryad.k3j9kd56r",
 "expected":{"clarkia_data_2017_kingsnake.csv":"a6b82a14761f1de282148195d9a2d56f6ba355812419253b8a2cb24014fa0568",
 "README.md":"5a6d112a0f17baaf193aca4ded1df4289855e13a3cc88e32e2bd6e4b59aac337",
 "codeToPublish_16Sep.R":"b073d870a540036d72dba4d483cb8985469c0243b37819694258129f80f60fb2"}}
}
# Correct silently? Never. Check all via pinned member SHA, reject mismatch.
def official_zip(version):
    url=f"https://datadryad.org/api/v2/versions/{version}/download"
    r=urllib.request.Request(url,headers={"User-Agent":"original-reproduction-stage-identification/0.1"})
    with urllib.request.urlopen(r,timeout=90) as res:data=res.read(750000)
    assert zipfile.is_zipfile(io.BytesIO(data))
    return zipfile.ZipFile(io.BytesIO(data))
def frame(raw):
    txt=raw.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
    return pd.read_csv(io.StringIO(txt),sep=None,engine="python",keep_default_na=True).rename(columns=lambda x:str(x).strip())
def levels(col):
    z=col.astype("string").str.strip()
    val=z.value_counts(dropna=False)
    return {str(k):int(v) for k,v in val.items()} if len(val)<=14 else None
def numeric(col):
    try:return pd.to_numeric(col,errors="coerce")
    except Exception:return pd.Series([],dtype=float)
def describe(df):
    cols=list(df.columns)
    out={"n_rows":len(df),"n_columns":len(cols),"headers":cols,"categories":{},"numerics":{},
        "id_uniqueness":{}, "missing":{}}
    for c in cols:
        ser=df[c];nmiss=int(ser.isna().sum())
        if nmiss:out["missing"][c]=nmiss
        # Do not export unique plant IDs, just counts.
        if c.lower() in ("id","plant","plant_id","plot","block","mother_id","candidate_father_id","offspring_id","trial"):
            out["id_uniqueness"][c]=int(ser.nunique(dropna=True))
            continue
        vals=levels(ser)
        if vals is not None:out["categories"][c]=vals
        if (pd.api.types.is_numeric_dtype(ser) or (len(ser)>0 and numeric(ser).notna().mean()>.9)) and c.lower() not in ["long","lat","latitude","longitude"]:
            x=numeric(ser)
            if x.notna().sum()>0:
                out["numerics"][c]={"n_nonmissing":int(x.notna().sum()),
                    "n_zero":int((x==0).sum()),"n_negative":int((x<0).sum()),
                    "median":round(float(x.median()),5),"max":round(float(x.max()),5)}
    return out
def main():
    result={"status":"failed","sources":{},
            "unit_warning":"The first pass is source audit only, no new causal fitness coefficient fitted.",
            "no_source_individual_or_plot_ids_exported":True}
    for name,info in STUDIES.items():
        z=official_zip(info["version"]);files={}
        for member in z.infolist():
            fn=Path(member.filename).name
            if fn not in info["expected"]:continue
            b=z.read(member)
            observed=hashlib.sha256(b).hexdigest()
            expected=info["expected"][fn]
            # Previous audit pins exact member SHA; stop on errors.
            if observed!=expected:
                raise ValueError(f"pinned source SHA mismatch {name}/{fn}: {observed}")
            if fn.lower().endswith((".md",".r")):
                txt=b.decode("utf-8-sig","replace")
                files[fn]={"sha256":observed,"bytes":len(b),
                   "text_excerpt_for_authorized_protocol_review":txt[:9500],
                   "text_truncated":len(txt)>9500}
            else:
                files[fn]={"sha256":observed,**describe(frame(b))}
        if set(files)!=set(info["expected"]):raise ValueError("missing source file "+name)
        result["sources"][name]={"doi":info["doi"],"version":info["version"],"files":files}
        print("SOURCE_AUDIT",name,json.dumps(result["sources"][name],ensure_ascii=False)[:28000],flush=True)
    result["status"]="source_sha_and_design_schema_verified"
    (OUT/"original_silene_clarkia_design_audit_v08.json").write_text(
         json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("SOURCE_AUDIT_COMPLETE",flush=True)
if __name__=="__main__":main()
