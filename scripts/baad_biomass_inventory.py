#!/usr/bin/env python3
"""Inventory direct individual dry biomass columns in BAAD v1.0.0 (2015).

Dataset Falster et al. (2015), Ecological Archives, CC0. This stage only
inspects source fields, counts and exact version provenance; no mass
allometric comparison or sexual-system inference yet.
"""
import csv,hashlib,io,json,urllib.request,zipfile,time
from collections import Counter
from pathlib import Path
OUT=Path("analysis_results_baad")
OUT.mkdir(exist_ok=True)
URL="https://github.com/dfalster/baad/releases/download/v1.0.0/baad.zip"

def load():
    for i in range(3):
        try:
            req=urllib.request.Request(URL,headers={"User-Agent":"angiosperm-biomass-comparison/0.1"})
            with urllib.request.urlopen(req,timeout=90) as h:b=h.read(20_000_001)
            if len(b)>20_000_000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise ValueError("BAAD downloaded content invalid")
            return b
        except Exception:
            if i==2:raise
            time.sleep(2**i)
def main():
    report={"source":"Falster et al. 2015 BAAD 1.0.0 doi:10.1890/14-1889.1",
       "versioned_release":URL,"license":"CC0",
       "status":"not_started","files":[],"tables":[]}
    try:
        b=load()
        report["zip_sha256"]=hashlib.sha256(b).hexdigest()
        report["zip_size"]=len(b)
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            for f in z.infolist():
                if f.is_dir():continue
                name=f.filename
                report["files"].append({"name":name,"size":f.file_size})
                if not name.lower().endswith(".csv"):continue
                if any(w in name.lower() for w in ("contacts","references")):continue
                text=z.read(name).decode("utf8-sig",errors="replace")
                reader=csv.DictReader(io.StringIO(text))
                columns=reader.fieldnames or []
                rows=list(reader)
                profile=dict(name=name,rows=len(rows),columns=columns)
                if any(w in name.lower() for w in ("data.csv","data_")):
                    spec_col=next((k for k in columns if k.lower() in ("species","speciesname")),None)
                    if spec_col:
                        profile["distinct_species"]=len(set(r.get(spec_col) for r in rows if r.get(spec_col)))
                    numeric_cols=[x for x in columns if "mass" in x.lower() or "weight" in x.lower()]
                    profile["biomass_candidate_columns"]=numeric_cols
                    profile["nonempty_mass_candidates"]={col:sum(bool(str(r.get(col,"")).strip()) for r in rows) for col in numeric_cols}
                    profile["first_three_species"]=[r.get(spec_col) for r in rows[:3]] if spec_col else []
                report["tables"].append(profile)
        report["status"]="success"
    except Exception as e:report["error"]=type(e).__name__+": "+str(e)
    (OUT/"baad_v1_inventory.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:22000])
    if report["status"]!="success":raise RuntimeError(report["error"])
if __name__=="__main__":main()
