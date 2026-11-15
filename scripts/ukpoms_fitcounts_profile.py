#!/usr/bin/env python3
"""Independent UK PoMS 2017–2023 FIT-count data audit before any ecological model.

10-minute 0.5m×0.5m target-flower quadrat. Public participant observations
and planned 1-km square surveys are kept distinct.
No exact sample location, surveyor ID or raw rows are republished.

UK Pollinator Monitoring Scheme (2026), EIDC,
doi:10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618.
Data Open Government Licence, programme attribution required.
"""
from __future__ import annotations
import csv, datetime as dt, hashlib, io, json, math, time, urllib.request, zipfile
from collections import Counter, defaultdict
from pathlib import Path

OUT=Path("analysis_results_ukpoms_profile")
OUT.mkdir(parents=True,exist_ok=True)
UUID="d706088c-9d48-4be5-b156-6ed6ad5e0618"
SOURCE=f"https://data-package.ceh.ac.uk/data/{UUID}.zip"
DATA={
    "public":"data/ukpoms_publicfitcountdata_2017-2023.csv",
    "random_square":"data/ukpoms_1kmfitcountdata_2017-2023.csv",
}
INSECTS=("bumblebees","honeybees","solitary_bees","wasps","hoverflies",
         "other_flies","butterflies_moths","beetles","insects_small","insects_other")
META=("date","year","floral_unit_count","flower_cover","flower_context","flower_structure",
      "target_flower_corrected","target_flower_family","habitat_type","habitat",
      "land_cover","wind_speed","cloud_cover","sunshine","x1km_square")

def download():
    for n in range(3):
        try:
            req=urllib.request.Request(SOURCE,headers={"User-Agent":"onishibari-PoMS-audit/0.2"})
            with urllib.request.urlopen(req,timeout=100) as resp:
                b=resp.read(25_000_001)
            if len(b)>25_000_000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise ValueError("unexpected source archive")
            return b
        except OSError:
            if n==2:raise
            time.sleep(2**n)
def parse_number(x):
    try:
        a=float(str(x).strip())
        return a if math.isfinite(a) else None
    except (ValueError,TypeError):return None
def full_dist(x):
    return dict(Counter(str(z).strip() for z in x).most_common(60))
def profile_rows(rows,panel,cols):
    # Do not aggregate missing individual guild cells to zero for statistical tests.
    summary={
        "panel":panel,"rows":len(rows),"unique_sample_ids":len({r.get("sample_id","") for r in rows}),
        "missing_sample_ids":sum(not str(r.get("sample_id","")).strip() for r in rows),
        "duplicate_sample_id_rows":len(rows)-len({r.get("sample_id","") for r in rows}),
        "columns":cols,
    }
    summary["column_missing_fraction"]={
        k:round(sum(str(r.get(k,"")).strip().lower() in ("","na","nan","null") for r in rows)/len(rows),5)
        for k in META if k in cols
    }
    summary["category_distributions"]={
        k:full_dist(r.get(k,"") for r in rows)
        for k in ("year","country","target_flower_corrected","habitat_type",
                  "land_cover","flower_cover","flower_context","flower_structure",
                  "floral_unit","sample_projection","recorder_type") if k in cols
    }
    summary["date_check"]={"year_in_date_mismatch":0,
                           "unparseable_or_missing_date":0,
                           "year_distribution":{}}
    by_year=Counter()
    for r in rows:
        try:
            raw=str(r.get("date","")).strip()
            if not raw:raise ValueError("missing")
            # Robust ISO, British and common format parsing
            date=None
            for p in ("%Y-%m-%d","%d/%m/%Y","%d-%m-%Y","%Y/%m/%d"):
                try:date=dt.datetime.strptime(raw[:10],p).date();break
                except ValueError:pass
            if date is None:raise ValueError(raw)
            intyear=int(float(r.get("year","")))
            by_year[intyear]+=1
            if intyear!=date.year:summary["date_check"]["year_in_date_mismatch"]+=1
        except Exception:
            summary["date_check"]["unparseable_or_missing_date"]+=1
    summary["date_check"]["year_distribution"]=dict(sorted(by_year.items()))
    metrics={}
    for fld in (*INSECTS,"all_insects_total","floral_unit_count"):
        if fld not in cols:continue
        valid=[parse_number(r.get(fld)) for r in rows]
        good=[x for x in valid if x is not None]
        metrics[fld]={
            "n_valid":len(good),"n_missing":len(valid)-len(good),
            "n_zero":sum(x==0 for x in good),
            "n_negative":sum(x<0 for x in good),
            "median":sorted(good)[len(good)//2] if good else None,
            "maximum":max(good) if good else None,
            "total":sum(good)}
    summary["metric_quality"]=metrics
    differences=Counter()
    complete=0
    for r in rows:
        a=[parse_number(r.get(t)) for t in INSECTS]
        total=parse_number(r.get("all_insects_total"))
        if total is not None and all(x is not None for x in a):
            complete+=1
            diff=total-sum(a)
            differences[round(diff,8)]+=1
    summary["complete_guild_rows"]=complete
    summary["all_insects_minus_sum_guilds"]={str(k):v for k,v in differences.most_common(15)}
    summary["n_unique_1km_cells"]=len({r.get("x1km_square","") for r in rows
                                      if str(r.get("x1km_square","")).strip()}) if "x1km_square" in cols else None
    return summary

def main():
    report={
        "source":SOURCE,"doi":f"10.5285/{UUID}","status":"starting",
        "source_license":"Open Government Licence with required UK PoMS partner acknowledgement",
        "unit":"One 10-minute visit count on one 50cm × 50cm target-flower quadrat",
        "no_locations_or_recorder_ids_exported":True,
        "quality_guards":["Public convenience sample vs allocated 1km square sample are separate panels",
          "Don't impute missing visitor guild cells as zeros",
          "Check all_insects_total against sum of guild counts",
          "Check whether floral_unit_count is numeric and comparable across flower structural types",
          "Unmeasured flower pollen deposition, reproductive biomass and seed output",
          "NDVI still not sampled and must not be inferred from this file"]}
    try:
        b=download()
        report["archive_sha256"]=hashlib.sha256(b).hexdigest()
        report["archive_bytes"]=len(b)
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            report["panels"]={}
            for k,name in DATA.items():
                raw=z.read(name)
                rd=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
                rows=list(rd)
                report["panels"][k]=profile_rows(rows,k,rd.fieldnames or [])
                report["panels"][k]["file_sha256"]=hashlib.sha256(raw).hexdigest()
        report["total_cases"]=sum(p["rows"] for p in report["panels"].values())
        report["status"]="success"
    except Exception as e:report["status"]="failed";report["error"]=type(e).__name__+": "+str(e)
    (OUT/"ukpoms_coverage_and_grain.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False,default=str)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False,default=str)[:48000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
