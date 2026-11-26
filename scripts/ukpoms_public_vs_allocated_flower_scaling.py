#!/usr/bin/env python3
"""UK PoMS standardized public-FIT vs allocated-square FIT floral scaling.

*Same monitoring scheme but two distinct sampling streams.* 13,457 public
volunteer observations are NOT a probability sample; allocated 1km-square
3,053 FITs are the primary more structured source. A robust sign across
streams is useful sensitivity, not an independent experimental replication.

Poisson pseudo-likelihood cluster robust by national grid cell, covariates
harmonized across panels excluding land_cover (not in public panel).
Flower units include flower heads, umbels and individual flowers.
No site IDs, exact coordinates or observer IDs exported.
"""
import csv,hashlib,io,json,math,time,urllib.request,zipfile
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import scipy.stats as ss

OUT=Path("analysis_results_flower_scaling_two_panels")
OUT.mkdir(exist_ok=True)
SRC="https://data-package.ceh.ac.uk/data/d706088c-9d48-4be5-b156-6ed6ad5e0618.zip"
PANEL_FILES={
 "public_volunteers":"data/ukpoms_publicfitcountdata_2017-2023.csv",
 "allocated_1km":"data/ukpoms_1kmfitcountdata_2017-2023.csv"
}
MODEL=(
 "log_floral_units + C(target_flower_corrected) + C(floral_unit)"
 " + C(flower_structure) + C(flower_cover) + C(flower_context)"
 " + C(year) + C(month) + hour + C(sunshine) + C(wind_speed)"
 " + C(cloud_cover) + C(habitat_type) + C(habitat)"
)
OUTCOMES=("all_insects_total","all_bees","hoverflies","bumblebees")

def load_source():
    for j in range(3):
        try:
            q=urllib.request.Request(SRC,headers={"User-Agent":"onishibari-floral-scaling-streams/0.2"})
            with urllib.request.urlopen(q,timeout=110) as r:blob=r.read(22_000_001)
            if len(blob)>22_000_000 or not zipfile.is_zipfile(io.BytesIO(blob)):
                raise ValueError("UK PoMS download failed")
            return blob
        except OSError:
            if j==2:raise
            time.sleep(2**j)

def number(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None

def parse_rows(raw):
    from datetime import datetime
    out=[];discard=Counter()
    rd=csv.DictReader(io.StringIO(raw.decode("utf-8-sig",errors="replace")))
    for r in rd:
        flu=number(r.get("floral_unit_count"))
        if flu is None or flu<=0:
            discard["nonpositive_missing_floral_units"]+=1;continue
        count=number(r.get("all_insects_total"))
        h=number(r.get("hoverflies"))
        b=[number(r.get(k)) for k in ("bumblebees","honeybees","solitary_bees")]
        if count is None or h is None or any(x is None for x in b):
            discard["missing_count"]+=1;continue
        if min(count,h,*b)<0:
            discard["negative_count"]+=1;continue
        square=r.get("x1km_square","")
        projection=r.get("sample_projection","")
        if not square or not projection:discard["missing_square"]+=1;continue
        date=r.get("date","")
        dateobj=None
        for fmt in ("%Y-%m-%d","%d/%m/%Y","%d-%m-%Y"):
            try:dateobj=datetime.strptime(date[:10],fmt);break
            except ValueError:pass
        if dateobj is None:discard["missing_date"]+=1;continue
        hm=r.get("count_start_time","")
        try:
            z=hm.split(":")
            hour=int(z[0])+int(z[1])/60 if len(z)>=2 else float(hm)
            if not 0<=hour<24:hour=float("nan")
        except (ValueError,TypeError):
            hour=float("nan")
        row={
            "square_id":projection+":"+square,
            "year":int(dateobj.year),"month":int(dateobj.month),
            "hour":hour,
            "floral_units":flu,
            "log_floral_units":math.log(flu),
            "all_insects_total":count,"hoverflies":h,
            "all_bees":sum(b),"bumblebees":b[0]
        }
        for field in ("target_flower_corrected","floral_unit",
                      "flower_structure","flower_cover","flower_context",
                      "sunshine","wind_speed","cloud_cover","habitat_type","habitat"):
            row[field]=str(r.get(field,"") or "").strip() or "__MISSING__"
        out.append(row)
    return pd.DataFrame(out),dict(discard)

def run(df,response,variant,family):
    D=df.copy()
    D["hour"]=D["hour"].fillna(D["hour"].median())
    X=patsy.dmatrix(MODEL,D,return_type="dataframe",NA_action="raise")
    y=D[response].to_numpy(float)
    model=sm.GLM(y,X,family=family)
    fit=model.fit(cov_type="cluster",cov_kwds={"groups":D["square_id"].to_numpy()},
                  maxiter=90,tol=1e-8)
    beta=float(fit.params["log_floral_units"])
    se=float(fit.bse["log_floral_units"])
    if not math.isfinite(beta) or not (math.isfinite(se) and se>0):
        raise ValueError("non-finite model")
    return {"variant":variant,"response":response,
      "family":family.__class__.__name__,
      "records":len(D),
      "unique_1km_grid_refs":int(D["square_id"].nunique()),
      "n_target_flower_groups":int(D["target_flower_corrected"].nunique()),
      "n_zero_insect_count":int(sum(y==0)),
      "n_visitor_count":int(sum(y)),
      "n_model_features":X.shape[1],
      "beta_ln_floral_units":round(beta,6),
      "cluster_95pct_CI":[round(beta-1.96*se,6),round(beta+1.96*se,6)],
      "p_beta_eq_0":round(float(2*ss.norm.sf(abs(beta/se))),8),
      "p_beta_eq_1":round(float(2*ss.norm.sf(abs((beta-1)/se))),8),
      "double_flowers_total_visits_multiplier":round(2**beta,6),
      "double_flowers_visits_per_unit_multiplier":round(2**(beta-1),6),
      "overdispersion_pearson":round(float(fit.pearson_chi2/fit.df_resid),5),
      "interpretation":"sub/super proportional counts per 10min fixed area, not effective pollen per individual flower"}

def main():
    report={"status":"failed","doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "premise":"Compare 2017–23 volunteer vs structured 1km survey streams for subproportional visitor count scaling",
      "models":"Poisson and fixed-alpha NB with 1km-square cluster robust errors; only shared source covariates",
      "do_not_claim":["two completely independent research studies","experimental flower manipulation",
        "measured biomass, individual floral dry mass or seed success",
        "flower unit comparability across different floral structures",
        "visit rate per flower-hour rather than observed insect landings per structural unit"],
      "privacy":"Only aggregate estimates; no raw national grid refs or recorder IDs"}
    try:
        archive=load_source()
        report["source_sha256"]=hashlib.sha256(archive).hexdigest()
        frames={}
        with zipfile.ZipFile(io.BytesIO(archive)) as z:
            for panel,name in PANEL_FILES.items():
                b=z.read(name)
                frames[panel],drop=parse_rows(b)
                report.setdefault("panels",{})[panel]={
                  "rows":len(frames[panel]),"missing_or_excluded_records":drop,
                  "source_file_sha256":hashlib.sha256(b).hexdigest(),
                  "squares":int(frames[panel]["square_id"].nunique()),
                  "target_flower_types":int(frames[panel]["target_flower_corrected"].nunique()),
                  "fraction_other_target":round(float(frames[panel]["target_flower_corrected"].str.contains("Other",case=False).mean()),5)}
        for panel,D in frames.items():
            if panel=="public_volunteers" and len(D)<12000:raise ValueError("too few public FITs")
            if panel=="allocated_1km" and len(D)<3000:raise ValueError("too few allocated FITs")
        estimates=[]
        for panel,D in frames.items():
            for y in OUTCOMES:
                est=run(D,y,panel,sm.families.Poisson())
                estimates.append(est)
                print("POISSON",panel,y,est["beta_ln_floral_units"],est["cluster_95pct_CI"],flush=True)
            for y in ("all_insects_total","all_bees","hoverflies"):
                est=run(D,y,panel,sm.families.NegativeBinomial(alpha=1))
                estimates.append(est)
                print("NB",panel,y,est["beta_ln_floral_units"],est["cluster_95pct_CI"],flush=True)
            name=panel+"_exclude_other"
            D2=D.loc[~D["target_flower_corrected"].str.contains("Other",case=False)].copy()
            estimates.append(run(D2,"all_insects_total",name,sm.families.Poisson()))
        report["estimates"]=estimates
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
        print("ERROR",report["error"],flush=True)
    (OUT/"two_poms_streams_floral_scaling.json").write_text(json.dumps(report,indent=2)+"\n")
    if report["status"]!="success":raise RuntimeError(report["error"])
if __name__=="__main__":main()
