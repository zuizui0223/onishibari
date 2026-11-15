#!/usr/bin/env python3
"""PoMS 2017-2023: held-out 1-km-square comparison of local flower access
versus broad habitat / landcover in 10-minute timed flower-insect counts.

Not NDVI extraction or plant reproductive biomass. Spatial group CV at
x1km_square prohibits train/test sharing the same allocated 1-km square.
Analyzes allocated square survey only (N=3053) for primary inference;
public volunteer convenience-sample data kept separate.
No site coords or recorder identifiers exported.
"""
from __future__ import annotations
import csv,datetime as dt,hashlib,io,json,math,time,urllib.request,zipfile
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from sklearn.model_selection import GroupKFold
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error,mean_squared_error

OUT=Path("analysis_results_ukpoms_predict")
OUT.mkdir(parents=True,exist_ok=True)
UUID="d706088c-9d48-4be5-b156-6ed6ad5e0618"
URL=f"https://data-package.ceh.ac.uk/data/{UUID}.zip"
FILE="data/ukpoms_1kmfitcountdata_2017-2023.csv"
CATS_BASE=["target_flower_corrected","flower_structure","floral_unit","cloud_cover","sunshine","wind_speed"]
NUM_BASE=["month","year","hour"]
CATS_LAND=["land_cover","habitat_type"]
CATS_FLOWER=["flower_cover","flower_context"]
NUM_FLOWER=["log_floral_units"]
MODELS={
 "base_plant_season_weather":(CATS_BASE,NUM_BASE),
 "base_plus_landcover":(CATS_BASE+CATS_LAND,NUM_BASE),
 "base_plus_flower_unit_number":(CATS_BASE,NUM_BASE+NUM_FLOWER),
 "base_plus_patch_context":(CATS_BASE+CATS_FLOWER,NUM_BASE),
 "base_plus_local_flowers":(CATS_BASE+CATS_FLOWER,NUM_BASE+NUM_FLOWER),
 "landcover_plus_flower_unit_number":(CATS_BASE+CATS_LAND,NUM_BASE+NUM_FLOWER),
 "landcover_plus_patch_context":(CATS_BASE+CATS_LAND+CATS_FLOWER,NUM_BASE),
 "base_plus_both":(CATS_BASE+CATS_LAND+CATS_FLOWER,NUM_BASE+NUM_FLOWER),
}
RESPONSES=["all_insects_total","bumblebees","hoverflies","all_bees"]
SOURCE_ACK=("UK Pollinator Monitoring Scheme (UK PoMS), UK Centre for Ecology & Hydrology "
"and Joint Nature Conservation Committee / government partners, EIDC, "
"doi:10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618.")

def get_data():
    for t in range(3):
        try:
            q=urllib.request.Request(URL,headers={"User-Agent":"onishibari-PoMS-spatial-model/0.2"})
            with urllib.request.urlopen(q,timeout=100) as r:
                b=r.read(20_000_001)
            if len(b)>20_000_000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise ValueError("EIDC invalid source archive")
            with zipfile.ZipFile(io.BytesIO(b)) as z:
                raw=z.read(FILE)
            rdr=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
            return list(rdr),{"url":URL,"zip_sha256":hashlib.sha256(b).hexdigest(),
                 "csv_sha256":hashlib.sha256(raw).hexdigest(),
                 "research_acknowledgement":SOURCE_ACK}
        except OSError:
            if t==2:raise
            time.sleep(2**t)

def val(x):
    try:
        y=float(str(x).strip())
        return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None

def day(s):
    s=str(s).strip()
    for fmt in ("%Y-%m-%d","%d/%m/%Y","%d-%m-%Y"):
        try:return dt.datetime.strptime(s[:10],fmt)
        except ValueError:pass
    raise ValueError("Bad date format")

def parse_time(s):
    try:
        t=str(s).strip().split(":")
        if len(t)>=2:
            hour=int(t[0])+int(t[1])/60
            return hour if 0<=hour<24 else None
        v=float(s)
        return v if 0<=v<=24 else None
    except (ValueError,TypeError):return None

def prepare():
    raw,provenance=get_data()
    invalid=Counter()
    cleaned=[]
    for r in raw:
        total=val(r.get("all_insects_total"))
        a=[val(r.get(k)) for k in ("bumblebees","honeybees","solitary_bees")]
        h=val(r.get("hoverflies"))
        n=val(r.get("floral_unit_count"))
        site=str(r.get("x1km_square","")).strip()
        if not site:
            invalid["missing_square"]+=1;continue
        if n is None or n<=0:
            invalid["bad_flower_units"]+=1;continue
        if total is None or h is None or any(x is None for x in a):
            invalid["missing_required_outcomes"]+=1;continue
        if min([total,h,*a])<0:
            invalid["negative_count"]+=1;continue
        try:d=day(r["date"])
        except Exception:
            invalid["bad_date"]+=1;continue
        record={
            "square_id":site,
            "month":d.month,"year":d.year,
            "hour":parse_time(r.get("count_start_time")),
            "log_floral_units":math.log(n),
            "all_bees":sum(a),
        }
        for name in RESPONSES:
            if name!="all_bees":record[name]=val(r.get(name))
        for col in set(CATS_BASE+CATS_LAND+CATS_FLOWER):
            record[col]=str(r.get(col,"") or "").strip() or "__MISSING__"
        cleaned.append(record)
    return raw,cleaned,invalid,provenance

def learn(cats,nums):
    cat=Pipeline([("impute",SimpleImputer(strategy="most_frequent")),
                  ("onehot",OneHotEncoder(handle_unknown="ignore",
                      min_frequency=10, sparse_output=False))])
    numeric=Pipeline([("impute",SimpleImputer(strategy="median")),
                      ("scale",StandardScaler())])
    pre=ColumnTransformer([("cat",cat,list(cats)),("numeric",numeric,list(nums))],
                          remainder="drop",sparse_threshold=0)
    return Pipeline([("prep",pre),("ridge",Ridge(alpha=20))])

def matrix(cleaned):
    try:import pandas as pd
    except Exception as exc:raise RuntimeError("Pandas required") from exc
    return pd.DataFrame(cleaned)

def cv_experiment(df,folds,resp):
    yy=df[resp].to_numpy(dtype=float)
    # Actual response is a count in fixed 10 minute quadrat survey.
    if any(y<0 for y in yy):raise ValueError("Nonpositive count")
    yt=np.log1p(yy)
    stored={name:np.full(len(df),np.nan) for name in [*MODELS,"intercept_only"]}
    split_records=[]
    for i,(train,test) in enumerate(folds,1):
        base=float(np.mean(yt[train]))
        stored["intercept_only"][test]=np.expm1(base)
        for label,(cats,nums) in MODELS.items():
            cols=cats+nums
            model=learn(cats,nums)
            model.fit(df.iloc[train][cols], yt[train])
            pred=np.clip(np.expm1(np.clip(model.predict(df.iloc[test][cols]),-5,9)),0,10000)
            stored[label][test]=pred
        split_records.append({"fold":i,"test_samples":len(test),
           "test_unique_1km_squares":int(df.iloc[test]["square_id"].nunique()),
           "test_mean_count":float(yy[test].mean())})
    scores={}
    for label,pred in stored.items():
        if not np.isfinite(pred).all():raise ValueError("Missing CV predicted rows: "+label)
        log_mae=float(np.mean(abs(np.log1p(pred)-yt)))
        mae=float(np.mean(abs(pred-yy)))
        r2=1-float(np.sum((yt-np.log1p(pred))**2)/np.sum((yt-yt.mean())**2))
        scores[label]={"MAE_visits_per_10min":round(mae,6),
            "MAE_log1p_visits":round(log_mae,6),
            "heldout_R2_log1p_visits":round(r2,6)}
    baseline=scores["base_plant_season_weather"]["MAE_log1p_visits"]
    for label,metrics in scores.items():
        metrics["MAE_log1p_improvement_vs_base_pct"]=round(100*(baseline-metrics["MAE_log1p_visits"])/baseline,5)
    return {"response":resp,"n":len(df),"n_zero":int(sum(yy==0)),
            "folds":split_records,"scores":scores,
            "unit":"recorded insect visitor count per 10-minute 50cm-quadrat"}

def within_target_spearman(df):
    from scipy.stats import spearmanr
    z={}
    for plant,n in df.groupby("target_flower_corrected").size().sort_values(ascending=False).items():
        if n<70:continue
        subset=df[df["target_flower_corrected"]==plant]
        x=subset["log_floral_units"].to_numpy(float)
        if len(set(x))<5:continue
        entry={"records":int(n),"floral_unit_values_are_not_comparable_between_different_target_flowers":True}
        for k in ("hoverflies","bumblebees","all_insects_total"):
            y=subset[k].to_numpy(float)
            if len(set(y))<2:continue
            rho,p=spearmanr(x,y)
            entry[k]={"rho":round(float(rho),4),"p_naive":round(float(p),6)}
        z[plant]=entry
    return z

def main():
    report={"status":"started","source_doi":f"10.5285/{UUID}",
       "primary_panel":"allocated 1-km-square FIT counts, NOT self-selected public FIT counts",
       "training_validation":"5-fold GroupKFold blocking ALL records from each x1km_square into one fold",
       "models":{name:{"categorical":ca,"numeric":nu} for name,(ca,nu) in MODELS.items()},
       "measurement":"10-min fixed-area insect landings, floral-unit count and cover; NOT individual pollination effectiveness",
       "limits":[
         "Separating floral-unit count from patch context is a feature-ablation prediction exercise, not a causal decomposition. Flower count units differ by target type (individual flowers vs heads/spikes/umbels). Flower type and floral_unit are covariates; no universal per-flower equivalence.",
         "More floral units increases available opportunities in a fixed quadrat mechanically; prediction does not establish ecological scarcity or selection.",
         "Land cover is a coarse category, NOT remote NDVI or GPP/AGB measurement.",
         "Visitor detection depends on flower target choice, time of day, weather, taxonomic identification and survey eligibility.",
         "The 1 km square is the heldout unit, not the 3,053 observed visits; spatial survey location metadata never exported.",
         "No abundance of available flowers outside quadrat nor lifetime floral biomass, seed set or reproductive allocation.",
         "Observational predictions not causal effects of land cover or flower availability.",
         "Flower type 'Other' may aggregate many taxa."]}
    try:
        raw,cleaned,invalid,provenance=prepare()
        report["source"]=provenance
        report["rows_original"]=len(raw)
        report["rows_retained"]=len(cleaned)
        report["excluded_records"]=dict(invalid)
        df=matrix(cleaned)
        report["n_distinct_1km_squares"]=int(df["square_id"].nunique())
        report["years"]=dict(Counter(int(y) for y in df["year"]).most_common())
        groups=df["square_id"].to_numpy()
        folds=list(GroupKFold(n_splits=5).split(df, groups=groups))
        report["cross_validation"]={r:cv_experiment(df,folds,r) for r in RESPONSES}
        report["within_target_flower_rank_correlations"]=within_target_spearman(df)
        report["status"]="success"
    except Exception as exc:
        report["status"]="failed"
        report["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"ukpoms_1km_spatial_cv.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,default=str)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("within_target_flower_rank_correlations",)},
                     ensure_ascii=False,indent=2,default=str)[:26000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
