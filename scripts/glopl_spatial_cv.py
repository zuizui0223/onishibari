#!/usr/bin/env python3
"""Fixed-specification geographic and family out-of-sample validation.

Pilot: GloPL species-median pollen limitation vs coarse 1970-2000 WorldClim
precipitation context. Not matched to satellite NDVI nor experimental year.
No biological causality: tests whether climate/sex labels add generalizable
predictive signal beyond broad latitude/elevation.
"""
import csv,json,math,statistics
from collections import Counter
from pathlib import Path
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error,mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

OUT=Path("analysis_results_spatial_cv")
OUT.mkdir(exist_ok=True)
INPUT=Path("input/species_median_PL_worldclim_10m.csv")

def num(v):
    try:return float(v) if math.isfinite(float(v)) else float("nan")
    except (ValueError,TypeError):return float("nan")

def build(rows):
    cols=["abs_lat","elevation","log_annual_rain","rain_monthly_cv","sex_dioecy","sex_monoecy"]
    X=[]
    for r in rows:
        lat=num(r["latitude"]); elevation=num(r["elevation_m_worldclim"])
        rain=num(r["rf_annual_mean_mm_1970_2000"])
        rcv=num(r["rf_monthly_cv_normal_1970_2000"])
        X.append([abs(lat),elevation,math.log1p(rain) if rain>=0 else float("nan"),
                 rcv,int(r["sex_system"]=="dioecy"),int(r["sex_system"]=="monoecy")])
    return np.array(X,float),cols

def execute_cv(rows,holdout):
    XX,cols=build(rows)
    yy=np.array([num(r["median_pl_effect"]) for r in rows],float)
    if not np.isfinite(yy).all():raise ValueError("Nonfinite response")
    idx={c:i for i,c in enumerate(cols)}
    if holdout=="spatial_20deg":
        group=np.array([f'{math.floor(num(r["latitude"])/20)}_{math.floor(num(r["longitude"])/20)}' for r in rows])
    elif holdout=="family":
        group=np.array([r["family"] or "unknown" for r in rows])
    else:raise ValueError(holdout)
    if len(set(group))<5:raise ValueError("Not enough groups")
    designs={
      "latitude_elevation":["abs_lat","elevation"],
      "plus_rainfall":["abs_lat","elevation","log_annual_rain","rain_monthly_cv"],
      "plus_sex":["abs_lat","elevation","sex_dioecy","sex_monoecy"],
      "rainfall_and_sex":cols,
    }
    cv=GroupKFold(n_splits=5)
    folds=[]
    all_pred={k:np.full(len(rows),np.nan) for k in designs}
    for fold,(train,test) in enumerate(cv.split(XX,yy,group),1):
        entry=dict(fold=fold,test_species=len(test),
          heldout_group_count=len(set(group[test])),
          test_mean_PL=float(yy[test].mean()))
        for name,fields in designs.items():
            mat=XX[:,[idx[x] for x in fields]]
            pipe=Pipeline([("impute",SimpleImputer(strategy="median")),
                           ("scale",StandardScaler()),("ridge",Ridge(alpha=10.0))])
            pipe.fit(mat[train], yy[train])
            p=pipe.predict(mat[test])
            all_pred[name][test]=p
            entry[name+"_MAE"]=round(float(mean_absolute_error(yy[test],p)),6)
        folds.append(entry)
    summary={}
    for key,arr in all_pred.items():
        base=float(mean_absolute_error(yy,arr))
        summary[key]={
            "overall_MAE":round(base,6),
            "global_heldout_R2":round(1-float(np.sum((yy-arr)**2)/np.sum((yy-yy.mean())**2)),6),
            "fold_MAE":[f[key+"_MAE"] for f in folds]}
    for key,val in summary.items():
        val["MAE_improvement_vs_lat_elev_percent"]=round(
          (summary["latitude_elevation"]["overall_MAE"]-val["overall_MAE"])/
           summary["latitude_elevation"]["overall_MAE"]*100,4)
    return dict(holdout=holdout,n_species=len(rows),
      n_holdout_groups=len(set(group)),model="Ridge alpha 10 with median imputation and standardized numerical variables",
      folds=folds,summary=summary)

def main():
    with INPUT.open() as f:allrows=list(csv.DictReader(f))
    rows=[r for r in allrows if math.isfinite(num(r["median_pl_effect"])) and math.isfinite(num(r["latitude"])) and math.isfinite(num(r["longitude"]))]
    result={
      "scientific_status":"descriptive prediction sensitivity only; not a pollinator or biomass test",
      "source_data":"GloPL × Wang × WorldClim coarse climatology from run 43693752526",
      "n_source_species":len(allrows),
      "n_eligible_species":len(rows),
      "source_year_not_observation_year":True,
      "pre_specification":"latitude and elevation baseline; add rainfall and/or sexual system; ridge alpha=10; 5-fold spatial 20-degree and leave-families-out",
      "caveats":["GloPL mixes fruit, seed and pollen limitation outcome types",
         "Median species effect not a meta-analytic inverse-variance weighted estimate",
         "Remote NDVI and direct plant biomass NOT measured",
         "Current climate normals might not correspond to pollination survey years",
         "Phylogenetic correlations within and across families not explicitly modeled",
         "No causal environment-to-sex-system or evolutionary conclusion"]}
    result["analyses"]={}
    for kind in ["spatial_20deg","family"]:
        result["analyses"][kind]=execute_cv(rows,kind)
    (OUT/"spatial_validation.json").write_text(json.dumps(result,indent=2)+"\n")
    with (OUT/"spatial_validation_models.csv").open("w",newline="") as f:
        fieldnames=["holdout","model","n_species","mae","heldout_r2","mae_gain_pct"]
        w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader()
        for group,item in result["analyses"].items():
            for model,metrics in item["summary"].items():
                w.writerow(dict(holdout=group,model=model,n_species=item["n_species"],
                    mae=metrics["overall_MAE"],heldout_r2=metrics["global_heldout_R2"],
                    mae_gain_pct=metrics["MAE_improvement_vs_lat_elev_percent"]))
    print(json.dumps(result["analyses"],indent=2)[:16000])

if __name__=="__main__":main()
