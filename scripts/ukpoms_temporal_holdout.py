#!/usr/bin/env python3
"""Independent temporal validation of UK PoMS flower quantity and landscape.

Train 2017-2022, predict year 2023. Unlike spatial GroupKFold, locations
may recur across years. Both tests answer different generalization questions.
Same pre-specified Ridge(log1p count) model configurations as spatial pilot.
No fine site identifiers, source coordinates, collector IDs or predictions exported.
"""
import json, math
from pathlib import Path
import numpy as np
from sklearn.metrics import mean_absolute_error
from ukpoms_flower_vs_landcover_cv import prepare,matrix,RESPONSES,MODELS,learn

OUT=Path("analysis_results_ukpoms_temporal")
OUT.mkdir(exist_ok=True)

def temporal_validate():
    raw,cleaned,invalid,source=prepare()
    df=matrix(cleaned)
    train=df["year"].to_numpy()<2023
    test=df["year"].to_numpy()==2023
    if train.sum()<1000 or test.sum()<300:raise ValueError("Unexpected year split")
    report={
        "description":"Hold out ALL 2023 observations; train 2017–2022 PoMS allocated 1km FITs",
        "data_source":source,"n_train_2017_2022":int(train.sum()),
        "n_test_2023":int(test.sum()),
        "n_2023_heldout_square_ids":int(df.loc[test,"square_id"].nunique()),
        "temporal_not_geographic_holdout":True,
        "not_actual_NDVI_or_biomass":True,
        "limits":[
            "2023 shares some 1km sites with prior years; temporal generalization only",
            "Causal claims not established from flower units and habitat type",
            "The floral-unit count is not comparable across plants without type controls",
            "Visit count per 10-min fixed-area observation, not pollen transfer per single flower",
            "Year used as numerical covariate, so 2023 extrapolates one year beyond training"]}
    for target in RESPONSES:
        ytrain=np.log1p(df.loc[train,target].to_numpy(float))
        ytest=np.log1p(df.loc[test,target].to_numpy(float))
        pred_mean=float(np.mean(ytrain))
        metrics={}
        for key,features in MODELS.items():
            cats,nums=features
            cols=cats+nums
            pipeline=learn(cats,nums)
            pipeline.fit(df.loc[train,cols],ytrain)
            pred=pipeline.predict(df.loc[test,cols])
            metrics[key]={
              "mae_log1p":round(float(mean_absolute_error(ytest,pred)),6),
              "MAE_count":round(float(mean_absolute_error(np.expm1(ytest),np.clip(np.expm1(pred),0,1e4))),6),
              "heldout_R2_log1p":round(1-float(np.sum((ytest-pred)**2)/np.sum((ytest-ytest.mean())**2)),6)
            }
        metrics["intercept_only"]={
          "mae_log1p":round(float(np.mean(np.abs(ytest-pred_mean))),6)
        }
        base=metrics["base_plant_season_weather"]["mae_log1p"]
        for val in metrics.values():
            val["mae_log1p_gain_vs_base_pct"]=round((base-val["mae_log1p"])/base*100,5)
        report.setdefault("response_validation",{})[target]=metrics
    (OUT/"ukpoms_2023_heldout_validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="data_source"},indent=2)[:26000])

if __name__=="__main__":temporal_validate()
