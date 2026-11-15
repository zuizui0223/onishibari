#!/usr/bin/env python3
"""Independent PoMS 1-km survey + 2023 MODIS spatially aggregated NDVI pilot.

UK PoMS 2017-23 EIDC 10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618.
Spatial precision of x1km_square is ~1 km: average 4 locations inside each
square, NOT a single 250m pixel treated as a precise sample point.
2023 May-Aug MOD13Q1 imagery are long-season vegetation context, not exactly
contemporaneous at every survey's date. No precise locations exported.
Compare NDVI added to recorded habitat + floral units in 2023 survey subset;
cross-validate by held-out 1km squares. Not biomass or causal evidence.
"""
import datetime as dt
import json,math,re,time
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from pyproj import Transformer
from pystac_client import Client
import planetary_computer
import rasterio

from ukpoms_flower_vs_landcover_cv import prepare,matrix,learn,MODELS

OUT=Path("analysis_results_ukpoms_ndvi2023")
OUT.mkdir(exist_ok=True)
LETTERS="ABCDEFGHJKLMNOPQRSTUVWXYZ"
COLLECTION="modis-13Q1-061"
TARGET_DATES={5:dt.date(2023,5,15),6:dt.date(2023,6,15),
              7:dt.date(2023,7,15),8:dt.date(2023,8,15)}
TRANSFORMERS={
    "OSGB":Transformer.from_crs("EPSG:27700","EPSG:4326",always_xy=True),
    "OSIE":Transformer.from_crs("EPSG:29903","EPSG:4326",always_xy=True),
}
USER_AGENT="ukpoms-modis-ndvi-1km-pilot/0.1"

def grid_to_points(ref,projection):
    s=re.sub(r"\s+","",str(ref).upper())
    p=str(projection).upper()
    if p=="OSGB":
        if not re.fullmatch("[A-HJ-Z]{2}[0-9]{4}",s):
            raise ValueError("invalid British 1km grid code format")
        l1=LETTERS.index(s[0]);l2=LETTERS.index(s[1])
        e100=(((l1-2)%5)*5+(l2%5))
        n100=19-5*(l1//5)-(l2//5)
        east=e100*100000+int(s[2:4])*1000
        north=n100*100000+int(s[4:6])*1000
    elif p=="OSIE":
        if not re.fullmatch("[A-HJ-Z][0-9]{4}",s):
            raise ValueError("invalid Irish 1km grid code format")
        l=LETTERS.index(s[0])
        east=(l%5)*100000+int(s[1:3])*1000
        north=(4-l//5)*100000+int(s[3:5])*1000
    else:raise ValueError("unknown grid projection")
    return [TRANSFORMERS[p].transform(east+offx,north+offy)
             for offx,offy in [(250,250),(750,250),(250,750),(750,750)]]

def scene_date(item):
    if item.datetime is not None:return item.datetime.date()
    for fld in ("start_datetime","end_datetime","datetime"):
        x=item.properties.get(fld)
        if x:return dt.date.fromisoformat(x[:10])
    raise ValueError("MODIS item without date: "+item.id)

def scenelist():
    cat=Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
    items=list(cat.search(collections=[COLLECTION],
      bbox=[-11.0,49.0,3.0,61.8],
      datetime="2023-05-01/2023-08-31",
      max_items=500).items())
    candidates={}
    for im in items:
        if not (im.id.startswith("MOD13Q1") or str(im.properties.get("platform","")).lower()=="terra"):
            continue
        matched=re.search(r"h[0-9]{2}v[0-9]{2}",im.id,re.I)
        if not matched:continue
        date=scene_date(im)
        for month,target in TARGET_DATES.items():
            delay=abs((target-date).days)
            if delay>13:continue
            key=(matched.group(0).lower(),month)
            if key not in candidates or delay<candidates[key][0]:
                candidates[key]=(delay,im)
    selected=[(m,im) for (tile,m),(delta,im) in sorted(candidates.items())]
    if len(selected)<4:raise RuntimeError("No suitable 2023 MOD13Q1 UK scenes")
    return selected,len(items)

def sample_months(points):
    """Return only square-key/masked image aggregate values, no lat/lon exported."""
    sites=list(points)
    all_geo=[pt for site in sites for pt in points[site]]
    selected,total_stac=scenelist()
    out=defaultdict(list)
    audit=[]
    for month,item in selected:
        signed=planetary_computer.sign(item)
        qa_url=signed.assets["250m_16_days_pixel_reliability"].href
        nd_url=signed.assets["250m_16_days_NDVI"].href
        validrows=[]
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",GDAL_HTTP_MAX_RETRY="2"):
            with rasterio.open(qa_url) as qa:
                xs,ys=rasterio.warp.transform("EPSG:4326",qa.crs,
                     [x for x,y in all_geo],[y for x,y in all_geo])
                eligible=[i for i,(x,y) in enumerate(zip(xs,ys))
                    if qa.bounds.left<=x<qa.bounds.right and qa.bounds.bottom<=y<qa.bounds.top]
                if not eligible:continue
                qvals=list(qa.sample([(xs[i],ys[i]) for i in eligible],masked=True))
            with rasterio.open(nd_url) as nd:
                nvals=list(nd.sample([(xs[i],ys[i]) for i in eligible],masked=True))
        nvalid=0
        for idx,q,n in zip(eligible,qvals,nvals):
            q=q[0];n=n[0]
            if np.ma.is_masked(q) or np.ma.is_masked(n):continue
            if int(q) not in (0,1):continue
            if not (-2000<=float(n)<=10000):continue
            sid=sites[idx//4]
            out[(sid,month)].append(float(n)*.0001)
            nvalid+=1
        audit.append({"image_id":item.id,"month":month,
           "source_date":scene_date(item).isoformat(),
           "n_valid_subpixel_readings":nvalid})
    summarized={}
    for s in sites:
        vals={}
        for m in TARGET_DATES:
            p=out.get((s,m),[])
            if p:vals[m]=float(np.median(p))
        if len(vals)>=2:
            summarized[s]={
                "ndvi_mean_may_aug_2023":float(np.mean(list(vals.values()))),
                "ndvi_amplitude_4month_sample":float(max(vals.values())-min(vals.values())),
                "n_valid_months":len(vals)}
    return summarized,audit,total_stac

def fixed_validation(df):
    from sklearn.model_selection import GroupKFold
    predictors_cat,predictors_num=MODELS["base_plus_both"]
    all_features=list(predictors_cat)+list(predictors_num)
    newfeatures=all_features+["ndvi_mean_may_aug_2023","ndvi_amplitude_4month_sample"]
    folds=list(GroupKFold(n_splits=5).split(df,groups=df["square_id"].to_numpy()))
    outcomes={}
    for target in ("all_insects_total","hoverflies","all_bees","bumblebees"):
        obs=np.log1p(df[target].to_numpy(float))
        saved={}
        for name,columns,numeric in [
            ("floral_plus_landcover",all_features,list(predictors_num)),
            ("floral_landcover_plus_NDVI",newfeatures,list(predictors_num)+
                 ["ndvi_mean_may_aug_2023","ndvi_amplitude_4month_sample"])]:
            fit=np.full(len(df),np.nan)
            for train,test in folds:
                model=learn(predictors_cat,numeric)
                model.fit(df.iloc[train][columns],obs[train])
                fit[test]=model.predict(df.iloc[test][columns])
            saved[name]={
                "mae_log1p":float(np.mean(abs(obs-fit))),
                "heldout_R2_log1p":float(1-np.sum((obs-fit)**2)/np.sum((obs-obs.mean())**2)),
            }
        ref=saved["floral_plus_landcover"]["mae_log1p"]
        cmp=saved["floral_landcover_plus_NDVI"]["mae_log1p"]
        saved["gain_pct_from_adding_2023_NDVI_context"]=round((ref-cmp)/ref*100,5)
        for i in saved:
            if isinstance(saved[i],dict):
                for k in saved[i]:saved[i][k]=round(saved[i][k],6)
        outcomes[target]=saved
    return outcomes

def main():
    report={"status":"failed",
      "ukpoms_doi":"10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618",
      "satellite":"NASA MOD13Q1 V061 QA-masked 250m NDVI, 2023 May–August",
      "spatial_unit":"1 km national grid squares, sample 4 quadrants (0.25,0.75km centers) and average, never treat as GPS precise",
      "privacy":"No site grid references or computed coordinates in output",
      "limits":[
        "Only 2023 PoMS allocated square surveys; not an independent 7-year panel model",
        "Four MODIS composites nearest 15th of May–August are NOT exact date-specific flowering exposure",
        "Satellite greenness includes all vegetation; NOT actual floral mass, understory light or pollinators",
        "NDVI 2023 is not a causal predictor of within-season plant sex-system evolution",
        "Spatial 1km squares may contain mixtures of land cover and unknown quadrat microhabitat",
        "Model response is ten-minute insect count, not per-flower pollen deposition or seed production",
        "Limited NDVI contrast within 101 UK sites may be collinear with recorded land cover"]}
    try:
        raw,cleaned,invalid,source=prepare()
        square_meta={}
        for r in raw:
            key=r.get("x1km_square","")
            proj=r.get("sample_projection","")
            if key in square_meta and square_meta[key]!=proj:
                raise ValueError("same grid code used with conflicting projection")
            square_meta[key]=proj
        points={s:grid_to_points(s,p) for s,p in square_meta.items()}
        all_lat=[lat for cells in points.values() for lon,lat in cells]
        all_lon=[lon for cells in points.values() for lon,lat in cells]
        if not (49<=min(all_lat)<=max(all_lat)<61.7 and -11<=min(all_lon)<=max(all_lon)<3):
            raise ValueError("Grid code conversion produced implausible UK/Irish coordinates")
        report["geocoded_square_count"]=len(points)
        matched,scenes,stac_count=sample_months(points)
        report["STAC_2023_scene_candidates"]=stac_count
        report["n_selected_tiles_dates_with_valid_pixels"]=len(scenes)
        report["n_squares_with_2plus_valid_months"]=len(matched)
        report["QA_scene_samples"]=scenes
        df=matrix(cleaned)
        df=df[(df["year"]==2023)&df["square_id"].isin(matched)].copy()
        report["2023_observations_with_NDVI"]=len(df)
        report["2023_squares_with_NDVI_and_FIT"]=int(df["square_id"].nunique())
        if len(df)<200 or df["square_id"].nunique()<40:
            raise ValueError("Insufficient QC NDVI-matched 2023 sample for spatial CV")
        for name in ["ndvi_mean_may_aug_2023","ndvi_amplitude_4month_sample"]:
            df[name]=[matched[x][name] for x in df["square_id"]]
        report["NDVI_group_median"]=round(float(np.median(df["ndvi_mean_may_aug_2023"])),5)
        report["heldout_1km_square_predictions"]=fixed_validation(df)
        report["status"]="success"
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"ukpoms_2023_modis_1km_pilot.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2)[:20000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
