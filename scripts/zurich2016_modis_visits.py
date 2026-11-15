#!/usr/bin/env python3
"""Spatially and temporally matched 2016 MODIS greenness versus visits.

Chacko, Moretti & Frey (2025): actual survey 2016-06-15--07-20,
24 Zurich private gardens, four phytometer plants.
Samples NASA MOD13Q1 250m Terra 16-day NDVI/EVI + QA. No exact private
coordinates or garden IDs are emitted. Does NOT infer pollination efficacy,
flower biomass, or cause of species-level mating system.
"""
from collections import defaultdict,Counter
import csv,datetime as dt,hashlib,io,json,math,time,urllib.request,zipfile
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from pystac_client import Client
import planetary_computer
import rasterio
from rasterio.warp import transform

OUT=Path("analysis_results_zurich_2016_modis")
OUT.mkdir(exist_ok=True)
ZIPURL=("https://www.envidat.ch/dataset/d8e61f31-c2e5-4052-8784-b1624dfbfc40/"
 "resource/cef4a88c-0ef6-4327-a4a0-fc7d33b468e5/download/"
 "rejichacko_etal_2025_envidat.zip")
COL="modis-13Q1-061"
SITE_FILE="03_site_data/garden_site_coordinates.csv"
VISIT_FILE="06_trait_data/species_temporal_flower_visitation_matrix.csv"
ASSETS={"NDVI":"250m_16_days_NDVI","EVI":"250m_16_days_EVI","QA":"250m_16_days_pixel_reliability"}
NON_VISIT={"Id","phytometer_plant","capture_date","capture_window"}

def fetch():
    for attempt in range(3):
        try:
            q=urllib.request.Request(ZIPURL,headers={"User-Agent":"Daphne-interaction-environment/0.2 research"})
            with urllib.request.urlopen(q,timeout=100) as r:data=r.read(80_000_001)
            if len(data)>80_000_000 or not zipfile.is_zipfile(io.BytesIO(data)):
                raise ValueError("bad data archive")
            return data
        except OSError:
            if attempt==2:raise
            time.sleep(2**attempt)
def csvfrom(z,name):
    return list(csv.DictReader(io.StringIO(z.read(name).decode("utf-8-sig",errors="replace"))))
def number(s):
    try:
        x=float(s)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def date(i):
    if i.datetime is not None:return i.datetime.date()
    for key in ("start_datetime","end_datetime","datetime"):
        val=i.properties.get(key)
        if val:return dt.date.fromisoformat(val[:10])
    raise ValueError("No scene date "+i.id)
def spearman_safe(x,y):
    if len(x)<6 or len(set(x))<2 or len(set(y))<2:return None
    s=spearmanr(x,y)
    return {"n_gardens":len(x),"rho":round(float(s.statistic),5),
      "p_unadjusted":round(float(s.pvalue),6)}
def main():
    report={"source_zip":ZIPURL,"study_DOI":"10.16904/envidat.676",
      "satellite_product":"NASA MOD13Q1 V061 Terra 250m 16-day pixel reliability and NDVI EVI",
      "year_matched":2016,"site_count_expected":24,
      "coordinates_not_exported":True,"status":"failed","warnings":[
        "One 250m greenness pixel is not garden biomass or floral biomass",
        "Pollinator catch events per matrix window not flowers per hour",
        "Within-window number of events includes potential nectar robbers / inefficient visitors",
        "Observational site associations, N=24 sites, no urbanization/canopy/flower density adjustment",
        "Spatial coarse pixel mixing and QA masking make true NDVI uncertainty larger",
        "Cannot generalize to sexual system evolution from four phytometer plant species"]}
    try:
        raw=fetch()
        report["source_zip_sha256"]=hashlib.sha256(raw).hexdigest()
        with zipfile.ZipFile(io.BytesIO(raw)) as zipf:
            site=csvfrom(zipf,SITE_FILE)
            vis=csvfrom(zipf,VISIT_FILE)
            # Four distinct reproductive endpoints per standardized phytometer.
            # No garden coordinates, plant IDs or site-specific outcomes exported.
            endpoint_sources={
                "carrot":("07_pollination_success/daucus_carota_seed_set.csv","n_seeds",None),
                "sainfoin":("07_pollination_success/onobrychis_viciifolia_fruit_set.csv","n_flowers_with_fruits","n_flowers_without_fruits"),
                "radish":("07_pollination_success/raphanus_sativus_fruit_set.csv","n_flowers_with_fruits","n_flowers_without_fruits"),
                "comfrey":("07_pollination_success/symphytum_officinale_fruit_set.csv","n_flowers_with_seeds","n_flowers_without_seeds")}
            reproductive={}
            for plant,(file,success_col,fail_col) in endpoint_sources.items():
                buckets=defaultdict(lambda:[0.0,0.0])
                for rec in csvfrom(zipf,file):
                    excluded=str(rec.get("exclude") or rec.get("excluded") or "").strip().lower()
                    if excluded in ("1","yes","true","excluded","exclude"):continue
                    garden=str(rec.get("Id","")).strip().lower()
                    success=number(rec.get(success_col))
                    if success is None or success<0:continue
                    if fail_col is None:
                        numerator,denominator=success,1.0
                    else:
                        failures=number(rec.get(fail_col))
                        if failures is None or failures<0 or success+failures<=0:continue
                        numerator,denominator=success,success+failures
                    buckets[garden][0]+=numerator
                    buckets[garden][1]+=denominator
                reproductive[plant]={g:y[0]/y[1] for g,y in buckets.items() if y[1]>0}
            coords={}
            for r in site:
                key=str(r["Id"]).strip().lower()
                x=number(r["longitude"]);y=number(r["latitude"])
                if x is None or y is None:raise ValueError("missing site coords")
                coords[key]=(x,y)
            if len(coords)!=24:raise ValueError("Expected 24 garden coordinates; found "+str(len(coords)))
            # visitor sampling dates and counts per recorded window, excluding no visitor as zeros.
            taxa=[k for k in vis[0] if k not in NON_VISIT]
            visits=defaultdict(list)
            days=set()
            for r in vis:
                d=dt.datetime.strptime(r["capture_date"],"%d.%m.%Y").date()
                days.add(d)
                site_id=str(r["Id"]).strip().lower()
                plant=r["phytometer_plant"].lower()
                values=[number(r[k]) for k in taxa]
                n=sum(v for v in values if v is not None and v>0)
                visits[(site_id,plant)].append(n)
            min_day=min(days);max_day=max(days)
            if not (dt.date(2016,6,1)<=min_day<=dt.date(2016,6,30) and
                    dt.date(2016,7,1)<=max_day<=dt.date(2016,8,1)):
                raise ValueError("Unexpected actual field dates")
            report["field_start"]=min_day.isoformat()
            report["field_end"]=max_day.isoformat()
            report["n_field_dates"]=len(days)
            report["n_visitor_windows"]=len(vis)
            report["plant_species"]=sorted(set(k[1] for k in visits))
            xs=[p[0] for p in coords.values()];ys=[p[1] for p in coords.values()]
            bbox=[min(xs)-.02,min(ys)-.02,max(xs)+.02,max(ys)+.02]
            # 4-5 MODIS dates covering flowering field period
            cat=Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
            items=list(cat.search(collections=[COL],bbox=bbox,
                datetime="2016-06-01/2016-08-01",max_items=80).items())
            terra=[i for i in items if i.id.startswith("MOD13Q1") or
               str(i.properties.get("platform")).lower()=="terra"]
            # One scene per date and tile: if border crosses MODIS tiles,
            # each site will sample all with QA and keep at most one valid scene per date.
            scene_ids=set()
            observations=defaultdict(dict)
            image_counts={}
            for item in sorted(terra,key=lambda i:(date(i),i.id)):
                # A MODIS composite extends about 16 days from its start.
                # Include only scenes overlapping the actual June15-July20 field season.
                if date(item)>max_day or date(item)+dt.timedelta(days=15)<min_day:
                    continue
                key=(str(date(item)),item.id)
                if key in scene_ids:continue
                scene_ids.add(key)
                signed=planetary_computer.sign(item)
                vals={}
                for var,asset in ASSETS.items():
                    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                                      GDAL_HTTP_MAX_RETRY="2"):
                        with rasterio.open(signed.assets[asset].href) as ds:
                            xy=transform("EPSG:4326",ds.crs,xs,ys)
                            if var=="NDVI":
                                report.setdefault("unique_250m_pixels_by_scene",{})[item.id]=len(set(
                                   ds.index(x,y) for x,y in zip(xy[0],xy[1])))
                            data=[a[0] for a in ds.sample(list(zip(xy[0],xy[1])),masked=True)]
                            vals[var]=[(None if np.ma.is_masked(v) else float(v)) for v in data]
                nvalid=0
                for ix,siteid in enumerate(coords):
                    qa=vals["QA"][ix]
                    ndvi=vals["NDVI"][ix];evi=vals["EVI"][ix]
                    if (qa in (0,1) and ndvi is not None and evi is not None and
                            -2000<=ndvi<=10000 and -2000<=evi<=10000):
                        observations[siteid][str(date(item))]=(ndvi*.0001,evi*.0001)
                        nvalid+=1
                image_counts[item.id]={"scene_start_date":str(date(item)),
                                        "n_garden_pixels_accepted":nvalid}
            report["source_scene_validity"]=image_counts
            report["n_images"]=len(image_counts)
            report["sites_with_any_good_ndvi"]=sum(bool(v) for v in observations.values())
            # No garden IDs/coordinates output; only group-level correlations.
            summary={}
            for plant in sorted({p for (siteid,p) in visits}):
                xx=[];yy=[];aa=[]
                ndvi_repro=[]; evi_repro=[]; visits_repro=[]; repro_out=[]
                for sid in coords:
                    visits_raw=visits.get((sid,plant))
                    scenes=list(observations.get(sid,{}).values())
                    if not visits_raw or not scenes:continue
                    ndvi=sum(v[0] for v in scenes)/len(scenes)
                    evi=sum(v[1] for v in scenes)/len(scenes)
                    visitors=sum(visits_raw)/len(visits_raw)
                    xx.append(ndvi);aa.append(evi)
                    yy.append(visitors)
                    if sid in reproductive.get(plant,{}):
                        ndvi_repro.append(ndvi)
                        evi_repro.append(evi)
                        visits_repro.append(visitors)
                        repro_out.append(reproductive[plant][sid])
                summary[plant]={
                  "n_gardens_qa_matched":len(xx),
                  "ndvi_median":round(float(np.median(xx)),5) if xx else None,
                  "visits_per_recorded_window_median":round(float(np.median(yy)),5) if yy else None,
                  "NDVI_vs_visit_count":spearman_safe(xx,yy),
                  "EVI_vs_visit_count":spearman_safe(aa,yy),
                  "reproduction_endpoint":("mean seeds per umbel" if plant=="carrot" else
                      "fruit set proportion" if plant in ("radish","sainfoin") else
                      "fraction flowers with seeds"),
                  "n_gardens_reproduction_matched":len(repro_out),
                  "NDVI_vs_reproductive_outcome":spearman_safe(ndvi_repro,repro_out),
                  "EVI_vs_reproductive_outcome":spearman_safe(evi_repro,repro_out),
                  "visitors_vs_reproductive_outcome":spearman_safe(visits_repro,repro_out)}
            # Exploratory eight planned tests (4 plants × 2 vegetation indices).
            # All shares the same 24 garden locations and tests are correlated.
            eligible=[]
            for p,x in summary.items():
                for index in ("NDVI","EVI"):
                    rec=x.get(index+"_vs_visit_count")
                    if rec and rec.get("p_unadjusted") is not None:
                        eligible.append((p,index,rec["p_unadjusted"]))
            ordered=sorted(eligible,key=lambda x:x[2])
            qvals={}
            qnext=1.0
            for idx in range(len(ordered)-1,-1,-1):
                pl,vi,pv=ordered[idx]
                qnext=min(qnext,pv*len(ordered)/(idx+1))
                qvals[(pl,vi)]=round(qnext,6)
            for p,x in summary.items():
                for ix in ("NDVI","EVI"):
                    z=x.get(ix+"_vs_visit_count")
                    if z:z["q_BH_across_eight_tests"]=qvals.get((p,ix))
            report["multiple_testing_method"]="Benjamini-Hochberg across 4 plants × 2 indices for visits association; study site predictors not independent"
            # Additional eight tests are a different outcome family. They are
            # exploratory and cannot establish mediation or adaptation.
            obs=[]
            for p, x in summary.items():
                for index in ("NDVI","EVI"):
                    rec=x.get(index+"_vs_reproductive_outcome")
                    if rec and rec.get("p_unadjusted") is not None:
                        obs.append((p,index,rec["p_unadjusted"]))
            ranked=sorted(obs,key=lambda x:x[2])
            qnext=1.0
            for j in range(len(ranked)-1,-1,-1):
                plant,index,pval=ranked[j]
                qnext=min(qnext,pval*len(ranked)/(j+1))
                summary[plant][index+"_vs_reproductive_outcome"]["q_BH_across_eight_repro_tests"]=round(qnext,6)
            report["additional_endpoint_tests"]="NDVI/EVI vs reproductive endpoint separately BH-adjusted across eight tests; visitor-reproductive correlations remain descriptive"
            report["by_plant"]=summary
            report["status"]="success" if report["sites_with_any_good_ndvi"]>=18 else "insufficient"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
    (OUT/"zurich2016_ndvi_visitation_summary.json").write_text(
      json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:30000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
