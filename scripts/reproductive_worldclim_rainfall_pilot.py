#!/usr/bin/env python3
"""Pilot geographic climate extraction for GloPL pollen limitation experiment sites.

Observational description ONLY. WorldClim 2.1 precipitation normals = 1970-2000,
10-arcminute climatology, NOT observed rainfall in each pollen experiment year.
GloPL "Year" is PUBLICATION YEAR, do not align NDVI by that year.
RF is interpreted as rainfall, not Random Forest.
Sexual system is Wang 2021 species label; no morphological sex measurement.
NDVI, GEDI biomass, LAI or pollinator abundance are not measured here.
"""
from __future__ import annotations
import csv, hashlib, io, json, math, os, statistics, sys, tempfile, time, urllib.request, zipfile
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import rasterio
from scipy.stats import spearmanr
sys.path.insert(0,str(Path(__file__).resolve().parent))
from reproductive_data_coverage import get, find_glopl, csvrows,clean_name

OUT=Path("analysis_results_rainfall_pilot")
OUT.mkdir(exist_ok=True)
WC="https://geodata.ucdavis.edu/climate/worldclim/2_1/base/"
WANG_API="https://zenodo.org/api/records/4895147"
SOURCES=["wc2.1_10m_prec.zip","wc2.1_10m_elev.zip"]
FIELDS=["species","sex_system","latitude","longitude","family",
        "pl_experiment_cases","median_pl_effect","rf_annual_mean_mm_1970_2000",
        "rf_monthly_cv_normal_1970_2000","elevation_m_worldclim","source_year_is_publication_year"]
ACCEPT={"hermaphroditism","dioecy","monoecy"}

def fetch_zip(filename):
    url=WC+filename
    filepath=OUT/filename
    # retry with verification for interrupted servers
    for j in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"onishibari-rainfall-pilot/0.1 research"})
            with urllib.request.urlopen(req,timeout=90) as r:
                data=r.read(150_000_001)
                if len(data)>150_000_000:raise ValueError("raster exceeds expected file ceiling")
            if not zipfile.is_zipfile(io.BytesIO(data)):
                raise ValueError("not zip response")
            filepath.write_bytes(data)
            return data,{"url":url,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}
        except Exception:
            if j==2:raise
            time.sleep(2**j)

def float_valid(text):
    try:
        n=float(text)
        return n if math.isfinite(n) else None
    except Exception:return None

def get_raster_point(zips,zip_file,tif_name,coordinates):
    with zipfile.ZipFile(io.BytesIO(zips[zip_file])) as z:
        names=z.namelist()
        if tif_name not in names:
            # Monthly WorldClim archives use zero-padded month names _01 to _12.
            stem,last=tif_name.rsplit("_",1)
            month_string=last.removesuffix(".tif")
            if not month_string.isdigit():
                raise ValueError("Missing raster "+tif_name+"; actual members "+str(names[:30]))
            corrected=stem+"_"+month_string.zfill(2)+".tif"
            if corrected not in names:
                raise ValueError("Missing raster "+tif_name+" or "+corrected+"; actual members "+str(names[:30]))
            tif_name=corrected
        b=z.read(tif_name)
    tmp=OUT/tif_name
    tmp.write_bytes(b)
    try:
        with rasterio.open(tmp) as ds:
            if ds.crs.to_epsg()!=4326:raise ValueError("Unexpected raster CRS")
            v=[]
            for r in ds.sample(coordinates,masked=True):
                x=r[0]
                if np.ma.is_masked(x) or not math.isfinite(float(x)):
                    v.append(None)
                else:
                    val=float(x)
                    if ds.nodata is not None and val==float(ds.nodata): v.append(None)
                    else:v.append(val)
        return v
    finally:
        tmp.unlink(missing_ok=True)

def main():
    report={
      "dataset":"GloPL published effect sizes, Wang 2021 species labels, WorldClim 2.1 rain normals and elevation",
      "question":"Do species with different mating sex systems / pollen limitation occupy different climatological precipitation regimes?",
      "evidence":{},
      "status":"failure",
      "not_measured":["NDVI","forest AGB","plant-level biomass","flower abundance","pollinator visits or population density","functional reproductive sex","historical rainfall aligned to experimental year"],
      "warnings":["RF interpreted as rainfall; remotely sensed NDVI/AGB deliberately separate",
        "WorldClim v2.1 10 arcmin precipitation is climatological normal (1970–2000), not satellite observation and not year-matched",
        "GloPL Year field refers to publication year, not necessarily measurement year",
        "Species-level medians still represent unequal geographic sampling and shared phylogenetic effects",
        "The output is not a causal sex-system interaction effect or phylogenetic model",
        "Climatology extracted at site point with low spatial resolution, not flower-scale vegetation",
        "Raster files temporarily downloaded only; use upstream license for redistribution"]}
    try:
        wangmeta=json.loads(get(WANG_API,3_000_000)[0])
        file=next(f for f in wangmeta["files"] if f.get("key")=="SplistWithSexSyst.csv")
        wraw,_=get(file["links"]["self"],6_000_000)
        grow,src=find_glopl()
        w=csvrows(wraw)[1];g=csvrows(grow)[1]
        system={}
        for r in w:
            name=clean_name(r.get("Genus_species"))
            cat=str(r.get("SexSyst","")).lower()
            if name and cat in ACCEPT:
                if name in system and system[name]!=cat:raise ValueError("sex system conflicts")
                system[name]=cat
        rows=[]
        invalid={"lat_lon":0,"not_selected":0,"no_effect":0,"no_wang_class":0}
        for r in g:
            if r.get("Included_Bennett_et_al_meta_analysis","").lower()!="yes":
                invalid["not_selected"]+=1;continue
            name=clean_name(r.get("Species_accepted_names"))
            cat=system.get(name)
            if cat not in ACCEPT:
                invalid["no_wang_class"]+=1;continue
            val=float_valid(r.get("PL_Effect_Size"))
            if val is None:
                invalid["no_effect"]+=1;continue
            lat=float_valid(r.get("Latitude"));lon=float_valid(r.get("Longitude"))
            if lat is None or lon is None or not (-90<=lat<=90 and -180<=lon<=180):
                invalid["lat_lon"]+=1;continue
            rows.append({"species":name,"sex_system":cat,"family":r.get("Family"),
                "latitude":lat,"longitude":lon,"pl":val,"study_author":r.get("Author"),
                "publication_year":r.get("Year"),"level":r.get("Level_of_Supplementation"),
                "effect_type":r.get("PL_Effect_Size_Type1")})
        if len(rows)<200:raise ValueError(f"Unexpected low source subset {len(rows)}")
        report["input_counts"]={"GloPL_source_cases":len(g),
            "valid_matched_with_coords":len(rows),"screened_cases":invalid,
            "unique_species":len({r["species"] for r in rows})}
        points=[(r["longitude"],r["latitude"]) for r in rows]
        data={};provenance={}
        for filename in SOURCES:
            b,meta=fetch_zip(filename)
            data[filename]=b;provenance[filename]=meta
        monthly=[]
        for i in range(1,13):
            monthly.append(get_raster_point(data,SOURCES[0],f"wc2.1_10m_prec_{i}.tif",points))
        elevation=get_raster_point(data,SOURCES[1],"wc2.1_10m_elev.tif",points)
        report["source_rasters"]=provenance
        sites=defaultdict(list)
        for r in rows:
            # 6 decimal places; actual source coordinate precision can be much poorer.
            site_key=f'{r["latitude"]:.6f},{r["longitude"]:.6f}'
            sites[site_key].append(r)
        with (OUT/"site_coordinates_for_ndvi_sampling.csv").open("w",newline="") as site_file:
            writer=csv.DictWriter(site_file,fieldnames=["site_id","longitude","latitude",
                "n_effect_size_cases","year_field_is_publication_not_observation","temporally_linkable_to_modis"])
            writer.writeheader()
            for site_key,values in sorted(sites.items()):
                r=values[0]
                writer.writerow({
                    "site_id":"glopl_"+hashlib.sha256(site_key.encode()).hexdigest()[:16],
                    "longitude":r["longitude"],"latitude":r["latitude"],
                    "n_effect_size_cases":len(values),
                    "year_field_is_publication_not_observation":True,
                    "temporally_linkable_to_modis":False})
        report["distinct_coordinate_pairs_for_ndvi"]=len(sites)
        for i,r in enumerate(rows):
            mm=[arr[i] for arr in monthly]
            if any(v is None or v<0 for v in mm):
                r["annual_rf_mm"]=None;r["monthly_rain_cv"]=None
            else:
                total=sum(mm)
                r["annual_rf_mm"]=total
                r["monthly_rain_cv"]=statistics.pstdev(mm)/statistics.mean(mm) if total>0 else None
            r["elevation_m"]=elevation[i]
        grouped=defaultdict(list)
        for r in rows:
            if r["annual_rf_mm"] is not None:
                grouped[r["species"]].append(r)
        meta_rows=[]
        for name,rs in sorted(grouped.items()):
            # each species one row; no pseudo-replication of site x year cases
            rain=[x["annual_rf_mm"] for x in rs if x["annual_rf_mm"] is not None]
            cv=[x["monthly_rain_cv"] for x in rs if x["monthly_rain_cv"] is not None]
            elev=[x["elevation_m"] for x in rs if x["elevation_m"] is not None]
            if not rain:continue
            meta_rows.append({
              "species":name,"sex_system":rs[0]["sex_system"],
              "family":rs[0]["family"],
              "latitude":statistics.median(x["latitude"] for x in rs),
              "longitude":statistics.median(x["longitude"] for x in rs),
              "pl_experiment_cases":len(rs),
              "median_pl_effect":statistics.median(x["pl"] for x in rs),
              "rf_annual_mean_mm_1970_2000":statistics.median(rain),
              "rf_monthly_cv_normal_1970_2000":statistics.median(cv) if cv else None,
              "elevation_m_worldclim":statistics.median(elev) if elev else None,
              "source_year_is_publication_year":True
            })
        report["climatically_valid_species"]=len(meta_rows)
        report["no_rain_raster_cells"]=len(rows)-sum(len(x) for x in grouped.values())
        report["study_group_descriptives"]={}
        for cat in sorted(ACCEPT):
            rs=[r for r in meta_rows if r["sex_system"]==cat]
            out={"species":len(rs)}
            for col in ["median_pl_effect","rf_annual_mean_mm_1970_2000",
                "rf_monthly_cv_normal_1970_2000","elevation_m_worldclim"]:
                vals=[r[col] for r in rs if r[col] is not None]
                out[col+"_median"]=round(statistics.median(vals),4) if vals else None
            report["study_group_descriptives"][cat]=out
        report["species_level_spearman"]={}
        for cat in ["all","hermaphroditism","dioecy","monoecy"]:
            rs=meta_rows if cat=="all" else [r for r in meta_rows if r["sex_system"]==cat]
            vals={}
            for predictor in ["rf_annual_mean_mm_1970_2000",
                 "rf_monthly_cv_normal_1970_2000","elevation_m_worldclim"]:
                filtered=[r for r in rs if r[predictor] is not None]
                if len(filtered)<8:continue
                rho,pv=spearmanr([r[predictor] for r in filtered],
                                [r["median_pl_effect"] for r in filtered])
                vals[predictor]={"species":len(filtered),"spearman_rho":round(float(rho),5),
                    "p_naive":round(float(pv),7),
                    "WARNING":"P is not phylogenetically/spatially adjusted and is hypothesis generating only"}
            report["species_level_spearman"][cat]=vals
        with (OUT/"species_median_PL_worldclim_10m.csv").open("w",newline="") as fh:
            wr=csv.DictWriter(fh,fieldnames=list(meta_rows[0]))
            wr.writeheader();wr.writerows(meta_rows)
        report["status"]="success"
    except Exception as ex:
        report["error"]=type(ex).__name__+": "+str(ex)
    (OUT/"pilot_summary.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,default=str)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False,default=str)[:30000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
