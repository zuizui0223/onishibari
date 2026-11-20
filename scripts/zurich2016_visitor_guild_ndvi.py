#!/usr/bin/env python3
"""2016 Zurich garden NDVI × visitor guild composition (exploratory).

Uses published 24-garden 2016 effort-standardized observation matrix and
published pollinator taxa checklist. Reads actual MOD13Q1 V061 2016 imagery
using Planetary Computer STAC. Emits NO exact garden coordinates or site IDs.

Ecological ceiling: 24 gardens, 4 phytometer plants, no direct pollination
effectiveness or biomass; visitation counts are observation-window rates.
Do not interpret any association as sex-system evolution or causality.
"""
from collections import Counter, defaultdict
import csv, datetime as dt, hashlib, io, json, math, time, urllib.request, zipfile
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from pystac_client import Client
import planetary_computer
import rasterio
from rasterio.warp import transform

OUT=Path("analysis_results_zurich_guild")
OUT.mkdir(exist_ok=True)
SRC=("https://www.envidat.ch/dataset/d8e61f31-c2e5-4052-8784-b1624dfbfc40/"
     "resource/cef4a88c-0ef6-4327-a4a0-fc7d33b468e5/download/"
     "rejichacko_etal_2025_envidat.zip")
NON_TAXON={"Id","phytometer_plant","capture_date","capture_window"}
GUILDS=("bees","hoverflies","wasps","beetles","other")

def get_source():
    for i in range(3):
        try:
            req=urllib.request.Request(SRC,headers={"User-Agent":"Daphne-pollination-guild-audit/0.2"})
            with urllib.request.urlopen(req,timeout=100) as resp:
                b=resp.read(80_000_001)
            if len(b)>80_000_000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise ValueError("Dataset missing/not ZIP")
            return b
        except OSError:
            if i==2:raise
            time.sleep(2**i)

def csvrows(z,path):
    return list(csv.DictReader(io.StringIO(z.read(path).decode("utf-8-sig",errors="replace"))))

def norm(s):
    return " ".join(str(s or "").replace("_"," ").strip().lower().split())

def number(x):
    try:
        v=float(x)
        if math.isfinite(v):return v
    except (ValueError,TypeError):pass
    return None

def imgdate(item):
    if item.datetime is not None:return item.datetime.date()
    for k in ("start_datetime","end_datetime","datetime"):
        v=item.properties.get(k)
        if v:return dt.date.fromisoformat(v[:10])
    raise ValueError("No image date for "+item.id)

def spearman_safe(x,y):
    if len(x)<6 or len(set(x))<=1 or len(set(y))<=1:return None
    s=spearmanr(x,y)
    return {"n_gardens":len(x),"rho":round(float(s.statistic),5),
            "p_unadjusted":round(float(s.pvalue),6)}

def p_adjust_bh(records):
    ranked=sorted([(k,v["p_unadjusted"]) for k,v in records.items()
                   if v and v.get("p_unadjusted") is not None], key=lambda t:t[1])
    q=1.0
    result={}
    for i in range(len(ranked)-1,-1,-1):
        key,p=ranked[i]
        q=min(q,p*len(ranked)/(i+1))
        result[key]=round(q,6)
    return result

def canonical_group(label):
    s=norm(label)
    if "bee" in s or s in ("apoidea","anthophila"):return "bees"
    if "hover" in s or "syrph" in s:return "hoverflies"
    if "wasp" in s or "vespid" in s:return "wasps"
    if "beetle" in s or "coleoptera" in s:return "beetles"
    return "other"

def build_guild_map(checklist,taxa):
    names={}
    for r in checklist:
        token=norm(r.get("taxon"))
        if token:
            names[token]=canonical_group(r.get("pollinator_group"))
    mapping={}
    missing=[]
    for taxon in taxa:
        matched=norm(taxon)
        if matched not in names:
            missing.append(taxon)
        else:
            mapping[taxon]=names[matched]
    return mapping,missing,names

def get_ndvi(site_coordinates):
    ids=list(site_coordinates)
    xs=[site_coordinates[i][0] for i in ids]
    ys=[site_coordinates[i][1] for i in ids]
    bounds=[min(xs)-.02,min(ys)-.02,max(xs)+.02,max(ys)+.02]
    client=Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
    items=list(client.search(collections=["modis-13Q1-061"],
        bbox=bounds,datetime="2016-06-01/2016-08-01",max_items=60).items())
    terra=[i for i in items if i.id.startswith("MOD13Q1") or
       str(i.properties.get("platform","")).lower()=="terra"]
    scenes=defaultdict(dict)
    seen=set()
    for item in sorted(terra,key=lambda i:(imgdate(i),i.id)):
        t=imgdate(item)
        if t>dt.date(2016,7,20) or t+dt.timedelta(days=15)<dt.date(2016,6,15):
            continue
        if (t,item.id) in seen:continue
        seen.add((t,item.id))
        signed=planetary_computer.sign(item)
        values={}
        for name,asset in (("QA","250m_16_days_pixel_reliability"),
                           ("NDVI","250m_16_days_NDVI"),
                           ("EVI","250m_16_days_EVI")):
            with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                              GDAL_HTTP_MAX_RETRY="2"):
                with rasterio.open(signed.assets[asset].href) as ds:
                    xx,yy=transform("EPSG:4326",ds.crs,xs,ys)
                    sample=list(ds.sample(list(zip(xx,yy)),masked=True))
                    values[name]=[
                        None if np.ma.is_masked(v[0]) else float(v[0])
                        for v in sample]
        for k,sid in enumerate(ids):
            qa=values["QA"][k];ndvi=values["NDVI"][k];evi=values["EVI"][k]
            if qa not in (0,1) or ndvi is None or evi is None:continue
            if not (-2000<=ndvi<=10000 and -2000<=evi<=10000):continue
            scenes[sid][t.isoformat()]=(ndvi*.0001,evi*.0001)
    return scenes,[i.id for i in terra if dt.date(2016,6,1)<=imgdate(i)<=dt.date(2016,7,20)]

def main():
    rep={"status":"failure","doi":"10.16904/envidat.676",
         "year":2016,"garden_count_expected":24,
         "site_coordinate_values_exported":False,
         "limitations":["Visitors per recorded window, not per individual flower-hour",
           "Taxon visitor grouping ≠ actual pollination effectiveness",
           "None of the four plant species are a broad evolutionary sex-system contrast",
           "One NDVI/EVI value per garden shared across four plant species",
           "No actual individual biomass or GPP, no causal habitat effects",
           "Spatial NDVI 250 m pixels may mix neighboring cover types",
           "Site-to-site effects may be due to canopy, local floral resources or urbanization"]}
    try:
        raw=get_source()
        rep["source_sha256"]=hashlib.sha256(raw).hexdigest()
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            coords={}
            for x in csvrows(z,"03_site_data/garden_site_coordinates.csv"):
                lon=number(x["longitude"]);lat=number(x["latitude"])
                if lon is None or lat is None:raise ValueError("Missing coordinate")
                coords[norm(x["Id"])]=(lon,lat)
            visitors=csvrows(z,"06_trait_data/species_temporal_flower_visitation_matrix.csv")
            check=csvrows(z,"04_taxonomic_data/taxa_checklist.csv")
        if len(coords)!=24:raise ValueError("Expected 24 garden sites")
        taxa=[k for k in visitors[0] if k not in NON_TAXON]
        mapping,missing,names=build_guild_map(check,taxa)
        counts=Counter()
        groupcounts=defaultdict(Counter)
        bytaxon=Counter()
        windows=defaultdict(int)
        windows_zero=defaultdict(int)
        window_labels=Counter()
        survey_days=set()
        for r in visitors:
            sid=norm(r["Id"]);plant=norm(r["phytometer_plant"])
            if sid not in coords:raise ValueError("Unknown garden ID")
            x=(sid,plant)
            date=dt.datetime.strptime(r["capture_date"],"%d.%m.%Y").date()
            survey_days.add(date)
            window_labels[r["capture_window"]]+=1
            windows[x]+=1
            total=0.0
            for taxon in taxa:
                n=number(r.get(taxon))
                if n is None or n<0:continue
                if n==0:continue
                bytaxon[taxon]+=n
                counts[plant]+=n
                group=mapping.get(taxon,"other")
                groupcounts[x][group]+=n
                total+=n
            if total==0:windows_zero[x]+=1
        mapped_count=sum(v for taxon,v in bytaxon.items() if taxon in mapping)
        rep["total_visitor_events"]=sum(bytaxon.values())
        rep["taxa_available_in_matrix"]=len(taxa)
        rep["taxa_found_in_taxonomic_checklist"]=len(mapping)
        rep["taxa_unmapped_in_taxonomic_checklist"]=missing[:25]
        rep["taxa_visits_matching_checklist_fraction"]=round(mapped_count/max(1,sum(bytaxon.values())),5)
        rep["recorded_windows"]=len(visitors)
        rep["distinct_window_labels"]=len(window_labels)
        rep["survey_day_count"]=len(survey_days)
        rep["dates"]=dict(first=min(survey_days).isoformat(),last=max(survey_days).isoformat())
        rep["total_visits_by_plant"]=dict(counts)
        if rep["taxa_visits_matching_checklist_fraction"]<.98:
            raise ValueError("Taxonomic checklist does not map >=98% of observed visits")
        ndvi,scene_ids=get_ndvi(coords)
        rep["scene_ids"]=sorted(scene_ids)
        rep["gardens_with_good_ndvi"]=sum(bool(x) for x in ndvi.values())
        aggregated={}
        tests={}
        for plant in sorted(counts):
            data=[]
            for sid in coords:
                r=(sid,plant)
                n=windows.get(r,0)
                if n==0 or not ndvi.get(sid):continue
                rec=groupcounts[r]
                gpp={g:rec.get(g,0)/n for g in GUILDS}
                gtotal=sum(gpp.values())
                richness=sum(rec.get(g,0)>0 for g in GUILDS)
                t=[rec.get(g,0) for g in GUILDS if rec.get(g,0)>0]
                shannon=math.exp(-sum((v/sum(t))*math.log(v/sum(t)) for v in t)) if t else 0.0
                vi=np.array(list(ndvi[sid].values()),float)
                data.append(dict(ndvi=float(vi[:,0].mean()),evi=float(vi[:,1].mean()),
                    visits_per_window=gtotal,zero_window_share=windows_zero[r]/n,
                    pollinator_group_richness=richness,
                    effective_visiting_group_number=shannon,**gpp))
            if len(data)!=24:raise ValueError("Expected 24 gardens in "+plant)
            aggregated[plant]={
                "n_gardens":len(data),
                "visitors_total":counts[plant],
                "total_observation_windows":sum(windows[(sid,plant)] for sid in coords),
                "visitor_guild_share":{g:round(sum(groupcounts[(sid,plant)].get(g,0) for sid in coords)/counts[plant],5) for g in GUILDS},
                "visitors_zero_window_fraction":round(sum(windows_zero[(sid,plant)] for sid in coords)/sum(windows[(sid,plant)] for sid in coords),5),
                "median_pollinator_group_richness_per_garden":float(np.median([d["pollinator_group_richness"] for d in data])),
                "median_effective_group_number":round(float(np.median([d["effective_visiting_group_number"] for d in data])),4)
                }
            for metric in ["visits_per_window","zero_window_share","pollinator_group_richness","effective_visiting_group_number"]+list(GUILDS):
                test=spearman_safe([x["ndvi"] for x in data],[x[metric] for x in data])
                tests[(plant,metric)]=test
        # BH adjustment across 4 plants × 9 metrics = 36 planned NDVI tests.
        adj=p_adjust_bh(tests)
        rep["visitor_composition"]=aggregated
        rep["NDVI_vs_visitor_metrics"]={
            p:{k:dict(tests[(p,k)],q_BH_across_all_tests=adj.get((p,k)))
               if tests[(p,k)] else None
               for k in ["visits_per_window","zero_window_share","pollinator_group_richness",
                         "effective_visiting_group_number"]+list(GUILDS)}
            for p in aggregated}
        rep["top_taxa_visits"]=[{"taxon":t,"total_visits":round(v)}
            for t,v in bytaxon.most_common(10)]
        rep["status"]="success"
    except Exception as e:
        rep["error"]=type(e).__name__+": "+str(e)
    (OUT/"guild_ndvi_summary.json").write_text(json.dumps(rep,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(rep,indent=2,ensure_ascii=False)[:30000])
    if rep["status"]!="success":raise RuntimeError(rep["error"])

if __name__=="__main__":main()
