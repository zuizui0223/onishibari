#!/usr/bin/env python3
"""Real MOD13Q1 V061 NDVI/EVI pixel probe at GloPL-referenced localities.

This tests remote sensing access and QA, NOT any ecological hypothesis. Dates
2013-02/08 are arbitrary habitat snapshots, NOT pollen experiment dates.
"""
from collections import defaultdict
import csv, datetime as dt, hashlib, json, math, os, time
from pathlib import Path
import planetary_computer
from pystac_client import Client
import rasterio
from rasterio.warp import transform

OUT=Path("analysis_results_modis_probe")
OUT.mkdir(exist_ok=True)
SOURCE=Path("input/species_median_PL_worldclim_10m.csv")
GROUPS=("dioecy","monoecy","hermaphroditism")
DATES=("2013-02-01","2013-08-01")
STAC="https://planetarycomputer.microsoft.com/api/stac/v1"
ASSETS={"NDVI":"250m_16_days_NDVI",
        "EVI":"250m_16_days_EVI",
        "QA":"250m_16_days_pixel_reliability"}

def sites():
    with SOURCE.open() as f:rows=list(csv.DictReader(f))
    by=defaultdict(list)
    for r in rows:
        try: lat=float(r["latitude"]);lon=float(r["longitude"])
        except (ValueError,KeyError):continue
        if r["sex_system"] in GROUPS and -60<=lat<=60 and -180<=lon<=180:
            r["_lat"]=lat;r["_lon"]=lon
            by[r["sex_system"]].append(r)
    chosen=[]
    for group in GROUPS:
        ranked=sorted(by[group],key=lambda r:hashlib.sha256(r["species"].encode()).hexdigest())
        if len(ranked)<3:raise ValueError("Few sites for "+group)
        chosen.extend(ranked[:3])
    return chosen,len(rows),{g:len(by[g]) for g in GROUPS}

def extract(item,asset,lon,lat):
    href=planetary_computer.sign(item).assets[asset].href
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      GDAL_HTTP_MAX_RETRY="2"):
        with rasterio.open(href) as ds:
            x,y=transform("EPSG:4326",ds.crs,[lon],[lat])
            a=next(ds.sample([(x[0],y[0])],masked=True))[0]
            return None if bool(getattr(a,"mask",False)) else int(a)

def observe(client,r,day):
    date=dt.date.fromisoformat(day)
    start=(date-dt.timedelta(days=21)).isoformat()
    end=(date+dt.timedelta(days=21)).isoformat()
    lon=r["_lon"];lat=r["_lat"]
    items=list(client.search(collections=["modis-13Q1-061"],
      bbox=[lon-0.015,lat-0.015,lon+0.015,lat+0.015],
      datetime=start+"/"+end,max_items=32).items())
    if not items:raise RuntimeError("STAC has no items for this point/date")
    def rank(i):
        terra=0 if i.id.startswith("MOD13Q1") or str(i.properties.get("platform")).lower()=="terra" else 1
        gap=abs((i.datetime.date()-date).days) if i.datetime else 1000
        return terra,gap,i.id
    item=min(items,key=rank)
    values={v:extract(item,asset,lon,lat) for v,asset in ASSETS.items()}
    good=values["QA"] in [0,1]
    return dict(species=r["species"],sex_system=r["sex_system"],
        GloPL_species_median_PL=r["median_pl_effect"],latitude=lat,longitude=lon,
        requested_date=day,image_date=str(item.datetime.date()),
        stac_item=item.id,platform=item.properties.get("platform"),
        QA_reliability=values["QA"],QA_acceptable=good,
        NDVI=(values["NDVI"]*0.0001 if good and values["NDVI"] is not None and -2000<=values["NDVI"]<=10000 else None),
        EVI=(values["EVI"]*0.0001 if good and values["EVI"] is not None and -2000<=values["EVI"]<=10000 else None),
        contemporaneous_with_pollination_experiment=False)

def main():
    selected,n,eligible=sites()
    report=dict(purpose="MODIS actual-pixel retrieval feasibility",
        source="NASA MOD13Q1 V061, Microsoft Planetary Computer STAC signed COG",
        source_species=n,eligible_by_group=eligible,n_species_selected=len(selected),
        requested_dates=list(DATES),n_attempts=len(selected)*len(DATES),
        date_matched_to_experiment=False,
        no_direct_biomass_or_pollinator_measurement=True,successful_samples=0,
        QA_accepted_samples=0,errors=[],status="started")
    good=[]
    catalog=Client.open(STAC)
    for r in selected:
        for date in DATES:
            for attempt in (0,1):
                try:
                    row=observe(catalog,r,date)
                    good.append(row)
                    print("SAMPLE",r["species"],date,row["stac_item"],row["NDVI"],row["QA_reliability"],flush=True)
                    break
                except Exception as ex:
                    if attempt==1:
                        report["errors"].append(dict(species=r["species"],date=date,
                            error=type(ex).__name__+": "+str(ex)[:300]))
                        print("FAILED",r["species"],date,type(ex).__name__,str(ex)[:150],flush=True)
                    else:time.sleep(2)
    report["successful_samples"]=len(good)
    report["QA_accepted_samples"]=sum(r["NDVI"] is not None and r["EVI"] is not None for r in good)
    if good:
        with (OUT/"modis_points.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(good[0]));w.writeheader();w.writerows(good)
    report["status"]="success" if report["QA_accepted_samples"]>=3 else "insufficient"
    (OUT/"summary.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,indent=2))
    if report["status"]!="success":raise RuntimeError("Point sampling failed or insufficient QA")

if __name__=="__main__":main()
