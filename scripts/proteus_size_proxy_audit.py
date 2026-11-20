#!/usr/bin/env python3
"""Exploratory size/allometry proxy audit in PROTEUS DiveRS (not biomass).

Maximumverticalheight != dry biomass; flowerSize != floral investment.
Morphologically dimorphic sexual system != strict functional dioecy.
Study sample of 364 species is not globally representative.
"""
import csv, hashlib, io, json, math, statistics, urllib.request
from collections import Counter
from pathlib import Path
from scipy.stats import spearmanr
from reproductive_data_coverage import get,csvrows,clean_name

OUT=Path("analysis_results_body_size_proxy")
OUT.mkdir(exist_ok=True)
URL="https://raw.githubusercontent.com/divers-it/rs-traitspace/main/supplementary_data/cleaned_data.csv"
FIELDS=["species","family","Woodiness","Lifespan","Mating","SexualSystem",
        "FlowerSex","Maximumverticalheight","flowerSize"]

def num(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except Exception:return None

def rankcorr(rows):
    x=[r["height"] for r in rows if r["height"] is not None and r["flower"] is not None and r["height"]>0]
    y=[r["flower"] for r in rows if r["height"] is not None and r["flower"] is not None and r["height"]>0]
    if len(x)<8:return {"n":len(x),"spearman_rho":None}
    rho,p=spearmanr(x,y)
    return {"n":len(x),"spearman_rho":round(float(rho),4),
      "naive_p":round(float(p),6),
      "note":"No phylogenetic, habitat, sampling or woodiness adjustment"}

def median_safe(values):
    v=[x for x in values if x is not None]
    return round(statistics.median(v),5) if v else None

def main():
    b,_=get(URL,limit=500_000)
    cols,raw=csvrows(b)
    if len(raw)!=364:raise ValueError("Unexpected PROTEUS cleaned n="+str(len(raw)))
    rows=[]
    for r in raw:
        species=clean_name(r.get(""))
        h=num(r.get("Maximumverticalheight"))
        f=num(r.get("flowerSize"))
        if h is not None and h<=0:h=None
        if f is not None and f<=0:f=None
        rows.append({"species":species,"height":h,"flower":f,
          "woody":r.get("Woodiness"),"sexual":r.get("SexualSystem"),
          "flowersex":r.get("FlowerSex"),"mating":r.get("Mating"),
          "lifespan":r.get("Lifespan")})
    report={"source":URL,"sha256":hashlib.sha256(b).hexdigest(),
      "n_species":len(rows),"n_height":sum(r["height"] is not None for r in rows),
      "n_flowerSize":sum(r["flower"] is not None for r in rows),
      "n_both":sum(r["height"] is not None and r["flower"] is not None for r in rows),
      "not_measured":["individual dry biomass","flower dry mass","seed biomass / floral investment","population biomass","NDVI","pollinator abundance","recruitment"],
      "allometry_spearman":{},"mating_height_groups":{},
      "flower_sex_height_groups":{},"sexual_system_height_groups":{}}
    for group in ("all","woody","herbaceous"):
        rr=rows if group=="all" else [r for r in rows if r["woody"]==group]
        report["allometry_spearman"][group]=rankcorr(rr)
    for field,out_key in [("mating","mating_height_groups"),
                          ("flowersex","flower_sex_height_groups"),
                          ("sexual","sexual_system_height_groups")]:
        for level in sorted({r[field] for r in rows if r[field]}):
            rr=[r for r in rows if r[field]==level]
            report[out_key][level]={
              "species":len(rr),
              "n_height":sum(r["height"] is not None for r in rr),
              "median_height_proxy":median_safe([r["height"] for r in rr]),
              "n_flower_size":sum(r["flower"] is not None for r in rr),
              "median_flower_size_proxy":median_safe([r["flower"] for r in rr]),
              "woody_count":sum(r["woody"]=="woody" for r in rr)}
    (OUT/"body_size_vs_reproductive_traits_descriptive.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
