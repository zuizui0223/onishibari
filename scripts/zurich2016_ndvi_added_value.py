#!/usr/bin/env python3
"""Does observed MODIS greenness add to the *published* garden habitat covariates?

Original experiment: 24 Zurich gardens in 2016, 4 phytometer species.
Predictor response: hoverfly visitation events per *recorded observation window*;
not per flower per hour or true effective pollen transfer.

Controlling study-published local floral species richness and urban impervious
area addresses a major confound in raw NDVI correlations. Descriptive only;
post hoc follow-up of 2016 carrot hoverfly greenness signal, not confirmation.
No exact garden positions or identifiers are written to results.
"""
from collections import defaultdict
import csv,hashlib,io,json,math,time,urllib.request
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneOut

from zurich2016_visitor_guild_ndvi import (
    get_source,csvrows,norm,number,build_guild_map,get_ndvi,NON_TAXON)

OUT=Path("analysis_results_zurich_ndvi_incremental")
OUT.mkdir(exist_ok=True)
BASE="https://raw.githubusercontent.com/mrejichacko/BetterBlooms/main/raw_data/"
COV1="explanatory_variables.txt"
COV2="plant_floristic_data.txt"

def read_github_semicolon(name):
    uri=BASE+name
    for t in range(3):
        try:
            req=urllib.request.Request(uri,headers={"User-Agent":"Daphne-garden-environment-audit/0.4"})
            with urllib.request.urlopen(req,timeout=45) as r:b=r.read(150000)
            data=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig")),delimiter=";"))
            if len(data)!=24:raise ValueError(f"Expected 24 gardens for {name}, got {len(data)}")
            return data,{"url":uri,"sha256":hashlib.sha256(b).hexdigest(),"rows":len(data)}
        except OSError:
            if t==2:raise
            time.sleep(2**t)

def partial_rank(x,y,cov):
    xx=rankdata(x)
    yy=rankdata(y)
    Z=np.asarray(cov,dtype=float)
    Z=(Z-Z.mean(axis=0))/(Z.std(axis=0)+1e-12)
    C=np.column_stack((np.ones(len(Z)),Z))
    xres=xx-C @ np.linalg.lstsq(C,xx,rcond=None)[0]
    yres=yy-C @ np.linalg.lstsq(C,yy,rcond=None)[0]
    return float(np.corrcoef(xres,yres)[0,1])

def spearman(x,y):
    if len(set(x))<=1 or len(set(y))<=1:return None
    z=spearmanr(x,y)
    return {"rho":round(float(z.statistic),5),"p_naive":round(float(z.pvalue),6)}

def cv_models(x,y):
    y=np.asarray(y,dtype=float)
    X=np.asarray(x,dtype=float)
    if X.ndim!=2:raise ValueError("Bad design")
    p=np.zeros(len(y))
    for train,test in LeaveOneOut().split(X):
        model=Pipeline([("z",StandardScaler()),("ridge",Ridge(alpha=5.0))])
        model.fit(X[train],y[train])
        p[test]=model.predict(X[test])
    return {"leave_one_garden_out_MAE":round(float(np.mean(np.abs(y-p))),6),
      "leave_one_garden_out_R2":round(1-float(np.sum((y-p)**2)/np.sum((y-y.mean())**2)),5)}

def main():
    report={"status":"failed","primary_target":"carrot hoverfly visits per observation window",
        "source_field_year":2016,"N_gardens_expected":24,
        "focus":"Does 2016 NDVI/EVI add predictive value beyond the same-study 500m impervious surface and local floral richness?",
        "prior_research_DOI":"10.1111/1365-2664.70384",
        "warnings":[
          "The authors already showed urban densification and flower-richness differences in the SAME original 24-garden experiment (2026); this audit assesses remote-sensing proxies, not independent novelty or validation",
          "Post hoc selected NDVI-hoverfly correlation; partial ranks or LOOCV do not turn this into causal inference",
          "All landscape predictors are garden-level, 24 effective independent environments; do not count 22 survey days as landscape replicates",
          "Neither 250m NDVI nor urban cover directly measures floral biomass, flower availability, or pollen deposition",
          "Linear adjustment for urban 500m + local floral richness, no spatial covariance model",
          "Very small N for detecting additive predictive value; negative LOOCV R2 can indicate no usable predictive skill",
          "Precise private garden site identifiers and coordinates omitted from outputs"]}
    try:
        raw=get_source()
        report["original_data_sha256"]=hashlib.sha256(raw).hexdigest()
        with __import__("zipfile").ZipFile(io.BytesIO(raw)) as z:
            coords={}
            for rec in csvrows(z,"03_site_data/garden_site_coordinates.csv"):
                coords[norm(rec["Id"])]=(number(rec["longitude"]),number(rec["latitude"]))
            visits=csvrows(z,"06_trait_data/species_temporal_flower_visitation_matrix.csv")
            taxa=csvrows(z,"04_taxonomic_data/taxa_checklist.csv")
        weather,weather_meta=read_github_semicolon(COV1)
        flora,flora_meta=read_github_semicolon(COV2)
        report["source_published_covariates"]={"impervious":weather_meta,"local_flowers":flora_meta}
        urban={norm(r["Id"]):number(r["Urban_500"]) for r in weather}
        flowers={norm(r["Id"]):number(r["SR_all_insect_pollinated_May_August"]) for r in flora}
        if set(coords)!=set(urban) or set(coords)!=set(flowers):
            raise ValueError("Garden ID mismatch across remote source files")
        if any(urban[k] is None or flowers[k] is None for k in coords):
            raise ValueError("Missing published urban/floral covariate")
        colnames=[k for k in visits[0] if k not in NON_TAXON]
        mapping,missing,_=build_guild_map(taxa,colnames)
        if missing:raise ValueError("Pollinator taxonomy incomplete")
        bysite=defaultdict(lambda:defaultdict(float))
        windows=defaultdict(int)
        for r in visits:
            site=norm(r["Id"]);plant=norm(r["phytometer_plant"])
            key=(site,plant)
            windows[key]+=1
            for name in colnames:
                v=number(r.get(name))
                if v is not None and v>0:
                    bysite[key][mapping[name]]+=v
        imagery,_=get_ndvi(coords)
        report["n_sites_good_ndvi"]=sum(bool(v) for v in imagery.values())
        dat={}
        for plant in sorted({p for s,p in windows}):
            items=[]
            for sid in coords:
                t=windows[(sid,plant)]
                if t==0 or not imagery.get(sid):continue
                iv=np.asarray(list(imagery[sid].values()),float)
                items.append({
                    "ndvi":float(iv[:,0].mean()),"evi":float(iv[:,1].mean()),
                    "urban500":urban[sid],"floral_species":flowers[sid],
                    "hoverflies":bysite[(sid,plant)].get("hoverflies",0)/t,
                    "bees":bysite[(sid,plant)].get("bees",0)/t})
            if len(items)!=24:raise ValueError("Unmatched garden predictor "+plant)
            x=np.array([r["ndvi"] for r in items])
            e=np.array([r["evi"] for r in items])
            u=np.array([r["urban500"] for r in items])
            f=np.array([r["floral_species"] for r in items])
            y=np.array([r["hoverflies"] for r in items])
            b=np.array([r["bees"] for r in items])
            # Response transforms reduce leverage of rare high-count windows.
            ytrans=np.log1p(y)
            baseline=np.column_stack((u,f))
            tests={
              "NDVI_vs_urban500":spearman(x,u),
              "NDVI_vs_flower_species_richness":spearman(x,f),
              "NDVI_vs_hoverflies_unadjusted":spearman(x,y),
              "NDVI_vs_bees_unadjusted":spearman(x,b),
              "NDVI_hoverflies_partial_rank_urban_flora":round(partial_rank(x,y,baseline),5)
                  if len(set(y))>1 else None,
              "EVI_hoverflies_partial_rank_urban_flora":round(partial_rank(e,y,baseline),5)
                  if len(set(y))>1 else None,
              "urban500_vs_hoverflies":spearman(u,y),
              "flower_richness_vs_hoverflies":spearman(f,y),
              "model_comparison":{}
            }
            for model_name,modelx in [
                ("urban500_floral_richness",baseline),
                ("urban500_floral_richness_plus_NDVI",np.column_stack((baseline,x))),
                ("urban500_floral_richness_plus_EVI",np.column_stack((baseline,e)))]:
                tests["model_comparison"][model_name]=cv_models(modelx,ytrans)
            base_mae=tests["model_comparison"]["urban500_floral_richness"]["leave_one_garden_out_MAE"]
            tests["gain_NDVI_LOO_MAE_pct"]=round(
                 (base_mae-tests["model_comparison"]["urban500_floral_richness_plus_NDVI"]["leave_one_garden_out_MAE"])
                 /base_mae*100,4)
            tests["gain_EVI_LOO_MAE_pct"]=round(
                 (base_mae-tests["model_comparison"]["urban500_floral_richness_plus_EVI"]["leave_one_garden_out_MAE"])
                 /base_mae*100,4)
            tests["garden_count"]=len(items)
            dat[plant]=tests
        report["NDVI_urban500_full_site_correlation"]=dat["carrot"]["NDVI_vs_urban500"]
        report["NDVI_flower_richness_full_site_correlation"]=dat["carrot"]["NDVI_vs_flower_species_richness"]
        report["by_plant"]=dat
        report["status"]="success"
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"garden_published_covariates_and_NDVI.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:30000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
