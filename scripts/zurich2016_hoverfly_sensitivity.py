#!/usr/bin/env python3
"""Stress test post-hoc carrot-hoverfly NDVI association in 24 Zurich gardens.

Result is an EXPLORATORY association and must not be interpreted as a
causal effect of plant biomass, NDVI on hoverflies, or reproductive fitness.
Output ONLY aggregated statistics; no coordinates or garden IDs.
"""
from collections import Counter,defaultdict
import hashlib,json,math,statistics,time,zipfile,io
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata
from zurich2016_visitor_guild_ndvi import (
    get_source,csvrows,norm,number,build_guild_map,get_ndvi,NON_TAXON)

OUT=Path("analysis_results_guild_spatial_sensitivity")
OUT.mkdir(exist_ok=True)
SEED=20261115
RNG=np.random.default_rng(SEED)

def calc_spearman(x,y):
    ans=spearmanr(x,y)
    return float(ans.statistic),float(ans.pvalue)

def residual_rank_correlation(x,y,covs,perm=5000):
    """Partial rank correlation of x,y after OLS on a priori covariates."""
    rx=rankdata(x);ry=rankdata(y)
    z=np.asarray(covs,dtype=float)
    z=(z-z.mean(axis=0))/(z.std(axis=0)+1e-9)
    C=np.column_stack((np.ones(len(z)),z))
    # Remove deterministic linear geographic trend and optionally total visits.
    ex=rx-C @ np.linalg.lstsq(C,rx,rcond=None)[0]
    ey=ry-C @ np.linalg.lstsq(C,ry,rcond=None)[0]
    rho=float(np.corrcoef(ex,ey)[0,1])
    permvals=[]
    for _ in range(perm):
        shuffled=RNG.permutation(ey)
        permvals.append(float(np.corrcoef(ex,shuffled)[0,1]))
    p=(1+sum(abs(v)>=abs(rho) for v in permvals))/(perm+1)
    return {"partial_rank_corr":round(rho,5),
      "residual_permutation_p_exploratory":round(p,6),
      "covariate_n":z.shape[1], "permutations":perm,
      "warning":"Rank residual OLS controls linear lat/lon trend only; nonspatial permutation does not remove spatial autocorrelation."}

def main():
    report={"source":"2025 Chacko et al EnviDat 10.16904/envidat.676; actual 2016 MOD13Q1 images",
       "hypothesis":"Higher surrounding NDVI predicts more carrot visits by hoverflies but not equivalently by bees",
       "status":"failure","post_hoc_test":True,
       "privacy":"Only aggregated results; original garden coordinates/IDs never exported",
       "source_spatial_units":"24 unique MODIS 250m garden pixels, three 2016 June/July scenes",
       "limitations":["one city and 24 sites; unmeasured local flower availability and urban impervious surface",
         "carrot floral visitation is not true pollen deposition or reproduction fitness",
         "group composition is a proxy, not pollinator functional redundancy",
         "geographic partial Spearman does not model arbitrary spatial correlation",
         "post hoc selected strongest signal after 27 tests, no independent replication"]}
    try:
        raw=get_source()
        report["source_zip_sha256"]=hashlib.sha256(raw).hexdigest()
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            coords={}
            for r in csvrows(z,"03_site_data/garden_site_coordinates.csv"):
                coords[norm(r["Id"])]=(number(r["longitude"]),number(r["latitude"]))
            visit=csvrows(z,"06_trait_data/species_temporal_flower_visitation_matrix.csv")
            check=csvrows(z,"04_taxonomic_data/taxa_checklist.csv")
        if len(coords)!=24 or any(x is None or y is None for x,y in coords.values()):
            raise ValueError("Missing garden sites")
        taxa=[k for k in visit[0] if k not in NON_TAXON]
        mapping,missing,names=build_guild_map(check,taxa)
        if missing:raise ValueError("Missing taxa in pollinator group map")
        garden_visits=defaultdict(Counter)
        windows=Counter()
        for r in visit:
            if norm(r["phytometer_plant"])!="carrot":continue
            sid=norm(r["Id"]);windows[sid]+=1
            for tx in taxa:
                n=number(r.get(tx))
                if n is not None and n>0:
                    garden_visits[sid][mapping[tx]]+=n
        ndvi,_=get_ndvi(coords)
        data=[]
        for sid in coords:
            if sid not in ndvi or not windows[sid]:continue
            xy=coords[sid]
            imag=np.asarray(list(ndvi[sid].values()),float)
            v=garden_visits[sid]
            b=v.get("bees",0);h=v.get("hoverflies",0);w=v.get("wasps",0)
            total=b+h+w
            data.append({"ndvi":float(imag[:,0].mean()),"evi":float(imag[:,1].mean()),
                  "lat":xy[1],"lon":xy[0],
                  "hoverfly_per_record":h/windows[sid],
                  "bee_per_record":b/windows[sid],
                  "total_per_record":total/windows[sid],
                  "hoverfly_bee_logratio":math.log((h+.5)/(b+.5))})
        if len(data)!=24:raise ValueError("Not all 24 gardens matched")
        nd=np.array([r["ndvi"] for r in data])
        hf=np.array([r["hoverfly_per_record"] for r in data])
        bees=np.array([r["bee_per_record"] for r in data])
        ratio=np.array([r["hoverfly_bee_logratio"] for r in data])
        total=np.array([r["total_per_record"] for r in data])
        cov=np.array([[r["lat"],r["lon"]] for r in data])
        report["n_gardens"]=len(data)
        report["baseline_spearman"]={
          "hoverflies":dict(zip(("rho","p_unadjusted"),map(lambda x:round(x,6),calc_spearman(nd,hf)))),
          "bees":dict(zip(("rho","p_unadjusted"),map(lambda x:round(x,6),calc_spearman(nd,bees)))),
          "hoverfly_to_bee_logratio":dict(zip(("rho","p_unadjusted"),map(lambda x:round(x,6),calc_spearman(nd,ratio)))),
        }
        report["geography_adjusted"]={
          "hoverflies":residual_rank_correlation(nd,hf,cov),
          "bees":residual_rank_correlation(nd,bees,cov),
          "hoverfly_bee_logratio":residual_rank_correlation(nd,ratio,cov),
          "hoverfly_controlling_total_visits_and_geography":
              residual_rank_correlation(nd,hf,np.column_stack((cov,np.log1p(total))))
        }
        loo=[]
        for j in range(24):
            keep=np.arange(24)!=j
            loo.append(calc_spearman(nd[keep],hf[keep])[0])
        report["leave_one_garden_out_rho_min_max"]=[round(float(min(loo)),5),round(float(max(loo)),5)]
        # Cluster by garden bootstrap, not by individual visit/observation window.
        boot=[]
        for _ in range(5000):
            indices=RNG.integers(0,24,24)
            if len(set(nd[indices]))<3 or len(set(hf[indices]))<3:continue
            boot.append(calc_spearman(nd[indices],hf[indices])[0])
        report["garden_bootstrap_percentile_95_interval"]=list(map(
            lambda x:round(float(x),5),np.quantile(boot,[0.025,0.975])))
        report["garden_bootstrap_samples"]=len(boot)
        report["NDVI_spatial_range"]=[round(float(nd.min()),5),round(float(nd.max()),5)]
        report["carrot_total_visitor_rate_median"]=round(float(np.median(total)),5)
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
    (OUT/"carrot_hoverfly_NDVI_sensitivity.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:20000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
