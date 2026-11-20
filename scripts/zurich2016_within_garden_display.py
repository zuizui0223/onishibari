#!/usr/bin/env python3
"""Within-garden same-day flower display to insect landing count elasticity.

Plant × garden × actual 2016 sampling date: n=285 active observations.
Prior public PoMS sublinearity came from different quadrats; this audit
instead separates within-garden-plant flowering change from between-site
display differences, adjusting offset log(volunteer observation minutes).

Neither visitor landings nor flower counts measure pollen deposition,
sexual function, seed set or direct plant biomass. Authors' work on this
same dataset (Chacko et al. 2026, JAE DOI 10.1111/1365-2664.70384)
already established guild and reproductive outcomes, so do not claim new
field discovery or infer causation from seasonal flower covariation.
"""
from pathlib import Path
import json,math
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
import patsy
import scipy.stats as ss
import statsmodels.api as sm
from zurich2016_flower_date_link import prepare,PLANTS

OUT=Path("analysis_results_zurich_within_display")
OUT.mkdir(exist_ok=True)
FLOWERS=tuple(PLANTS)
BASE_TERMS=(
 "C(plant) + C(site) + within_carrot + within_radish"
 " + within_sainfoin + within_comfrey"
 " + between_carrot + between_radish + between_sainfoin + between_comfrey")
FORMULAS={
 "no_day_adjustment":BASE_TERMS,
 "actual_survey_calendar_date":BASE_TERMS+" + C(day_str)",
 "plant_specific_smooth_season":BASE_TERMS+" + C(plant):bs(doy,df=3)",
}
def row_table(records):
    dat=pd.DataFrame(records).copy()
    dat["day_str"]=dat["day"].map(lambda x:x.isoformat())
    dat["doy"]=dat["day_of_year"]
    dat["site"]=dat["site"].astype(str)
    dat["plant"]=dat["plant"].astype(str)
    dat["visits"]=dat["all_visitor_count"].astype(float)
    dat["logflower"]=np.log(dat["flower_units"].to_numpy(float))
    dat["offset"]=np.log(dat["effort_minutes"].to_numpy(float))
    dat["siteplant"]=dat["site"]+"_"+dat["plant"]
    dat["meansiteplant"]=dat.groupby("siteplant")["logflower"].transform("mean")
    dat["within"]=dat["logflower"]-dat["meansiteplant"]
    for plant in FLOWERS:
        dat["within_"+plant]=np.where(dat["plant"]==plant,dat["within"],0.0)
        dat["between_"+plant]=np.where(dat["plant"]==plant,dat["meansiteplant"],0.0)
    return dat

def estimate(d,formula,label,family):
    X=patsy.dmatrix(formula,d,return_type="dataframe",NA_action="raise")
    m=sm.GLM(d["visits"].to_numpy(),X,
         family=sm.families.Poisson() if family=="poisson" else sm.families.NegativeBinomial(alpha=1.0),
         offset=d["offset"].to_numpy())
    fit=m.fit(maxiter=90,tol=1e-8,cov_type="cluster",
              cov_kwds={"groups":d["site"].astype(str).to_numpy()})
    result={"model":label,"count_family":family,
       "n_site_plant_dates":len(d),
       "n_gardens":int(d["site"].nunique()),
       "n_site_plant_groups":int(d["siteplant"].nunique()),
       "design_matrix_columns":X.shape[1],
       "pearson_overdispersion":round(float(fit.pearson_chi2/fit.df_resid),4),
       "outcome":"visitors per calendar-day flower display with sampling-effort minutes offset",
       "within_garden_plant_elasticities":{},
       "between_siteplant_elasticities":{}}
    for p in FLOWERS:
        for mode, dest in [("within","within_garden_plant_elasticities"),
                           ("between","between_siteplant_elasticities")]:
            name=mode+"_"+p
            beta=float(fit.params[name])
            se=float(fit.bse[name])
            if not (math.isfinite(beta) and math.isfinite(se) and se>0):
                raise ValueError("Coefficient invalid "+name)
            result[dest][p]={
                "beta":round(beta,6),
                "se_cluster_garden":round(se,6),
                "ci95":[round(beta-1.96*se,6),round(beta+1.96*se,6)],
                "p_vs_beta_zero":round(float(2*ss.norm.sf(abs(beta/se))),7),
                "p_vs_beta_one":round(float(2*ss.norm.sf(abs((beta-1)/se))),7),
                "doubling_F_total_visitors_multiplier":round(2**beta,6),
                "doubling_F_visits_per_flower_unit_multiplier":round(2**(beta-1),6)}
    return result

def main():
    report={"status":"failed",
      "analysis":"2016 Zurich same garden × date × plant flowering display vs visit counts",
      "source_study":"2016 field study published Reji Chacko et al. JAE 2026 doi:10.1111/1365-2664.70384",
      "response":"per-date visitor landings, offset volunteer effort minutes",
      "critical_distinction":"within site/plant seasonal flower change vs differences between garden-plant means",
      "limitations":[
        "Flower unit = carrot umbel versus individual/flowering unit in other phytometer species",
        "All entries are observational; season, plant maturation and resource competition may affect both flowers and insect visits",
        "Per-date flower count sampled once during survey round; exact minute-level flowering exposure unknown",
        "Poisson overdispersion and 24 garden-cluster asymptotics may invalidate precise nominal confidence intervals",
        "Per-floral unit visitor count is not compatible pollen deposition or per-genet reproductive output",
        "Visitor matrix includes some source-recorded nectar robbers/ineffective visitors",
        "Original authors already measured fruit and seed set and guild effects: not a novel biological mechanism",
        "No date-specific actual plant dry biomass, no long-term reproductive strategy evolution"]}
    try:
        rows,qc,provenance=prepare()
        report["source"]=provenance
        report["source_join_qa"]=qc
        dat=row_table(rows)
        if len(dat)!=285:raise ValueError("unexpected source join size "+str(len(dat)))
        report["n_plant_site_date"]=len(dat)
        report["siteplants_with_within_flower_variation"]={}
        for p in FLOWERS:
            d=dat[dat["plant"]==p]
            s=d.groupby("siteplant")["logflower"].agg(["size","nunique"])
            report["siteplants_with_within_flower_variation"][p]={
                "total_siteplant_groups":len(s),
                "groups_with_multiple_days":int(sum(s["size"]>=2)),
                "groups_with_multiple_flower_values":int(sum(s["nunique"]>=2))}
        fits=[]
        for label,formula in FORMULAS.items():
            a=estimate(dat,formula,label,"poisson")
            fits.append(a)
            print("ESTIMATE",label,"within",
                {p:(a["within_garden_plant_elasticities"][p]["beta"],
                    a["within_garden_plant_elasticities"][p]["ci95"])
                for p in FLOWERS},flush=True)
        # Negative binomial as a guard against extreme overdispersion.
        fits.append(estimate(dat,FORMULAS["actual_survey_calendar_date"],
                             "actual_survey_calendar_date","negbin_alpha1"))
        report["estimates"]=fits
        report["status"]="success"
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        print("FAILED",report["error"],flush=True)
    (OUT/"date_matched_within_garden_display_scaling.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
