#!/usr/bin/env python3
"""Do observation-hour and floral-unit-hour visitor rates predict seed/fruit set?

Source Chacko et al. (2025 EnviDat, 2026 J Applied Ecology), 24 Zurich gardens.
This is a LOW N (only 20-23 gardens per phytometer) ecological sensitivity:
the original paper already modeled pollinator groups and reproductive success
using individual plants/umbels, so this does NOT recreate or supplant GLMMs.
It asks whether garden-level count normalization changes conclusions.

For each garden×plant the observer effort (minutes) and per-date flower count
come from author-published daily observations, with source GloPL-like raw data
from EnviDat. Two visit rates:
- visitors / observing hours
- visitors / SUM_D[focal flower units(D) × observing hours(D)]
The latter is a PROXY; flower number was enumerated once per round, not every
moment, and can represent umbels vs individual flowers.
Distinct reproduction endpoints are kept in within-species units, not pooled.
Source geography and private garden ID are never saved in artifacts.
"""
from __future__ import annotations
from collections import defaultdict,Counter
import csv,hashlib,io,json,math,statistics,zipfile
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,rankdata
from zurich2016_flower_date_link import prepare,load,STUDY

OUT=Path("analysis_results_zurich_flower_to_seed")
OUT.mkdir(exist_ok=True)
BASE="https://raw.githubusercontent.com/mrejichacko/BetterBlooms/main/"
CATS={
 "carrot":{
  "visits":"A_Syrphidae_Carrot",
  "file":"07_pollination_success/daucus_carota_seed_set.csv",
  "outcome":"mean_seeds_per_umbel",
  "success":"n_seeds","fail":None,
  "exclude":set(["39"])},
 "radish":{
  "visits":"A_socialBees_Radish",
  "file":"07_pollination_success/raphanus_sativus_fruit_set.csv",
  "outcome":"fruit_set_proportion",
  "success":"n_flowers_with_fruits","fail":"n_flowers_without_fruits",
  "exclude":set(["39"])},
 "sainfoin":{
  "visits":"A_allBees_Sainfoin",
  "file":"07_pollination_success/onobrychis_viciifolia_fruit_set.csv",
  "outcome":"fruit_set_proportion",
  "success":"n_flowers_with_fruits","fail":"n_flowers_without_fruits",
  "exclude":set(["19","28","39","52"])},
 "comfrey":{
  "visits":"A_Bombus_Comfrey",
  "file":"07_pollination_success/symphytum_officinale_seed_set.csv",
  "outcome":"flower_with_seeds_proportion",
  "success":"n_flowers_with_seeds","fail":"n_flowers_without_seeds",
  "exclude":set(["39"])},
}
REPS=2000
GEN=np.random.default_rng(20261115)

def csv_from_bytes(data, delimiter=","):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig","replace")),delimiter=delimiter))

def num(x):
    try:
        z=float(x)
        return z if math.isfinite(z) else None
    except (TypeError,ValueError):return None

def pearson_resid_rank(x,y,z):
    R=np.column_stack([rankdata(x),rankdata(y),*([rankdata(c) for c in z])]).astype(float)
    C=np.column_stack([np.ones(len(R)),R[:,2:]])
    resid1=R[:,0]-C@np.linalg.lstsq(C,R[:,0],rcond=None)[0]
    resid2=R[:,1]-C@np.linalg.lstsq(C,R[:,1],rcond=None)[0]
    return float(np.corrcoef(resid1,resid2)[0,1])

def calc_rho(x,y):
    if len(set(x))<2 or len(set(y))<2:return (None,None)
    z=spearmanr(x,y)
    return float(z.statistic),float(z.pvalue)

def analysis(records,category):
    # all records now include one rate per garden / phytometer and garden outcome
    m=[]
    for key,kind in (("visits_per_hour","per_observed_hour"),
                     ("visits_per_flower_hour_proxy","per_flower_unit_hour")):
        y=np.array([r["reproductive_output"] for r in records],float)
        x=np.array([r[key] for r in records],float)
        covs=[np.array([r[q] for r in records],float) for q in ("urban500","flowering_species_richness")]
        raw,p=calc_rho(x,y)
        partial=pearson_resid_rank(x,y,covs) if raw is not None else None
        bootstrap=[]
        for i in range(REPS):
            ii=GEN.integers(0,len(records),size=len(records))
            r,_=calc_rho(x[ii],y[ii])
            if r is not None and math.isfinite(r):bootstrap.append(r)
        limits=[round(float(t),5) for t in np.quantile(bootstrap,[.025,.975])] if bootstrap else None
        m.append(dict(label=kind,n_gardens=len(records),
            spearman_rho=round(raw,5) if raw is not None else None,
            p_unadjusted=round(p,6) if p is not None else None,
            bootstrap_95pct_percentile_interval=limits,
            urban500_flowering_richness_adjusted_partial_rank_rho=
                round(partial,5) if partial is not None else None,
            unadjusted_p_is_not_source_mixed_model_p=True))
    mx=np.array([r["weighted_flower_count_proxy"] for r in records])
    vy=np.array([r["visits_per_hour"] for r in records])
    d,_=calc_rho(mx,vy)
    return {
       "visit_rate_comparisons":m,
       "same_garden_flower_count_vs_expected_guild_visits_per_hour_rho":round(d,5) if d is not None else None,
       "garden_count":len(records),
       "range_effective_flower_hours":list(map(lambda x:round(float(x),3),
          [min(r["summed_flower_unit_hours"] for r in records),
           max(r["summed_flower_unit_hours"] for r in records)])),
       "no_study_garden_identifiers_exported":True}

def main():
    report={"status":"failed",
      "source":"Reji Chacko et al 2026 JAE DOI 10.1111/1365-2664.70384 and EnviDat 2025 10.16904/envidat.676",
      "question":"Does floral unit-hour normalization of guild landings change correlation with later seed or fruit outcome?",
      "units":["per observation hour","per approximate floral unit × observation hour"],
      "source_plant_guild":{
         p:{"insect_group_column":d["visits"],"reproductive_endpoint":d["outcome"],
            "excluded_garden_count":len(d["exclude"])} for p,d in CATS.items()},
      "limitations":[
        "The original paper already analyzed guild-specific pollinators and seed success using hierarchical plant-level models. This descriptive garden-level reanalysis is NOT a new mechanistic demonstration.",
        "Flower number estimated once per source garden-day (may be umbels rather than individual flowers); not repeated each observer minute.",
        "Normalization by time × flower number changes estimand, not necessarily improves proxy for effective pollination.",
        "Source visitor counts may include non-effective pollen transporters, but chosen guilds track author's selected functional groups.",
        "End-of-season seed/fruit production is subject to soil/watering/resources and genetics as well as visitor frequency.",
        "A garden-level 20-23 observation association has low power and cannot identify mediation.",
        "Site identities and coordinates excluded from exported artifacts.",
        "The outcome seed count and binomial fruit fraction differ; never pool raw values across species."]}
    try:
        field_rows,qa,source=prepare()
        report["field_qa"]=qa
        result_by_garden=defaultdict(lambda:[0.0,0.0])
        for r in field_rows:
            k=(r["site"],r["plant"])
            effort=r["effort_minutes"]/60
            result_by_garden[k][0]+=effort
            result_by_garden[k][1]+=r["flower_units"]*effort
        b=load(STUDY,80_000_000)
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            outcomes={p:csv_from_bytes(z.read(v["file"])) for p,v in CATS.items()}
        rg=load(BASE+"cleaned_data/pollinator_abundance_aggregated_garden_phytometer.txt",500000)
        grads=csv_from_bytes(rg,delimiter=";")
        urn=load(BASE+"raw_data/explanatory_variables.txt",500000)
        urb=csv_from_bytes(urn,delimiter=";")
        frn=load(BASE+"raw_data/plant_floristic_data.txt",500000)
        flor=csv_from_bytes(frn,delimiter=";")
        report["source_sha256"]={"visitor_guild_abundance":hashlib.sha256(rg).hexdigest(),
             "urban":hashlib.sha256(urn).hexdigest(),
             "flowering_species":hashlib.sha256(frn).hexdigest(),
             "reproductive_data_zip":hashlib.sha256(b).hexdigest()}
        ab={str(r["Id"]):r for r in grads}
        urban={str(r["Id"]):num(r["Urban_500"]) for r in urb}
        richness={str(r["Id"]):num(r["SR_all_insect_pollinated_May_August"]) for r in flor}
        report["by_plant"]={}
        all_tests=[]
        for p,spec in CATS.items():
            buckets=defaultdict(lambda:[0.,0.])
            for r in outcomes[p]:
                site=str(r.get("Id","")).strip()
                if site in spec["exclude"]:continue
                if str(r.get("exclude",r.get("excluded",""))).lower() in ("1","yes","true"):
                    continue
                success=num(r.get(spec["success"]))
                if success is None or success<0:continue
                if spec["fail"] is None:
                    n=1.
                else:
                    fail=num(r.get(spec["fail"]))
                    if fail is None or fail<0 or (success+fail)<=0:continue
                    n=success+fail
                buckets[site][0]+=success
                buckets[site][1]+=n
            garden_records=[]
            for site,(numer,denom) in buckets.items():
                if denom<=0 or site not in ab:continue
                eff,Fh=result_by_garden.get((site,p),(0.,0.))
                guild=num(ab[site].get(spec["visits"]))
                if guild is None or guild<0 or eff<=0 or Fh<=0:continue
                u=urban.get(site);f=richness.get(site)
                if u is None or f is None:continue
                garden_records.append({
                    "reproductive_output":numer/denom,
                    "visits_per_hour":guild/eff,
                    "visits_per_flower_hour_proxy":guild/Fh,
                    "weighted_flower_count_proxy":Fh/eff,
                    "summed_flower_unit_hours":Fh,
                    "urban500":u,
                    "flowering_species_richness":f})
            if len(garden_records)<12:raise ValueError("Too few valid gardens "+p+":"+str(len(garden_records)))
            res=analysis(garden_records,p)
            res["germination_or_viability_not_measured"]=True
            res["source_reproduction_rows"]=len(outcomes[p])
            res["source_reproduction_gardens"]=len(buckets)
            res["endpoint_unit"]=spec["outcome"]
            report["by_plant"][p]=res
            for x in res["visit_rate_comparisons"]:
                if x["p_unadjusted"] is not None:
                    all_tests.append((p,x["label"],x["p_unadjusted"]))
            print("PLANT",p,res,flush=True)
        # Benjamini-Hochberg across eight descriptive comparisons; tests
        # are correlated and not independently preregistered research claims.
        ordered=sorted(all_tests,key=lambda x:x[2])
        q=1.
        for i in range(len(ordered)-1,-1,-1):
            p,label,pval=ordered[i]
            q=min(q,pval*len(ordered)/(i+1))
            for x in report["by_plant"][p]["visit_rate_comparisons"]:
                if x["label"]==label:x["q_BH_over_eight_correlated_descriptives"]=round(q,6)
        report["status"]="success"
    except Exception as e:
        report["error"]=type(e).__name__+": "+str(e)
        print("FAIL",report["error"],flush=True)
    (OUT/"visitation_effort_flower_exposure_seed_outcome.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
