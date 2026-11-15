#!/usr/bin/env python3
"""2016 Zurich flower visitation temporal concentration against Poisson null.

Condition on each GARDEN×PLANT total 2016 insect landings. Under a null where
each garden/plant's visits occur at a constant rate PER FLOWER-UNIT ×
OBSERVER-MINUTE across its sampled days, allocation to days is multinomial
with p_day ∝ floral_units(day)*observer_minutes(day).

This is a *conditional effort and flower-exposure null*, not a mechanistic
process model. Clustering beyond null may reflect weather, developmental
phenology, observer effort quality, flower nectar rewards, species flight
phenology or independent insect behaviors; cannot attribute it to pollinator
collapse or evolutionary selection without extra controls.

Important: original time-series visit matrix includes some illegitimate
insects; visitor events are NOT necessarily compatible pollen deposition.
No original garden IDs or precise locations released in output.
"""
from collections import defaultdict
import json
import math
import numpy as np
from pathlib import Path
from zurich2016_flower_date_link import prepare,PLANTS

OUT=Path("analysis_results_zurich_temporal_reliability")
OUT.mkdir(exist_ok=True)
N_SIM=10000

def evaluate(groups,seed=20261115,floral_exposure=True):
    RNG=np.random.default_rng(seed)
    case=[]
    for (site,plant),days in sorted(groups.items()):
        days=sorted(days,key=lambda x:x["day"])
        if len(days)<2:continue
        n=int(round(sum(r["all_visitor_count"] for r in days)))
        if n<=0:continue  # all-zero gardens can't inform temporal allocation
        counts=np.array([int(round(r["all_visitor_count"])) for r in days])
        minutes=np.array([r["effort_minutes"] for r in days],float)
        flowers=np.array([r["flower_units"] for r in days],float)
        if not np.all(minutes>0) or not np.all(flowers>0):
            raise ValueError("Nonpositive exposure in observed day")
        exposure=minutes*flowers if floral_exposure else minutes
        p=exposure/exposure.sum()
        sim=RNG.multinomial(n,p,size=N_SIM)
        obs_zero=int(np.sum(counts==0))
        sim_zero=np.sum(sim==0,axis=1)
        obs_hhi=float(np.sum((counts/n)**2))
        sim_hhi=np.sum((sim/n)**2,axis=1)
        case.append({"plant":plant,"n_days":len(days),
             "n_visits":n,"observed_zero":obs_zero,
             "sim_zero":sim_zero,"observed_HHI":obs_hhi,"sim_HHI":sim_hhi})
    output={}
    for plant in PLANTS:
        active=[x for x in case if x["plant"]==plant]
        if not active:continue
        obs_zero=sum(x["observed_zero"] for x in active)
        sim_zero=sum((x["sim_zero"] for x in active))
        obs_hhi=float(np.mean([x["observed_HHI"] for x in active]))
        sim_hhi=np.mean(np.stack([x["sim_HHI"] for x in active]),axis=0)
        p_excess=(1+np.count_nonzero(sim_zero>=obs_zero))/(1+N_SIM)
        p_deficit=(1+np.count_nonzero(sim_zero<=obs_zero))/(1+N_SIM)
        p_conc=(1+np.count_nonzero(sim_hhi>=obs_hhi))/(1+N_SIM)
        output[plant]={
            "siteplant_groups_with_2plus_days_and_nonzero_total":len(active),
            "total_observation_site_days":sum(x["n_days"] for x in active),
            "visitor_events_in_siteplants":sum(x["n_visits"] for x in active),
            "observed_zero_visit_site_days":obs_zero,
            "null_expected_zero_visit_days":round(float(np.mean(sim_zero)),5),
            "null_zero_visit_days_95pct_interval":list(map(float,np.quantile(sim_zero,[.025,.975]))),
            "one_sided_excess_zero_p_conditional":round(float(p_excess),6),
            "one_sided_deficit_zero_p_conditional":round(float(p_deficit),6),
            "observed_mean_within_siteplant_visit_HHI":round(obs_hhi,6),
            "null_mean_within_siteplant_visit_HHI":round(float(np.mean(sim_hhi)),6),
            "one_sided_excess_HHI_p_conditional":round(float(p_conc),6),
            "n_null_simulations":N_SIM,
            "null":"within each garden/plant total visits fixed; visit dates allocated multinomial with p∝exposure",
        }
    return output

def main():
    report={"status":"failed",
       "name":"Temporal reliability / excess-zero conditional null in flowering site-days",
       "source":"Zurich 2016 original dated flower and visitor data, fieldwork metadata in Chacko et al. 2026",
       "scientific_units":"site x plant x real observation day, count of visitors; no direct effective pollen",
       "source_spatial_privacy":"No actual garden identifiers/coordinates exported",
       "strong_caveats":[
          "Excess visits concentrated in fewer dates does not imply pollinator shortage; weather and plant blooming phenology vary",
          "Observation effort can vary by volunteer and weather even at equal recorded minutes",
          "Null constant per-floral-unit intensity within a garden across time is deliberately simple and may be unrealistic",
          "Null is conditional on observed total events per garden/plant, so it cannot estimate absolute pollinator population size",
          "Insect visits on flowering target may include ineffective pollen delivery or nectar robbing",
          "Sites with zero visitors over all days do not contribute to conditional timing inference",
          "No paternity, flower biomass, survival, seed set or animal population density inferred"]}
    try:
        rows,qa,src=prepare()
        report["source_qa"]=qa
        report["source"]=src
        gr=defaultdict(list)
        for r in rows:
            gr[(r["site"],r["plant"])].append(r)
        report["flower_unit_times_effort_null"]=evaluate(gr,seed=20261115,floral_exposure=True)
        report["effort_only_null_sensitivity"]=evaluate(gr,seed=20261116,floral_exposure=False)
        report["status"]="success"
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        print("FAILED",report["error"],flush=True)
    (OUT/"date_conditional_zero_concentration_null.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2)[:20000])
    if report["status"]!="success":raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
