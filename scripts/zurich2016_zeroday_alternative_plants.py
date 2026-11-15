#!/usr/bin/env python3
"""Source-checked diagnostic of *why* carrot had excess zero-visit dates.

Query directly observed other flowering phytometers in the same garden/day.
This is an opportunity proxy, NOT a controlled test of pollinator scarcity;
different flowers attract different pollinator guilds, and zero-visit dates
may reflect source observation effort. No site/date identifiers exported.
"""
from collections import defaultdict
import json,math,statistics
from pathlib import Path
from zurich2016_flower_date_link import prepare,PLANTS

OUT=Path("analysis_results_zeroday_diagnostics")
OUT.mkdir(exist_ok=True)

def median_safe(x):
    return round(float(statistics.median(x)),5) if x else None

def main():
    result={"status":"failed",
        "unit":"garden×actual survey date×flowering phytometer with observer minutes positive",
        "source":"EnviDat Reji Chacko et al. 2025, 2026 JAE site/day phytometer and visitor data",
        "privacy":"No garden IDs or dates of particular gardens emitted",
        "limits":["Other plant visit count positive shows insects are active at SAME garden/date, not necessarily compatible carrot pollinators.",
                  "Other plants may bloom at different abundance and attract different functional guilds.",
                  "We do not know individual flower pollen deposition, visit efficiency or male/female reproductive function.",
                  "One year of sampling cannot estimate between-year reliability."]}
    try:
        rows,qa,src=prepare()
        groups=defaultdict(dict)
        for r in rows:
            groups[(r["site"],r["day"])][r["plant"]]=r
        plant_stats={}
        for plant in PLANTS:
            zeros=[];positives=[]
            for r in rows:
                if r["plant"]!=plant:continue
                daygroup=groups[(r["site"],r["day"])]
                other=[z for p,z in daygroup.items() if p!=plant]
                item={
                    "floral_unit_count":r["flower_units"],
                    "effort_minutes":r["effort_minutes"],
                    "n_other_plants_positive_flower_and_effort":len(other),
                    "other_plants_visits_sum":sum(x["all_visitor_count"] for x in other),
                    "any_other_plant_has_visitors":any(x["all_visitor_count"]>0 for x in other),
                    "no_other_plant_observation":not other,
                    "other_radish_visits":(
                        daygroup["radish"]["all_visitor_count"] if "radish" in daygroup and plant!="radish" else None),
                    "visits":r["all_visitor_count"]}
                (zeros if r["all_visitor_count"]==0 else positives).append(item)
            plant_stats[plant]={
                "n_active_flowering_garden_dates":len(zeros)+len(positives),
                "n_zero_visitor_garden_dates":len(zeros),
                "zero_day_flower_count_median":median_safe([x["floral_unit_count"] for x in zeros]),
                "nonzero_day_flower_count_median":median_safe([x["floral_unit_count"] for x in positives]),
                "zero_day_observation_effort_min_median":median_safe([x["effort_minutes"] for x in zeros]),
                "nonzero_day_observation_effort_min_median":median_safe([x["effort_minutes"] for x in positives]),
                "zero_days_with_other_flowering_plant_visit":sum(x["any_other_plant_has_visitors"] for x in zeros),
                "zero_days_with_some_other_flowering_plant_observed":sum(not x["no_other_plant_observation"] for x in zeros),
                "median_other_plant_visits_on_target_zero_days":median_safe([x["other_plants_visits_sum"] for x in zeros]),
                "median_other_plant_visits_on_target_positive_days":median_safe([x["other_plants_visits_sum"] for x in positives]),
                "zero_days_radish_observed":sum(x["other_radish_visits"] is not None for x in zeros),
                "zero_days_radish_positive_visits":sum(x["other_radish_visits"] is not None and x["other_radish_visits"]>0 for x in zeros),
            }
        result["by_plant"]=plant_stats
        result["status"]="success"
    except Exception as ex:
        result["error"]=type(ex).__name__+": "+str(ex)
    (OUT/"same_day_cross_phytometer_zero_visit_diagnostic.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,indent=2,ensure_ascii=False)[:14000])
    if result["status"]!="success":raise RuntimeError(result.get("error"))
if __name__=="__main__":main()
