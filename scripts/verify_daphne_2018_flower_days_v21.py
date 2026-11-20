#!/usr/bin/env python3
"""Verify published-source v21 female versus hermaphrodite flower-days result.
This protects denominators, source identities, missingness handling and no causal claims.
"""
import json,math
from pathlib import Path
FILE=Path("analysis_results_daphne_flower_days_v21/daphne_2018_flower_opportunity_auc_source_v21.json")
REPORT=Path("results/2026-11-16_daphne_flower_days_sex_opportunity_v21.md")
def close(a,b):return math.isclose(float(a),float(b),abs_tol=0.000015,rel_tol=0.000015)
def main():
    r=json.loads(FILE.read_text())
    assert r["status"]=="source_original_2018_flower_days_v21"
    assert r["source_zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    assert r["source_table_sha256"]=="140aeae0064ff9e13ff4199794abe7dc533d5ac4f3fa0e929fc815eb0d4df017"
    assert r["n_source_censuses"]==408 and r["n_inflorescences"]==51 and r["n_unique_ramets"]==14
    assert r["n_missing_open_flower_counts"]==9
    assert r["census_day_offsets"]==[0,2,5,7,12,15,20,24]
    expected={
      "F|bag":(5,0,22.253333,22.253333),
      "F|natural":(5,0,22.733333,22.733333),
      "F|outcross":(5,5,18.54,22.14),
      "H|emasculate_bag":(9,0,15.188404,15.188404),
      "H|emasculate_outcross":(9,4,13.74709,15.358201),
      "H|natural":(9,0,16.084656,16.084656),
      "H|self":(9,0,13.570304,13.570304)
    }
    assert set(r["groups"])==set(expected)
    for k,(n,missing,lo,hi) in expected.items():
        p=r["groups"][k]
        assert p["n_ramets"]==n and p["n_missing"]==missing,(k,p)
        assert close(p["mean_auc_low"],lo) and close(p["mean_auc_high"],hi),(k,p)
        assert p["mean_auc_high"]>=p["mean_auc_low"]
    sex=r["primary_unpollinated_sex_comparison"]
    assert sex["F_n"]==5 and sex["H_n"]==9 and sex["n_missing_primary"]==0
    assert close(sex["F_minus_H_auc_days"],7.064929)
    assert close(sex["F_to_H_auc_ratio"],1.465153)
    assert sex["sex_group_boot_CI95"][0]>0
    assert sex["not_causal_sex_effect"] and sex["H_emasculation_is_not_identical_manipulation_to_F_bagging"]
    labels={x["contrast"]:x for x in r["same_plant_treatment_label_sensitivities"]}
    assert labels["self minus emasculate_bag"]["matched_ramets"]==9
    assert labels["self minus emasculate_bag"]["n_treat_smaller"]==8
    assert close(labels["self minus emasculate_bag"]["interpolated_mean"],-1.6181)
    assert labels["outcross minus bag"]["n_matched_with_any_missing"]==1
    assert r["measurement_boundaries"]["not_exact_flower_individual_lifetime"]
    assert r["measurement_boundaries"]["not_measured_resource_maintenance_mass_or_CO2"]
    assert r["measurement_boundaries"]["not_fitness_or_compatible_pollination"]
    assert r["measurement_boundaries"]["not_focal_Daphne_pseudomezereum"]
    assert REPORT.exists() and "22.2533" in REPORT.read_text() and "15.1884" in REPORT.read_text()
    print(json.dumps({"status":"PASS_ORIGINAL_DAPHNE_2018_SEX_FLOWER_DAYS_V21",
      "original_records":408,"independent_ramets":14,
      "n_flower_days_groups":7,"literal_missing_not_zero":9,
      "focal_untested":True,"no_construction_cost_or_plant_fitness_claim":True
    }),flush=True)
if __name__=="__main__":main()
