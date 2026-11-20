#!/usr/bin/env python3
"""Verify original source v0.36 loss-of-support diagnosis and safe interpretation."""
import json,math
from pathlib import Path
FILE=Path("analysis_results_daphne_lagged_fruiting_identification_v36/daphne_2018fruit_2020flower_common_support_v36.json")
REPORT=Path("results/2026-11-17_daphne_residual_female_sex_dynamics_cost_identification_v36.md")
def close(x,y):return math.isclose(float(x),float(y),rel_tol=.00001,abs_tol=.00001)
def main():
    j=json.loads(FILE.read_text())
    assert j["status"]=="source_original_Daphne_cost_identification_common_support_v36"
    groups={q["sex_exact_both_years"]:q for q in j["by_sex"]}
    assert set(groups)=={"H","F"}
    h=groups["H"];f=groups["F"]
    assert h["n_genets"]==45 and h["2018_fruit_positive_n"]==6
    assert h["2018_flower_number_if_fruited"]["min"]==51
    assert h["2018_flower_number_if_fruited"]["max"]==360
    assert h["2018_flower_number_if_not_fruited"]["min"]==1
    assert h["2018_flower_number_if_not_fruited"]["max"]==114
    assert h["common_support_2018_flower_range"]==[51,114]
    assert h["n_genets_common_support"]==11
    assert h["n_fruiting_genets_common_support"]==1
    assert h["n_nonfruiting_genets_common_support"]==10
    assert h["common_support_regression_coef"] is None
    assert h["nearest_baseline_flower_match_diagnostic"]["n_with_nearest_control_differing_by_lt_20pct_in_baseline_flower_number"]==3
    assert close(h["nearest_baseline_flower_match_diagnostic"]["mean_matched_difference_log_flower2020_with_replacement_no_uncertainty"],.445848)
    assert close(h["baseline_adjusted_forward_fruit2018_to_flower2020_coef"],-1.201643)
    assert f["n_genets"]==19 and f["2018_fruit_positive_n"]==9
    assert f["n_fruiting_genets_common_support"]==5
    assert f["n_nonfruiting_genets_common_support"]==8
    assert j["interpretation"]["observational_forward_association_not_direct_resource_cost"]
    assert j["interpretation"]["unobserved_2019_and_survivor_selection_cannot_be_fixed_by_matching"]
    assert j["interpretation"]["not_focal_D_pseudomezereum"]
    assert REPORT.exists()
    txt=REPORT.read_text()
    assert all(t in txt for t in ("51–114","1/6","−1.2016","+0.446","v0.35","v0.36"))
    print(json.dumps({"status":"PASS_DAPHNE_2018_2020_FRUIT_COST_IDENTIFICATION_V36",
      "H_genets":45,"H_2018_fruited":6,"common_support_fruited":1,
      "same_size_causal_cost_not_identified":True,
      "reverse_time_and_unobserved_2019_not_causal":True,
      "no_focal_Onishibari_fitness_claim":True}))
if __name__=="__main__":main()
