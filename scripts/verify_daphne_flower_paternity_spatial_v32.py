#!/usr/bin/env python3
"""CI audit of numerical and interpretation boundaries of Daphne v0.32.
Checks H male display scaling under original-source, exact-F, top-donor,
and geographical-block holdout scenarios. Does NOT establish causation.
"""
import json,math
from pathlib import Path
FILE=Path("analysis_results_daphne_fatherhood_robust_v32/daphne_H_flower_paternity_adversarial_v32.json")
REPORT=Path("results/2026-11-16_daphne_male_display_paternity_female_residual_v32.md")

def close(x,y):
    return math.isclose(float(x),float(y),rel_tol=1e-4,abs_tol=1e-4)

def main():
    d=json.loads(FILE.read_text())
    assert d["status"]=="v32_original_Daphne_male_flower_count_paternity_source_robustness"
    assert d["original_data"]=="10.5281/zenodo.8253921" and d["lambda_metres_preselected"]==10
    y={v["year"]:v for v in d["years"]}
    assert set(y)=={2018,2020}
    assert y[2018]["source_matched_H_children"]==180
    assert y[2020]["source_matched_H_children"]==206
    assert y[2018]["source_H_candidate_genets"]==105
    assert y[2020]["source_H_candidate_genets"]==71
    assert y[2018]["source_mothers"]==39 and y[2020]["source_mothers"]==33
    assert y[2018]["mothers_exact_F_code"]==35 and y[2020]["mothers_exact_F_code"]==32
    assert close(y[2018]["baseline_beta"],.625023)
    assert close(y[2020]["baseline_beta"],.895633)
    for year in (2018,2020):
        row=y[year]
        opts={z["scenario"]:z for z in row["subset_sensitivity"]}
        assert len(opts)==4
        for item in opts.values():
            assert item["beta"]>0
            assert item["beta_mother_boot_CI95"][0]>0
            assert item["per_child_gain_distance_only"]>0
        assert opts["only_exact_F_mothers"]["remaining_sampled_maternal_families"]==row["mothers_exact_F_code"]
        assert opts["exclude_top_10pct_flower_H_genets"]["remaining_confident_H_paternity_assignments"] in (123,122)
        assert len(row["spatial_block_holdout"])==4
        for b in row["spatial_block_holdout"]:
            assert b["spatial_block_holdout_log_score_gain_per_child"]>0
            assert b["n_blocks"]>=5
            assert b["n_test_assigned_offspring"]==row["source_matched_H_children"]
    assert all(d["limits"][key] for key in (
      "spatial_blocks_do_not_remove_environmental_confounding",
      "not_full_population_paternal_seed_production",
      "not_focal_D_pseudomezereum"))
    report=REPORT.read_text()
    assert all(x in report for x in (
        "0.6250","0.8956","10","2018","2020",
        "3.59","0.624","0.908","10.1002/ajb2.16239"))
    print(json.dumps({
       "status":"PASS_DAPHNE_ORIGINAL_H_FLORAL_SIRING_V32",
       "year_2018_H_assigned":180,
       "year_2020_H_assigned":206,
       "all_control_beta_95pct_intervals_positive":True,
       "all_spatial_block_holdout_scores_positive":True,
       "original_published_source_not_new_D_pseudomezereum_data":True
    }))
if __name__=="__main__":main()
