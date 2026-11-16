#!/usr/bin/env python3
"""Source-derived consistency checks for original 2015 Daphne flower opportunity v23.

Guard: F/H definitions, 15 undefined HF excluded, 166 independent plants,
biological distinction between normalized opportunity and absolute display,
and no interpretation as costs, pollen fertilization or fitness.
"""
import json,math
from pathlib import Path
DATA=Path("analysis_results_daphne_2015_flower_days_v23/daphne_2015_166_FH_seasonal_open_days_v23.json")
REPORT=Path("results/2026-11-16_daphne_sex_number_vs_duration_2015_2018_v23.md")
def near(a,b):return math.isclose(float(a),float(b),abs_tol=.00005,rel_tol=.00005)
def main():
    x=json.loads(DATA.read_text())
    assert x["status"]=="PASS_original_2015_Daphne_window_exposure_v23"
    assert x["n_source_rows"]==905 and x["n_individual_ramets"]==181
    assert x["census_35day_offsets"]==[0,9,18,27,35]
    assert x["sex_counts"]=={"F":76,"H":90,"HF_undefined_excluded_from_FH":15}
    assert x["zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    assert x["csv_sha256"]=="1c3daed5ad7b61af2175ae02194c948e91fb0e3cbfc0ed0ab1c49b5d5e608968"
    q=x["sex_group_summary"]
    assert near(q["F"]["mean_total_flower_per_plant"],23.6842105)
    assert near(q["H"]["mean_total_flower_per_plant"],54.2222222)
    assert near(q["F"]["mean_normalized_open_days"],21.8742129)
    assert near(q["H"]["mean_normalized_open_days"],12.2014751)
    assert near(q["F"]["mean_absolute_flower_days"],500.2302632)
    assert near(q["H"]["mean_absolute_flower_days"],668.5444444)
    assert q["F"]["n_with_last_census_positive"]==61
    assert q["H"]["n_with_last_census_positive"]==16
    c=x["F_vs_H"]
    assert near(c["normalized_open_days"]["F_minus_H_mean"],9.6727378)
    assert near(c["normalized_open_days"]["F_div_H_mean_ratio"],1.7927515)
    assert c["normalized_open_days"]["diff_ci95"][0]>0
    assert near(c["total_flowers"]["F_div_H_mean_ratio"],.436799)
    assert c["absolute_open_flower_days"]["diff_ci95"][0]<0<c["absolute_open_flower_days"]["diff_ci95"][1]
    assert x["unmeasured_true_flower_lifespan"]
    assert x["abs_flower_days_are_display_not_paternal_fitness"]
    assert x["no_Daphne_pseudomezereum_functional_sex_inferred"]
    assert REPORT.exists() and "21.8742" in REPORT.read_text() and "12.2015" in REPORT.read_text()
    print(json.dumps({"status":"PASS_SOURCE_2015_DAPHNE_SEX_OPPORTUNITY_V23",
      "n_independent_F":76,"n_independent_H":90,"n_HF_unclassified":15,
      "normalized_opportunity_F_gt_H":True,
      "absolute_F_minus_H_CI_crosses_zero":True,
      "not_male_or_female_genetic_fitness":True
    }),flush=True)
if __name__=="__main__":main()
