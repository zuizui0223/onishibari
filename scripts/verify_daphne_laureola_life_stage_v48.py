#!/usr/bin/env python3
"""Validate 2005 published Daphne laureola source table and no false fitness inference."""
import json,math
from pathlib import Path
P=Path("analysis_results_daphne_life_stage_filters_v48/daphne_laureola_2005_5pop_mature_seed_to_adult_filter_v48.json")
R=Path("results/2026-11-17_daphne_early_fruit_vs_late_genetic_recruitment_v48.md")
def near(x,y):return math.isclose(float(x),float(y),rel_tol=1e-5,abs_tol=1e-5)
def main():
    x=json.loads(P.read_text(encoding="utf8"))
    assert x["status"]=="PASS_2005_DAPHNE_LAUREOLA_SOURCE_MATING_STAGE_AUDIT_V48"
    assert x["doi"]=="10.1038/sj.hdy.6800550"
    assert x["n_geographic_populations"]==5
    assert x["n_H_maternal_seed_families"]==88
    assert x["n_H_genotyped_mature_seeds"]==1678
    assert x["n_F_genotyped_mature_seeds"]==265
    assert x["n_adult_plants_scored_morph_sex_not_genotyped"]==486
    assert x["n_adult_plants_genotyped_for_FIS_years_2003_2004"]==138
    assert x["n_populations_authors_report_adult_FIS_below_mating_equilibrium"]==4
    a={r["site"]:r for r in x["site_original_aggregates_and_derived_equilibrium_scenarios"]}
    assert set(a)=={"CV","RH","NC","CT","CE"}
    t={"CV":(.024,.024,40.3,0.,.399539),
       "RH":(.026,.020,58.3,6.9,.4726456),
       "NC":(.001,.001,84.7,11.1,.5364535),
       "CT":(.125,None,73.6,11.1,.4673439),
       "CE":(.026,None,31.9,1.4,.6201731)}
    for site,(s,m,f,h,expect) in t.items():
        row=a[site]
        assert near(row["H_t_s_single_locus_as_reported"],s)
        assert row["H_t_m_multilocus_as_reported_or_null"]==m
        assert near(row["stigma_without_any_pollen_F_percent"],f)
        assert near(row["stigma_without_any_pollen_H_percent"],h)
        assert near(row["pop_mating_equilibrium_Fexp_scenario"],expect)
        assert row["not_an_adult_survival_probability"]
    assert a["NC"]["source_published_t_estimate_bootstrap_SD"]["H_t_m"]==.142
    assert a["RH"]["F_t_m_as_reported"]==1.168
    assert a["RH"]["F_t_used_for_biologically_clipped_calculation"]==1.0
    assert a["NC"]["F_t_m_as_reported"] is None
    assert a["CE"]["F_t_s_as_reported"] is None
    assert a["CE"]["published_adult_FIS_vs_equilibrium_classification"]=="similar"
    assert x["original_source_t_estimate_parentheses_are_bootstrap_SD_not_independent_seed_SE"]
    limits=x["inference_constraints"]
    for item in (
        "printed_site_mating_estimates_not_raw_seed_genotypes",
        "2_populations_H_t_m_absent_due_to_one_variable_locus",
        "2_populations_F_outcrossing_unmeasured_assumed_1_for_scenario",
        "FIS_points_approximate_figure_do_not_use_as_exact_numeric_data",
        "no_causal_late_inbreeding_depression_directly_tracked_to_adults",
        "no_direct_Onishibari_reproduction_experiment"):
        assert limits[item],item
    assert R.exists()
    rr=R.read_text(encoding="utf8")
    for val in ("1,678","138株","0.142","0.400","0.620","同じ種子を生存するまで追った直接試験ではない"):
        assert val in rr,val
    print(json.dumps({"status":"PASS_SOURCE_DAPHNE_LAUREOLA_FIVE_POP_STAGE_GATE_V48",
      "independent_published_study":1,
      "five_populations":5,"H_genotyped_seeds":1678,
      "adult_genets":138,
      "Fexp_is_conditional_not_individual_survival":True,
      "no_Onishibari_maternal_function_inferred":True},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
