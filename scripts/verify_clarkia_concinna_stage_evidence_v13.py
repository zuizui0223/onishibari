#!/usr/bin/env python3
"""Validate 16 v0.13 source-dependent Clarkia concinna evidence rows.

Runs AFTER the five-study source original script AND the independent
specification stress test in SAME GitHub Actions job. Does not pretend rows
from one source are independent replications and does not transform pollen
surrogate into genetic mating success.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
S1=Path("analysis_results_clarkia_concinna_stages_v13/clarkia_concinna_original_five_stage_gates_v13.json")
S2=Path("analysis_results_clarkia_concinna_robustness_v13/clarkia_study5_pollen_supply_patch_interaction_sensitivity.json")
LED=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_clarkia_flower_signal_pollen_origin_patch_v13.md")

def main():
    a=json.loads(S1.read_text());b=json.loads(S2.read_text())
    rows=list(csv.DictReader(LED.open(encoding="utf-8",newline="")))
    d={r["record"]:r for r in rows}
    assert a["status"]=="success_original_five_study_stage_evidence_audit"
    assert b["status"]=="success_secondary_specification_robustness_audit"
    assert a["original_source_doi"]==b["original_source_DOI"]=="10.5061/dryad.8pk0p2nz2"
    assert a["original_source_counts"]=={"1":1626,"2":137,"3":130,"4":121,"5":204}
    assert a["self_vs_neighbor_pigment"]["source_manipulated_plant_stigmas"]==136
    assert a["focal_maternal_seed"]["after_complete_focal_outcome_n"]==105
    assert len(rows)>=141 and len(rows)==len(d)
    assert sum(z.startswith("v13_clarkia_") for z in d)==16
    assert REPORT.exists() and "omnibus" in REPORT.read_text().lower()
    ncheck=0
    def check(record,n,est,ci=None,p=None):
        nonlocal ncheck
        if record not in d:raise AssertionError("missing record "+record)
        v=d[record]
        if int(v["n"])!=int(n):raise AssertionError("incorrect independent unit count "+record)
        if v["measured_floral_biomass"]!="False":
            raise AssertionError("wrong true floral biomass label "+record)
        if v["primary_doi"]!="10.3120/0024-9637-72.1.a1":
            raise AssertionError("source DOI mismatch "+record)
        for field,target in [("estimate",est),("ci_lower",ci[0] if ci else None),
              ("ci_upper",ci[1] if ci else None),("p_unadjusted",p)]:
            if target is None:continue
            value=float(v[field])
            if not math.isclose(value,float(target),rel_tol=1e-5,abs_tol=1e-5):
                raise AssertionError(f"{record}/{field} source {target} evidence {value}")
            ncheck+=1
    q=a["consecutive_probe_movement"]
    c=q["overall_NB_log_display_to_additional_probes"]
    check("v13_clarkia_withinplant_additional_probes_per_log_display",
          q["unique_insect_tag_groups_genus_number"],c["beta"],c["CI_95"])
    for key,rec in (
       ("self_pigment_arrival_vs_male_display","v13_clarkia_self_pigment_prob_log_male_display"),
       ("neighbor_pigment_arrival_vs_male_display","v13_clarkia_neighbor_pigment_prob_log_male_display")):
        c=a["self_vs_neighbor_pigment"][key]["ln_male"]
        check(rec,136,c["beta"],c["CI_95"],c["p"])
    c=a["focal_maternal_seed"]["plant_covariate_adjusted_male_display_to_focal_fruit_seedset_beta"]
    check("v13_clarkia_male_display_to_focal_seed_fraction",105,c["beta"],c["CI_95"],c["p"])
    c=a["patch_20minute_visitation"]["ln_patch_to_total_20min_visits_NB"]["ln_patch"]
    check("v13_clarkia_patchsize_to_20min_flower_visits",121,c["beta"],c["CI_95"],c["p"])
    original=a["display_vs_pollen_presence_patch_factorial"]
    for record,key in (
      ("v13_clarkia_signal_only_large_patch_6e_vs_1_grains","signal_only_emasculated6_vs_single_large_patch"),
      ("v13_clarkia_signal_only_small_patch_6e_vs_1_grains","signal_only_emasculated6_vs_single_small_patch"),
      ("v13_clarkia_pollenflower_vs_emasculated_large_patch_grains","pollen_supplying6_vs_emasculated6_large_patch"),
      ("v13_clarkia_pollenflower_vs_emasculated_small_patch_grains","pollen_supplying6_vs_emasculated6_small_patch"),
      ("v13_clarkia_patch_intact_vs_emasculated_did","pollen_supply_patch_interaction")):
        c=original["direct_three_arm_contrasts_on_log1p_pollen_grains"][key]
        check(record,204,c["estimate_log1p"],c["CI_95"],c["p"])
    for record,key in (
      ("v13_clarkia_patch_intact_vs_emasculated_did_site_adjusted","original_scale_log1p_date_and_location_covariate"),
      ("v13_clarkia_patch_intact_vs_emasculated_did_NB","raw_counts_NB_alpha_1_date_covariate")):
        c=b["models"][key]
        check(record,204,c["estimate_log_scale"],c["HC3_95"],c["p_two_sided_exploratory"])
    c=b["date_cluster_bootstrap_on_unadjusted_one_df_contrast"]
    check("v13_clarkia_patch_intact_vs_emasculated_did_date_bootstrap",13,
          c["observed_unadjusted_difference_in_differences_log1p"],c["date_resampling95"])
    c=b["within_day_patch_treatment_permutation"]
    check("v13_clarkia_patch_intact_vs_emasculated_did_permutation",204,
          c["observed_DiD"],p=c["two_sided_stratified_label_p"])
    c=original["ln_time_to_first_pollen_simple_OLS_not_survival_model"]["C(patch)[T.small]"]
    check("v13_clarkia_small_patch_to_first_pollen_delay",203,c["beta"],c["CI_95"],c["p"])
    c=original["fraction_of_flower_life_remaining_after_first_pollen_simple_OLS"]["C(patch)[T.small]"]
    check("v13_clarkia_small_patch_to_remaining_flower_lifespan",203,c["beta"],c["CI_95"],c["p"])
    assert ncheck>=48,ncheck
    assert b["source_randomized_display_but_patch_partly_natural_and_thinned"] is True
    assert original["pollen_grains_are_not_compatible_pollen_tubes"] is True
    assert b["reported_original_omnibus_interaction_p_approx"]>=.05
    print(json.dumps({"status":"PASS_Clarkia_concinna_v13_source_and_ledger",
       "stage_evidence_rows":len(rows),
       "new_nonindependent_v13_evidence_rows":16,
       "numerical_comparisons":ncheck,
       "source_study5_stigmas":204,
       "source_study2_stigmas":136,
       "source_study3_complete_fruits":105,
       "original_omnibus_p_reported_approx":b["reported_original_omnibus_interaction_p_approx"],
       "posthoc_did_not_claim_confirmatory":True,
       "source_result1_SHA256":hashlib.sha256(S1.read_bytes()).hexdigest(),
       "source_result2_SHA256":hashlib.sha256(S2.read_bytes()).hexdigest(),
       "ledger_SHA256":hashlib.sha256(LED.read_bytes()).hexdigest()},indent=2))
if __name__=="__main__":main()
