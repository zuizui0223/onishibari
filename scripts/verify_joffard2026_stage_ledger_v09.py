#!/usr/bin/env python3
"""Fail CI if v0.9 original-data results differ from committed evidence ledger.

Runs AFTER source script. Checks exact published file input and model output;
no source identifiers exported; does not make causal claims.
"""
import csv,hashlib,json,math
from pathlib import Path
PATH=Path("analysis_results_joffard2026_original_refit_v09/joffard2026_stage_pollen_vs_genetic_compatibility_v09.json")
LEDGER=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_joffard_pollen_carryover_paternity_original_v09.md")
def main():
    x=json.loads(PATH.read_text())
    rows=list(csv.DictReader(LEDGER.open(newline="",encoding="utf8")))
    assert x["status"]=="success_original_pollen_paternity_two_experiment_reanalysis"
    assert len(rows)>=84
    assert len({r["record"] for r in rows})==len(rows)
    assert REPORT.exists()
    assert x["pollen_visitation_sequence"]["n_stigmas"]==186
    assert x["pollen_visitation_sequence"]["n_independent_bumblebees"]==25
    assert x["paternal_diversity_display_experiment"]["n_assigned_offspring"]==1150
    assert x["paternal_diversity_display_experiment"]["n_independent_mothers"]==21
    lookup={r["record"]:r for r in rows}
    checks=[]
    def test(rid,key,value):
        observed=float(lookup[rid][key])
        if not math.isclose(observed,float(value),rel_tol=5e-5,abs_tol=5e-5):
            raise ValueError(f"ledger drift {rid} {key}, source value={value}, ledger={observed}")
        checks.append(rid+":"+key)
    p=x["pollen_visitation_sequence"]
    for rid,src in [
       ("pollen_deposition_sequence_bee_FE_log1p",p["bee_FE_log1p_grain_per_additional_position"]),
       ("pollen_deposition_sequence_remove_first",p["bee_cluster_influence_sensitivity"]),
       ("pollen_deposition_sequence_outlier_trim",p["top_1pct_grains_trimmed_log1p_FE"])
    ]:
        if "remove_first" in rid:
            est=src["drop_first_visit_beta"];ci=src["drop_first_visit_95"]
        elif "outlier_trim" in rid:
            est=src["beta"];ci=src["95"]
        else:
            est=src["beta"];ci=src["cluster95"]
        test(rid,"estimate",est)
        test(rid,"ci_lower",ci[0])
        test(rid,"ci_upper",ci[1])
    test("pollen_deposition_sequence_bee_FE_log1p","p_unadjusted",p["bee_FE_log1p_grain_per_additional_position"]["within_bee_shuffle_p"])
    ppois=p["bee_FE_Poisson_grain_count_sensitivity"]
    test("pollen_deposition_sequence_Poisson_IRR","estimate",ppois["IRR"])
    test("pollen_deposition_sequence_Poisson_IRR","ci_lower",ppois["cluster95_IRR"][0])
    test("pollen_deposition_sequence_Poisson_IRR","ci_upper",ppois["cluster95_IRR"][1])
    for rec,k in [("pollen_deposition_first_visit_mean","1"),("pollen_deposition_tenth_visit_mean","10")]:
        test(rec,"estimate",p["positions"][k]["mean_grains"])
    m=x["paternal_diversity_display_experiment"]
    for rec,key in [
      ("female_display_to_distinct_fathers_per_fruit","mean_sire_count"),
      ("female_display_to_unbiased_simpson_fathers","mean_unbiased_simpson"),
      ("female_display_to_effective_fathers","mean_effective_sire_number"),
      ("female_display_to_12seed_rarefied_fathers","mean_rarefied_sires12"),
      ("female_display_to_within_mother_sire_profile_dissimilarity","within_mother_mean_fruit_sire_dissimilarity")
    ]:
        t=m["treatment_effects_mother_as_unit"][key]
        test(rec,"estimate",t["manipulated_minus_control"])
        test(rec,"ci_lower",t["bootstrap95"][0])
        test(rec,"ci_upper",t["bootstrap95"][1])
        test(rec,"p_unadjusted",t["randomization_permutation_two_sided_p"])
    print(json.dumps({"status":"PASS_JOFFARD_ORIGINAL_TO_LEDGER_CI",
        "matched_metrics":len(checks),
        "ledger_rows":len(rows),
        "original_output_sha256":hashlib.sha256(PATH.read_bytes()).hexdigest(),
        "committed_ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
        "scientific_limit":"within-bee sequence is not within-genet flower number; genetic sire diversity is not total paternal offspring"},indent=2))
if __name__=="__main__":main()
