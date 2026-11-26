#!/usr/bin/env python3
"""Defensive checks for HF paternal origin and local-fruit model sensitivity.
Protect exact source article counts, tiny HF parent sample warnings and no
permutation of source identities or release of row-level geolocation/genotypes.
"""
import json,math
from pathlib import Path
FILE=Path("analysis_results_daphne_HF_true_nearby_v29/daphne_2018_2020_HF_true_paternity_5m_sensitivity_v29.json")
REPORT=Path("results/2026-11-16_daphne_HF_true_paternity_local_fruit_v29.md")
def eq(x,y,tol=1e-5):return math.isclose(float(x),float(y),rel_tol=tol,abs_tol=tol)
def main():
    data=json.loads(FILE.read_text())
    assert data["status"]=="Daphne_2018_2020_HF_true_father_distance_sensitivity_v29"
    assert data["not_a_same_flower_mediation_estimate"]
    assert data["no_individual_ID_or_location_released"]
    assert len(data["by_year"])==2
    years={x["year"]:x for x in data["by_year"]}
    for y,n,hf_genets,HF_child,hf_lt10 in [(2018,50,1,3,1),(2020,56,3,8,7)]:
        v=years[y]
        assert v["n_source_F_matched_plant"]==n
        assert v["n_true_HF_sire_genets"]==hf_genets
        q=v["father_distances_by_code"]["HF_exact"]
        assert q["assigned_seeds"]==HF_child
        assert q["distinct_father_genets"]==hf_genets
        assert q["n_dist_lt10"]==hf_lt10
        p=v["spearman_F_fruit_rate_vs_5m_flowers"]
        assert p["H_plus_all_HF"]>=p["H_only"]
        assert p["confirmed_HF_increment_plant_boot_95%"][0]<0<p["confirmed_HF_increment_plant_boot_95%"][1]
        assert p["all_HF_increment_plant_boot_95%"][0]<0<p["all_HF_increment_plant_boot_95%"][1]
        assert v["limitations"]["observed_sire_filter_is_outcome_conditioned"]
        assert v["limitations"]["no_D_pseudomezereum_data"]
        assert len(v["grid_mean_sensitivity"])==4
        assert max(vv["n_blocks_with_females"] for vv in v["grid_mean_sensitivity"])<=9
    s=years[2020]["spearman_F_fruit_rate_vs_5m_flowers"]
    assert eq(s["H_only"],0.2356452)
    assert eq(s["H_plus_only_HF_with_confirmed_paternity"],0.2604005)
    assert eq(s["H_plus_all_HF"],0.3746982)
    assert eq(s["all_HF_increment_delta_rho"],0.139053)
    text=REPORT.read_text()
    assert all(x in text for x in ("0.2356","0.2604","0.3747","7/8"))
    print(json.dumps({"status":"PASS_DAPHNE_HF_GENETIC_FATHER_NEIGHBORHOOD_V29",
        "years":[2018,2020],"HF_paternal_offspring_total":11,
        "HF_genet_fathers_by_year":[1,3],
        "all_HF_increment_CI_crosses_zero":True,
        "confirmed_HF_increment_CI_crosses_zero":True,
        "individual_source_private":True,
        "D_pseudomezereum_not_observed":True}),flush=True)
if __name__=="__main__":main()
