#!/usr/bin/env python3
"""Validate published-original Daphne donor-location conditioned null v0.37.

Checks 4 exposure definitions, 7 permutations per exposure, original H/H+HF
father sets and tail direction; never treat conditional shuffling as a true
manipulated floral selection experiment or a new focal-species observation.
"""
import json,math
from pathlib import Path
SRC=Path("analysis_results_daphne_spatial_donor_permutation_v37/daphne_2018_2020_H_flower_father_spatial_label_permutation_v37.json")
REPORT=Path("results/2026-11-17_daphne_male_father_spatial_label_falsification_v37.md")
def close(a,b,tol=0.002):return math.isclose(float(a),float(b),abs_tol=tol,rel_tol=tol)
def main():
    j=json.loads(SRC.read_text())
    assert j["status"]=="source_verified_v37_conditional_donor_spatial_label_permutations"
    assert j["source_doi"]=="10.5281/zenodo.8253921"
    rows=j["years_and_sex_definitions"]
    assert len(rows)==4
    expect={
       (2018,"H"):(105,180,39,.625, .286347),
       (2018,"H,HF"):(108,183,40,.6278,.289146),
       (2020,"H"):(71,206,33,.8956,.447479),
       (2020,"H,HF"):(79,214,35,.8979,.402947),
    }
    for x in rows:
        k=(x["year"],",".join(x["donor_sex_codes"]))
        assert k in expect
        n,offspring,mothers,beta,gain=expect[k]
        assert x["n_candidate_genets"]==n
        assert x["n_ssr_assigned_children_in_conditioned_father_risk_set"]==offspring
        assert x["n_recipient_mothers"]==mothers
        assert x["lambda_metres"]==10
        assert close(x["beta_observed"],beta)
        assert close(x["observed_ll_gain_per_child"],gain)
        p=x["permutation_results"]
        assert len(p)==7
        assert p[0]["permutation_scenario"]=="global_shuffle_no_spatial_constraint"
        for v in p:
            assert v["n_permutations"]==1000
            assert v["permutations_changing_a_flower_label"]>=990
            assert v["n_candidate_H_in_multigenet_cells"]>=n*.8
            assert v["tail_fraction_beta_at_least_observed_conditional"]<=.003
            assert v["tail_fraction_gain_at_least_observed_conditional"]<=.003
            assert v["null_beta_ci95"][1]<x["beta_observed"]
            assert v["not_exact_randomization_test"]
        assert len([v for v in p if "within_grid" in v["permutation_scenario"]])==6
    for b in ("spatial_candidate_label_permutations_not_true_random_assignment",
              "one_original_research_population_not_independent_replicates",
              "no_evidence_for_focal_Daphne_pseudomezereum"):
        assert j["limits"][b]
    assert REPORT.exists()
    txt=REPORT.read_text()
    assert all(x in txt for x in ("0.6251","0.8960",".0010","元の花数","実際の花数をランダム割付した実験ではない"))
    print(json.dumps({"status":"PASS_SOURCE_DAPHNE_FATHER_SPATIAL_LABEL_V37",
       "source_years":[2018,2020],"sex_definitions":["H","H_plus_HF"],
       "permutations_per_scenario":1000,"n_spatial_permutation_scenarios_per_year":6,
       "no_unsupported_focal_or_causal_claim":True}),flush=True)
if __name__=="__main__":main()
