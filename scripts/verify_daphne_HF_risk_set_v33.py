#!/usr/bin/env python3
"""Enforce source/sex-code HF sensitivity constraints and typed source numbers v33."""
import json,math
from pathlib import Path
P=Path("analysis_results_daphne_HF_coding_father_model_v33/daphne_2018_2020_HF_inclusion_father_choice_v33.json")
R=Path("results/2026-11-16_daphne_HF_coding_father_choice_v33.md")
def near(a,b):
    return math.isclose(float(a),float(b),abs_tol=1e-4,rel_tol=1e-4)
def main():
    x=json.loads(P.read_text())
    assert x["status"]=="source_verified_Daphne_H_plus_HF_father_choice_v33"
    assert x["source_doi"]=="10.5281/zenodo.8253921"
    assert len(x["yearly"])==2
    o={y["year"]:y for y in x["yearly"]}
    for yr,ng,ns,nd,beta,fruit,flow in [
        (2018,108,183,40,.627779,67,4917),
        (2020,79,214,35,.897867,96,2939)
    ]:
        z=o[yr]
        assert z["paternity_source_data_counts_H_plus_HF"]["eligible_source_H_flowering_genets"]==ng
        assert z["paternity_source_data_counts_H_plus_HF"]["assigned_H_eligible_for_conditional_model"]==ns
        assert z["paternity_source_data_counts_H_plus_HF"]["independent_sampled_mothers_with_eligible_H_offspring"]==nd
        assert z["adult_fruit_rate_sex_code_sensitivity"]["H_plus_HF"]["fruit_sum"]==fruit
        assert z["adult_fruit_rate_sex_code_sensitivity"]["H_plus_HF"]["flower_sum"]==flow
        assert len(z["flower_to_father_choice_model"])==3
        by_lambda={m["lambda_m"]:m for m in z["flower_to_father_choice_model"]}
        assert set(by_lambda)=={5,10,20}
        assert near(by_lambda[10]["H_plus_HF_beta"],beta)
        for m in by_lambda.values():
            assert m["H_plus_HF_beta"]>0
            assert m["H_plus_HF_cluster_CI95"][0]>0
            assert m["H_plus_HF_mother_holdout_gain_per_child"]>0
            assert m["n_assigned_offspring_H_plus_HF"]==ns
        assert z["adult_fruit_rate_sex_code_sensitivity"]["F_exact"]["pooled_fruit_rate"]>z["adult_fruit_rate_sex_code_sensitivity"]["H_plus_HF"]["pooled_fruit_rate"]
    assert o[2020]["paternity_source_data_counts_H_plus_HF"]["assigned_H_excluded_due_to_no_positive_flower_or_coordinate"]==21
    assert all(x["limits"][k] for k in (
        "HF_sex_code_may_be_assigned_using_prior_fruiting_outcome",
        "comparing_H_vs_HplusHF_conditions_on_different_sire_sets",
        "not_Daphne_pseudomezereum_experiment"))
    assert R.exists()
    report=R.read_text()
    assert all(s in report for s in ("0.6278","0.8979","1.36%","3.27%","21"))
    print(json.dumps({"status":"PASS_DAPHNE_HF_CODE_RISK_SET_V33",
      "2018_conditional_sire_offspring":183,"2020_conditional_sire_offspring":214,
      "2018_H_plus_HF_genets":108,"2020_H_plus_HF_genets":79,
      "all_conditional_floral_beta_CI_lower_positive":True,
      "HF_fruiting_label_not_treated_as_causally_independent":True}))
if __name__=="__main__":main()
