#!/usr/bin/env python3
"""Verify source original father transfer & sire-exclusion numbers and caveats v38."""
import json,math
from pathlib import Path
P=Path("analysis_results_daphne_crossyear_paternity_v38/daphne_fatherhood_year_transfer_and_sire_exclusion_v38.json")
REPORT=Path("results/2026-11-17_daphne_year_transfer_paternity_sire_influence_v38.md")
def close(x,y):return math.isclose(float(x),float(y),rel_tol=.001,abs_tol=.001)
def main():
    x=json.loads(P.read_text())
    assert x["status"]=="v38_Daphne_2018_2020_crossyear_father_choice_and_sire_influence"
    assert x["source_doi"]=="10.5281/zenodo.8253921"
    assert x["distance_lambda_m"]==10
    assert len(x["groups"])==2
    expect={
      "H":{
        (2018,2020):(.409776,206,33,.180479),
        (2020,2018):(.238285,180,39,.037354)
      },
      "H,HF":{
        (2018,2020):(.3699,214,35,.149005),
        (2020,2018):(.241328,183,40,.039036)
      },
    }
    influence={"H":{2018:(56,.596941,.677755),2020:(35,.715913,.982648)},
               "H,HF":{2018:(57,.601653,.679089),2020:(38,.698786,.996968)}}
    for group in x["groups"]:
        sex=",".join(group["donor_sex_codes"])
        assert sex in expect
        for tr in group["cross_year_transfers"]:
            k=(tr["train_year"],tr["test_year"])
            gain,child,mother,cilow=expect[sex][k]
            assert close(tr["gain_over_distance_only_nats_per_child"],gain)
            assert tr["n_assigned_offspring"]==child
            assert tr["n_heldout_maternal_families"]==mother
            assert close(tr["fixed_train_beta_test_mother_boot_95%"][0],cilow)
            assert tr["train_only_model_applied_without_test_year_refit"]
            assert tr["not_an_independent_site_or_species_holdout"]
        for v in group["one_observed_sire_genet_leaveout_sensitivity"]:
            sirecount,lo,hi=influence[sex][v["year"]]
            assert v["n_influential_source_sire_genets_tested"]==sirecount
            assert v["n_different_sire_exclusion_estimates_positive"]==sirecount
            assert close(v["lowest_beta_after_dropping_one_sire"],lo)
            assert close(v["highest_beta_after_dropping_one_sire"],hi)
            assert v["lowest_beta_after_dropping_one_sire"]>0
            assert v["not_independent_sire_resampling"] and v["sire_exclusion_changes_conditioned_risk_set"]
    assert all(x["limits"][k] for k in (
      "cross_year_transfer_is_retrospective_not_pre_registered_forecast",
      "same_local_population_both_years_not_external_biogeographic_replication",
      "unmeasured_sire_spatial_resource_and_compatibility_confounding",
      "not_focal_Daphne_pseudomezereum"))
    assert REPORT.exists()
    txt=REPORT.read_text()
    for kw in ("0.4098","0.2383","0.5969","0.7159","同じ調査集団","遺伝的適合性"):
        assert kw in txt
    print(json.dumps({"status":"PASS_DAPHNE_ORIGINAL_CROSSYEAR_AND_SIRE_V38",
       "2018_H_father_genets_tested":56,
       "2020_H_father_genets_tested":35,
       "both_cross_year_directions_positive":True,
       "all_single_sire_exclusion_betas_positive":True,
       "never_claim_external_species_validation":True}),flush=True)
if __name__=="__main__":main()
