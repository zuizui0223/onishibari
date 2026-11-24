#!/usr/bin/env python3
"""Offline invariant gate on v65 SHA-pinned original-source output."""
import argparse,json
from pathlib import Path
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,default=Path("data/2026-11-24_daphne_zero_state_models_v65.json"))
    x=json.loads(p.parse_args().actual.read_text())
    assert x["original_zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    assert x["source_original_2018_2020_adult_genets"]==492
    expected={(2018,"H"):(105,4765,53,14),
              (2018,"F"):(63,619,255,37),
              (2020,"H"):(71,2573,28,9),
              (2020,"F"):(57,719,205,32)}
    for z in x["sex_year"]:
        k=z["year"],z["morph_exact"]
        assert (z["genets"],z["flowers"],z["fruits"],z["fruit_positive_genets"])==expected[k]
        cv=z["genet_level_crossvalidation"]
        assert cv["folds"]==5 and cv["repeats"]==2
        assert cv["beta_binomial"]["avg_log_predictive_score_per_genet"]>cv["binomial"]["avg_log_predictive_score_per_genet"]
        if z["morph_exact"]=="H":
            assert cv["beta_binomial"]["avg_log_predictive_score_per_genet"]>=cv["zero_inflated_beta_binomial"]["avg_log_predictive_score_per_genet"]-1e-5
        assert cv["beta_binomial"]["heldout_genets"]==z["genets"]
    a=x["previous_microscopy_published_thesis_stage_evidence"]
    assert a["genets_H"]==6 and a["pistils_H"]==55 and a["genets_F"]==11 and a["pistils_F"]==11
    for z in x["within_source_between_year_morph_parameter_transfer"]["H"]:
        pred=z["models"]["beta_binomial"]
        assert abs(pred["target_predicted_zero_fraction"]-pred["target_observed_zero_fraction"])<.005
    print("v65 original-source SHA, 4-model CV, H/F census, source-microscopy and year-transfer aggregate gate PASS")
