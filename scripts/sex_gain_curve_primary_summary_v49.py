#!/usr/bin/env python3
"""v0.49: verify PUBLICATION-SUMMARY anchor values, never claim raw-data refitting.

Input: fixed primary-published exponents and previous source-verified
Daphne v0.30 reanalysis. Output: 2**b descriptive translations only.
Original dryad and Zenodo raw rows are NOT retrieved by this audit.
"""
import csv
import json
from math import isclose
from pathlib import Path

SOURCE = Path("data/2026-11-20_sex_gain_curve_primary_anchors_v49.csv")
OUT = Path("analysis_results_sex_gain_curve_v49")
EXPECTED = {
    "ambrosia2022_male": (0.85, 0.25, 1.46),
    "pulsatilla2025_male": (0.50, 0.00, 1.01),
    "pulsatilla2025_female": (1.26, 1.06, 1.46),
    "daphne2018_male": (0.6250, 0.4583, 0.8326),
    "daphne2020_male": (0.8956, 0.5236, 1.2621),
}
def main():
    with SOURCE.open(encoding="utf8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len(EXPECTED), "unexpected source-anchor row count"
    assert {r["record_id"] for r in rows} == set(EXPECTED)
    assert len({r["study_population_cluster"] for r in rows}) == 3
    assert len({r["primary_article_doi"] for r in rows}) == 3
    output = []
    for row in rows:
        record_id = row["record_id"]
        reported = tuple(float(row[k]) for k in ("estimate_b", "ci_lower", "ci_upper"))
        assert all(isclose(x, y, abs_tol=1e-9) for x,y in zip(reported, EXPECTED[record_id]))
        b, low, high = reported
        assert low <= b <= high and all(v >= 0 for v in reported)
        assert row["estimand"] in {"published_power_exponent", "secondary_conditional_choice_exponent"}
        output.append({
            "record_id": record_id,
            "source_population_cluster": row["study_population_cluster"],
            "estimand": row["estimand"],
            "b_published_or_existing_verified": b,
            "ci95_published_or_existing_verified": [low, high],
            "twofold_multiplier_point": round(2**b, 5),
            "twofold_multiplier_interval_monotone_transform": [round(2**low, 5), round(2**high, 5)],
            "source": row["primary_article_doi"],
            "raw_data_refit": False,
        })
    assert isclose(8.05/4.44, 1.813063063063, rel_tol=1e-12)
    assert output[0]["twofold_multiplier_point"] == 1.8025
    assert 2**1.26 > 2.0 > 2**0.50
    result = {
        "status": "publication_summary_only_v49",
        "original_source_rows_retrieved": False,
        "original_model_refit": False,
        "new_independent_publication_anchors": 2,
        "independent_population_study_clusters_total": 3,
        "unique_parent_article_dois": 3,
        "parent_articles_are_not_independent_exponent_rows": True,
        "ambrosia_adjusted_high_vs_low_treatment_sired_seed_means": [8.05,4.44],
        "ambrosia_adjusted_high_low_point_ratio_no_CI": round(8.05/4.44,5),
        "rows": output,
        "daphne_choice_weight_not_direct_sired_seed_multiplier": True,
        "no_D_pseudomezereum_data": True,
    }
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"sex_gain_curve_published_anchor_transform_v49.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("SEX_GAIN_CURVES_PUBLISHED_ANCHORS_V49", json.dumps(result,ensure_ascii=False))
if __name__ == "__main__":
    main()
