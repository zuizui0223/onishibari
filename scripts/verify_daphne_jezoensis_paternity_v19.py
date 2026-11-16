#!/usr/bin/env python3
"""Contract checks for source-only v0.19 Daphne male distance evidence.
No independent novel experiment inferred from source reproductions.
"""
import json, math
from pathlib import Path

P=Path("analysis_results_daphne_paternity_family_v19/shibata_kudo2023_paternity_maternal_cluster_v19.json")
def main():
    j=json.loads(P.read_text())
    assert j["status"]=="original_zenodo_2023_paternity_secondary_finished"
    a=j["joint_source_integrity"];b=j["coordinate_quality"]
    assert a["n_all_seeds"]==584
    assert a["n_assigned_95"]==470
    assert b["n_95_assigned_with_both_coordinates"]==470
    assert a["n_unmatched_source_treatment"]==0
    assert a["n_dnaid_duplicate_rows"]==0
    expected={
      (2018,"Natural"):(132,30,88),
      (2018,"Day open"):(32,11,26),
      (2018,"Night open"):(45,10,41),
      (2020,"Natural"):(163,32,116),
      (2020,"Day open"):(52,7,42),
      (2020,"Night open"):(46,5,44)
    }
    obs={(p["year"],p["info"]):p for p in j["per_group"]}
    assert set(obs)==set(expected)
    for key,(seeds,mothers,near) in expected.items():
        p=obs[key]
        assert p["assigned_seed_with_distance"]==seeds,(key,p)
        assert p["distinct_mothers"]==mothers,(key,p)
        assert math.isclose(p["offspring_equal_frac_within_10m"],near/seeds,rel_tol=1e-12)
    pair={p["year"]:p for p in j["paired_day_night_maternal_families"]}
    assert pair[2018]["matched_mothers"]==9
    assert pair[2020]["matched_mothers"]==4
    assert pair[2018]["n_night_higher"]==1 and pair[2018]["n_day_higher"]==2 and pair[2018]["n_equal"]==6
    assert pair[2020]["n_night_higher"]==3 and pair[2020]["n_day_higher"]==0 and pair[2020]["n_equal"]==1
    assert j["not_a_direct_fitness_mediation_from_v18_fruit_stages"]
    assert j["no_anonymized_geolocations_genotypes_or_plant_ids_released"]
    print(json.dumps({"status":"PASS_DAPHNE_PARENTAGE_SOURCE_V19",
      "source_assigned_seeds":470,"source_total_seeds":584,
      "reproduced_author_short_distance_group_counts":len(expected),
      "paired_mothers_2018":9,"paired_mothers_2020":4,
      "effect_estimates_are_secondary_and_conditional":True
    }))
if __name__=="__main__":main()
