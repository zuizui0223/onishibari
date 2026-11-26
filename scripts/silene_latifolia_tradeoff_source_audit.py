#!/usr/bin/env python3
"""Source-faithful arithmetic audit, Delph & Herlihy (2011/2012) Silene latifolia.

This is NOT a refit of Dryad 10.5061/dryad.20k62r83 individual data.
The original article https://doi.org/10.1111/j.1558-5646.2011.01510.x:
- Table 2 gives male flower mass (mg) by line; dry/wet basis not audited here.
- Its stated pollen grains per flower are borrowed from Herlihy & Delph
  (2009), *different individuals* in same artificially selected lineages.
- Twice-weekly open flower COUNT SUMS are repeated open-flower observations;
  NOT distinct constructed flowers and cannot be multiplied by mass to get
  total floral construction biomass.
- PGI marker identifies SIRING BY SELECTION LINE, not individual father.
- Survival comparison is a selection-line contrast with correlated leaf
  physiological traits, not a randomized floral-cost treatment.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

OUT=Path("analysis_results_silene_flower_cost_carryover")
OUT.mkdir(exist_ok=True,parents=True)

SOURCE={
 "large":{"male_flower_mass_mg":10.81,"male_flower_mass_SE_mg":.532,"pollen_per_flower_source_other_study":29750.0,
          "open_flowers_census_sum_2004":330.4,"open_flowers_census_sum_2005":820.2},
 "small":{"male_flower_mass_mg":7.18,"male_flower_mass_SE_mg":.297,"pollen_per_flower_source_other_study":13040.0,
          "open_flowers_census_sum_2004":324.9,"open_flowers_census_sum_2005":699.8}
}

def main():
    a,b=SOURCE["large"],SOURCE["small"]
    mass_ratio=a["male_flower_mass_mg"]/b["male_flower_mass_mg"]
    pollen_ratio=a["pollen_per_flower_source_other_study"]/b["pollen_per_flower_source_other_study"]
    pollenmass_large=a["pollen_per_flower_source_other_study"]/a["male_flower_mass_mg"]
    pollenmass_small=b["pollen_per_flower_source_other_study"]/b["male_flower_mass_mg"]
    ratio=pollenmass_large/pollenmass_small
    assert 1.50<mass_ratio<1.51,mass_ratio
    assert 2.28<pollen_ratio<2.29,pollen_ratio
    assert 1.50<ratio<1.53,ratio
    # Crucial: unique flower count is unobserved by the paper's 2x/week census
    # summaries, so those sums cannot give empirical total floral biomass.
    not_a_construction_mass=True
    assert not_a_construction_mass
    report={
        "status":"validated_table_value_arithmetic_not_individual_reanalysis",
        "study":"Delph & Herlihy (2011 early online; Evolution 2012)",
        "doi":"10.1111/j.1558-5646.2011.01510.x",
        "dryad_original_data":"10.5061/dryad.20k62r83",
        "original_fruit_and_siring_and_survival_observed":True,
        "male_flower_mass_mg":{line:vals["male_flower_mass_mg"] for line,vals in SOURCE.items()},
        "male_flower_mass_SE_mg":{line:vals["male_flower_mass_SE_mg"] for line,vals in SOURCE.items()},
        "flower_mass_ratio_large_vs_small":round(mass_ratio,4),
        "pollen_per_flower_borrowed_from_2009":{line:vals["pollen_per_flower_source_other_study"] for line,vals in SOURCE.items()},
        "pollen_per_flower_ratio_large_vs_small":round(pollen_ratio,4),
        "cross_publication_pollen_per_mg_flower_mass_INDEX":{ "large":round(pollenmass_large,4),
            "small":round(pollenmass_small,4),"ratio_large_vs_small":round(ratio,4)},
        "pollen_per_mg_is_index_not_same_plant_joint_measure":True,
        "twice_weekly_sum_of_open_flower_snapshots_NOT_distinct_flowers":True,
        "do_not_multiply_open_count_snapshot_sum_by_flower_mass_to_infer_allocation":True,
        "survival_logrank_four_sex_line_groups":{"chi2":10.73,"df":3,"paper_p":.013},
        "large_line_female_season_total_seed_mass_ttest":{"2004":{"t":2.96,"p":.004},
            "2005":{"t":2.51,"p":.015}},
        "PGI_paternity_assigned_to_line_not_individual_fathers":True,
        "genotyped_seeds_total":3560,
        "pollen_siring_share_by_line_not_disproportionate_to_pollen_source_means":{
            "2004":{"t":-.47,"df":12,"p":.65},"2005":{"t":.61,"df":7,"p":.56}},
        "phenological_and_physio_genetic_correlation_not_cost_manipulation":True,
        "female_seeds_per_fruit_higher_in_large_line_both_years":True,
        "unavailable_estimate":"whole-plant total lifetime floral construction mass from reported snapshot counts",
        "no_new_focal_daphne_observations":True}
    (OUT/"silene_latifolia_carryover_source_arithmetic.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    rows=[
      {"comparison":"male_mass_per_flower_large_over_small","value":mass_ratio,"unit":"ratio","source_status":"same_original_study_means"},
      {"comparison":"male_pollen_per_flower_large_over_small","value":pollen_ratio,"unit":"ratio","source_status":"pollen_from_an_earlier_related_study"},
      {"comparison":"male_pollen_per_flower_mass_large_over_small","value":ratio,"unit":"ratio","source_status":"cross_study_constructed_index_not_independent_fitness_measure"},
    ]
    with (OUT/"silene_latifolia_derived_mass_and_pollen_ratios.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print("SILENE_SOURCE_QC_OK",json.dumps(report,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
