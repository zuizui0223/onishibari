#!/usr/bin/env python3
"""v0.62: Sharp worst-case treatment-contrast identification regions for original
2022 *Pulsatilla alpina* one-flower bisexual randomized stamen-removal cohorts.

Based on v61 exact 81 enrolled-flower ID cohort and fixed-source SHA hashes.
No source original individual identifiers, genotype rows or exact localities.

Observed typed flowers have direct mature-seed counts and source inferred selfing
fraction. Followed-but-untyped flowers have OBSERVED mature seed counts, but
self-origin fractions unknown within [0,1]. Untracked flowers have unknown
mature seed counts AND unknown self-origin fractions: conditionally bound each
nonnegative "effective maternal seed equivalent" y in [0,K] where K is an
explicit HYPOTHETICAL cap, not a measured per-flower ovule limit. Under this
rectangular support and treatment-specific independent unknown values, min and
max group contrasts are sharp partial-identification intervals, not statistical
confidence intervals. Changing K, delta answers what assumptions would identify
the sign of the ORIGINAL 54-vs-12 assigned-cohort contrast.

Source 2022 population.csv does NOT contain pistil counts for the 17
untracked H flowers, hence no individually validated physical upper caps.
Pulsatilla uniovulate pistils provide theoretical y<=pistils, but missing
pistil counts make it inapplicable here without more data.
Known 2025 113/126 mature seed outcomes retained exactly, not imputed to zero.

No bootstrap here: the parent v61 flower-level + Beta typed seed CI is separate,
and these are deterministic *identification* bounds conditional on typed values.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

SOURCE=Path("data/2026-11-24_pulsatilla_selection_tipping_v61.json")
DEFAULT=Path("analysis_results_pulsatilla_identification_v62/pulsatilla_original_cohort_partial_identification_v62.json")
DELTA=(0.,.25,.5,.75,.95,1.)
CAPS=(0,5,10,15,20,25,30,50,75,100,150,200,400)
SRC_SIG={
 "census":"f03cfd1f7f7d4c829dc21691298fa7cbb0930f3719cf5cadeb4b8fc6bd85cbc2",
 "phenotypes":"e49fef29cefa6e1ed419500a0e312a49524f04876cca5eb21dda14657fe11a98",
 "fitness":"2cb1c78efe86128255dabed72d01799cfe71e10208386f55319f708a8535c141",
}

def load():
    v=json.loads(SOURCE.read_text())
    for key,hash in SRC_SIG.items():
        assert v["provenance"][key]["sha256"]==hash
    assert v["source_single_flower_bisexual_original_assigned"]==81
    assert v["source_genotyped_bisexual_confirmed_original_census_ID"]==59
    assert v["source_genotyped_bisexual_without_exact_original_census_ID"]==1
    for arm,n,n_typed,n_known,n_lost in [
      ("intact",54,42,3,9),("half",15,11,1,3),("all",12,6,1,5)
    ]:
        g=v["cohort_flow_by_arm"][arm]
        assert (g["original_assigned_single_flower_bisexual"],
                g["genotyped_mature_outcome"],
                g["mature_outcome_known_but_selfing_untyped"],
                g["no_matching_mature_outcome_in_source"])==(n,n_typed,n_known,n_lost)
    return v

def terms(v,delta):
    d=v["tipping_by_delta_original_assigned"][str(float(delta))]
    g=v["cohort_flow_by_arm"]
    assert d["typed_n"]["intact"]==42 and d["typed_n"]["all"]==6
    return {
      a:{
        "enrolled_n":g[a]["original_assigned_single_flower_bisexual"],
        "typed_effective_mature_seed_sum":d["typed_sum_"+a],
        "known_untyped_mature_seed_sum":g[a]["known_mature_outcome_untyped_harvest_seeds_sum"],
        "lost_no_mature_outcome_n":g[a]["no_matching_mature_outcome_in_source"]
      } for a in ("all","intact")
    }

def sharp(v,delta,cap_lost_all,cap_lost_intact):
    q=terms(v,delta)
    a,i=q["all"],q["intact"]
    # Partial identification of individual unknown selfing fraction
    # and missing mature seed totals. delta in [0,1] ensures nonnegativity.
    y_a_min=(a["typed_effective_mature_seed_sum"]+
       (1-delta)*a["known_untyped_mature_seed_sum"])/a["enrolled_n"]
    y_a_max=(a["typed_effective_mature_seed_sum"]+
       a["known_untyped_mature_seed_sum"]+
       a["lost_no_mature_outcome_n"]*cap_lost_all)/a["enrolled_n"]
    y_i_min=(i["typed_effective_mature_seed_sum"]+
       (1-delta)*i["known_untyped_mature_seed_sum"])/i["enrolled_n"]
    y_i_max=(i["typed_effective_mature_seed_sum"]+
       i["known_untyped_mature_seed_sum"]+
       i["lost_no_mature_outcome_n"]*cap_lost_intact)/i["enrolled_n"]
    low=y_a_min-y_i_max
    high=y_a_max-y_i_min
    assert low<=high+1e-9
    status="positive_for_all_admissible_missing" if low>1e-10 else (
      "negative_for_all_admissible_missing" if high< -1e-10 else "sign_not_identified")
    return {
      "delta":delta,
      "hypothetical_cap_effective_seeds_per_untracked_all_flower":cap_lost_all,
      "hypothetical_cap_effective_seeds_per_untracked_intact_flower":cap_lost_intact,
      "sharp_treatment_mean_difference_lower":low,
      "sharp_treatment_mean_difference_upper":high,
      "sign_under_assumed_missingness_cap":status,
      "all_assigned_mean_range":[y_a_min,y_a_max],
      "intact_assigned_mean_range":[y_i_min,y_i_max],
    }

def boundary(v,delta):
    q=terms(v,delta)
    a,i=q["all"],q["intact"]
    # Suppose all lost all-removal flowers could produce at least 0; the
    # worst-case lower contrast assigns them zero. Positive sign guaranteed
    # only if intact lost flower effective seed output is <= this threshold.
    amin=(a["typed_effective_mature_seed_sum"]+
        (1-delta)*a["known_untyped_mature_seed_sum"])/a["enrolled_n"]
    imax_at_zero=(i["typed_effective_mature_seed_sum"]+
        i["known_untyped_mature_seed_sum"])/i["enrolled_n"]
    if i["lost_no_mature_outcome_n"]==0:
        return None
    kpos=(amin-imax_at_zero)*i["enrolled_n"]/i["lost_no_mature_outcome_n"]
    # If all untracked intact flower successes are 0 and all known untyped
    # are maximally selfed (minimal effective seeds), negative sign is
    # guaranteed when all-removed missing cap is sufficiently low.
    amax_at_zero=(a["typed_effective_mature_seed_sum"]+
        a["known_untyped_mature_seed_sum"])/a["enrolled_n"]
    imin=(i["typed_effective_mature_seed_sum"]+
        (1-delta)*i["known_untyped_mature_seed_sum"])/i["enrolled_n"]
    kneg=(imin-amax_at_zero)*a["enrolled_n"]/a["lost_no_mature_outcome_n"]
    return {
      "delta":delta,
      "lower_bound_at_both_missing_caps_0":amin-imax_at_zero,
      "upper_bound_at_both_missing_caps_0":amax_at_zero-imin,
      "max_intact_missing_effective_seeds_for_guaranteed_positive":kpos if kpos>=0 else None,
      "max_all_missing_effective_seeds_for_guaranteed_negative":kneg if kneg>=0 else None,
      "guaranteed_positive_requires_untracked_intact_cap_below_criterion":kpos>0,
      "guaranteed_negative_requires_untracked_all_cap_below_criterion":kneg>0,
    }

def main():
    v=load()
    scenarios={}
    for d in DELTA:
        key=str(float(d))
        rows=[sharp(v,d,c,c) for c in CAPS]
        scenarios[key]={
          "boundaries":boundary(v,d),
          "symmetric_missing_outcome_hypothetical_cap_scenarios":rows,
          "asymmetric_lost_outcome_cap_grid_selected":[
            sharp(v,d,cap_a,cap_i)
            for cap_a in (0,25,50,100) for cap_i in (0,10,20,25,50,100)
          ],
        }
    assert scenarios["0.95"]["boundaries"]["max_intact_missing_effective_seeds_for_guaranteed_positive"]>21
    assert scenarios["0.95"]["boundaries"]["max_intact_missing_effective_seeds_for_guaranteed_positive"]<22
    assert sharp(v,.95,0,0)["sharp_treatment_mean_difference_lower"]>0
    assert sharp(v,.95,25,25)["sharp_treatment_mean_difference_lower"]<0
    assert sharp(v,.95,25,25)["sharp_treatment_mean_difference_upper"]>0
    return {
      "status":"v62_original_enrolled_Pulsatilla_female_treatment_partial_identification",
      "original_sources":v["provenance"],
      "published_source_census_has_pistil_ovule_counts_for_unmatched_flowers":False,
      "enrollment_and_outcome_link":{
        a:{k:v["cohort_flow_by_arm"][a][k] for k in (
            "original_assigned_single_flower_bisexual","genotyped_mature_outcome",
            "mature_outcome_known_but_selfing_untyped","no_matching_mature_outcome_in_source",
            "known_mature_outcome_untyped_harvest_seeds_sum"
            )} for a in ("intact","half","all")
      },
      "δ_sensitivity":scenarios,
      "interpretation":{
        "bounds_are_sharp_only_with_rectangular_untracked_flower_caps":True,
        "hypothetical_K_is_not_biological_upper_bound_measured_from_original_census":True,
        "known_mature_seed_counts_untyped_preserved_as_observed":True,
        "untyped_selfed_fraction_in_closed_0_1_interval":True,
        "untracked_outcomes_not_assumed_zero_except_explicit_K_0_case":True,
        "delta_is_assumed_not_new_survival_measure":True,
        "no_statistical_confidence_interval_for_enrolled_treatment_contrast":True,
        "v61_typed_flower_beta_bootstrap_is_separate":True,
        "this_is_not_daphne_residual_pistil_measurement":True,
        "no_original_plant_ID_alleles_or_coordinates_output":True,
      }
    }
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=DEFAULT)
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result=main()
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    summary={k:{"boundary":v["boundaries"],
       "caps":[{"K":x["hypothetical_cap_effective_seeds_per_untracked_all_flower"],
                "lower":x["sharp_treatment_mean_difference_lower"],
                "upper":x["sharp_treatment_mean_difference_upper"],
                "classification":x["sign_under_assumed_missingness_cap"]}
               for x in v["symmetric_missing_outcome_hypothetical_cap_scenarios"]]
    } for k,v in result["δ_sensitivity"].items()}
    print("PULSATILLA_V62_PARTIAL_IDENTIFICATION",json.dumps(summary,allow_nan=False,ensure_ascii=False),flush=True)
