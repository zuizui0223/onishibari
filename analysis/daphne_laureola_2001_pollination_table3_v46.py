#!/usr/bin/env python3
"""v0.46: Published original 2001 Daphne laureola factorial pollen experiment table.

Source: Alonso & Herrera (2001), American Journal of Botany 88:1016–1024,
doi:10.2307/2657083, original full article Table 3 (printed p.1020, PDF p.5)
https://www.plant-animal.es/pdfs/Alonso%26Herrera2001AmJBot.pdf

Each row is flowers nested in 5 hermaphrodite (H) or 3 female (F) plants.
Percentages are first-party reported to 1 decimal (including 17.0 and 5.0).
Counts are reconstructed as the unique integer compatible with rounding
and are NEVER independent plants or a new study.

The original authors report marginal self-pollen emasculation contrast
p=0.058, cross-pollen emasculation contrast p=0.58, and emphasize
emasculation side-effects, stigma contamination in cut control flowers,
and a lack of major sex differences after emasculation.
No same-sex (within-/different-flower) self-pollination contrast was
directly applied in the 2001 paper: 'self' donor pollen came from a
different flower on the same individual. This is not a direct measure
of Daphne pseudomezereum pistil viability or genetic paternity.
"""
from __future__ import annotations
from pathlib import Path
import json
from math import log

OUT=Path("analysis_results_daphne_pollen_clogging_prior_v46")
SRC_ARTICLE="https://www.plant-animal.es/pdfs/Alonso%26Herrera2001AmJBot.pdf"
DOI="10.2307/2657083"

# Hand transcribed from published Table 3, with group statuses and N.
PUBLISHED=[
 ("F","not_applicable","no_pollen",44,6.8),
 ("F","not_applicable","hand_cross",45,37.8),
 ("H","intact_cut_sham","no_pollen",88,17.0),
 ("H","intact_cut_sham","hand_cross",83,34.9),
 ("H","intact_cut_sham","hand_geitonogamous_self",88,18.2),
 ("H","emasculated_cut","no_pollen",80,5.0),
 ("H","emasculated_cut","hand_cross",77,40.3),
 ("H","emasculated_cut","hand_geitonogamous_self",98,33.7),
 ("H","bagged_untouched","no_pollen",41,0.0),
]

def compatible_count(n,percent):
    # Original first-party table prints percentage with 1-decimal precision.
    ks=[k for k in range(n+1) if abs(round(100*k/n,1)-percent)<1e-10]
    assert len(ks)==1,(n,percent,ks)
    return ks[0]

def main():
    OUT.mkdir(exist_ok=True)
    groups=[]
    for sex,emas,pollen,n,p in PUBLISHED:
        k=compatible_count(n,p)
        groups.append({
          "sex":sex,
          "emasculation":emas,
          "pollination":pollen,
          "flower_n":n,
          "reported_fruit_percent":p,
          "inferred_fruit_n_exact_given_table_rounding":k,
          "mature_fruit_rate":round(k/n,8),
          "unit_note":"Flowers nested within 5 H or 3 F mother plants; one article"
        })
    idx={(x["sex"],x["emasculation"],x["pollination"]):x for x in groups}
    def pick(emas,pollen):return idx[("H",emas,pollen)]
    comparisons={}
    for pollen in ("hand_cross","hand_geitonogamous_self"):
        h0=pick("intact_cut_sham",pollen)
        h1=pick("emasculated_cut",pollen)
        r0=h0["mature_fruit_rate"]
        r1=h1["mature_fruit_rate"]
        comparisons[pollen]={
          "fruit_success_intact":[h0["inferred_fruit_n_exact_given_table_rounding"],h0["flower_n"]],
          "fruit_success_emasculated":[h1["inferred_fruit_n_exact_given_table_rounding"],h1["flower_n"]],
          "descriptive_risk_difference_emas_minus_intact":round(r1-r0,8),
          "descriptive_relative_risk_emas_over_intact":round(r1/r0,8),
          "published_plant_adjusted_emasculation_contrast_p":
            0.58 if pollen=="hand_cross" else .058,
          "not_a_replicated_independent_flower_level_test":True
        }
    diff_self=comparisons["hand_geitonogamous_self"]["descriptive_risk_difference_emas_minus_intact"]
    diff_cross=comparisons["hand_cross"]["descriptive_risk_difference_emas_minus_intact"]
    out={
      "status":"source_verified_original_AJBot_Daphne_laureola_table3_v46",
      "doi":DOI,
      "source_article_pdf":SRC_ARTICLE,
      "original_table":"Table 3, printed page 1020, PDF zero-index page 4",
      "n_independent_H_plants":5,
      "n_independent_F_plants":3,
      "n_arm_rows":len(groups),
      "n_total_source_flower_rows_including_bag_control":sum(x["flower_n"] for x in groups),
      "treatment_groups":groups,
      "within_H_emasculation_comparisons":comparisons,
      "difference_of_descriptive_risk_differences_self_minus_cross":
         round(diff_self-diff_cross,8),
      "published_emas_vs_intact_self_p":.058,
      "published_emas_vs_intact_cross_p":.58,
      "conditional_interpretation":{
        "self_pollen_source_is_DIFFERENT_FLOWER_SAME_INDIVIDUAL":True,
        "no_within_flower_same_flower_self_hand_pollen_treatment":True,
        "published_original_plant_level_model_not_reproduced":True,
        "pollination_results_are_not_independent_of_hand_cut_and_contamination":True,
        "emasculation_can_change_reproductive_physiology_and_pollen_supply":True,
        "nonrandom_sex_assignment_no_causal_H_vs_F_comparison":True,
        "self_reproductive_success_not_measured_genetic_adult_fitness":True,
        "self_vs_cross_difference_can_reflect_inbreeding_depression":True,
        "no_original_Daphne_pseudomezereum_measurement":True,
        "not_new_original_experiment_or_independent_replication":True
      }
    }
    assert out["n_total_source_flower_rows_including_bag_control"]==644
    assert comparisons["hand_geitonogamous_self"]["fruit_success_intact"]==[16,88]
    assert comparisons["hand_geitonogamous_self"]["fruit_success_emasculated"]==[33,98]
    assert comparisons["hand_cross"]["fruit_success_intact"]==[29,83]
    assert comparisons["hand_cross"]["fruit_success_emasculated"]==[31,77]
    assert 0<diff_cross<diff_self
    (OUT/"daphne_laureola_2001_emasculation_fruit_table3_v46.json").write_text(
       json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_2001_POLLEN_CLOGGING_FACTORIAL_V46",json.dumps({
       "status":out["status"],"group_counts":len(groups),
       "plants_H":5,"plants_F":3,
       "self_result":comparisons["hand_geitonogamous_self"],
       "cross_result":comparisons["hand_cross"],
       "diff_in_diffs_descriptive":out["difference_of_descriptive_risk_differences_self_minus_cross"],
       "no_new_D_pseudomezereum_data":True},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
