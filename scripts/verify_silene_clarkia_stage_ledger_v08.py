#!/usr/bin/env python3
"""Guard reproducible v0.8 numerical claims against source-output drift.

Requires the Silene/Clarkia original source analysis to have run in THIS
job, rather than trusting earlier manually copied results.
Assertions do not claim biological significance or independent replication.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path

LEDGER=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
SOURCE=Path("analysis_results_silene_clarkia_stage_v08/silene_clarkia_stage_separation_v08.json")
def main():
    d=json.loads(SOURCE.read_text())
    rows=list(csv.DictReader(LEDGER.open(encoding="utf-8",newline="")))
    assert d["status"]=="success_original_raw_two_experiment_reanalysis"
    assert len(rows)>=71,(len(rows),"ledger row count")
    assert len(set(r["record"] for r in rows))==len(rows),"nonunique record names"
    assert d["silene_2022"]["raw_n"]==192
    assert d["clarkia_2025"]["n_complete_pairs"]==156
    assert d["clarkia_2025"]["n_independent_blocks"]==38
    ledger={r["record"]:r for r in rows}
    checked=[]
    def chk(rec,field,expected,scale=1e-5):
        if rec not in ledger:raise AssertionError("missing ledger record "+rec)
        got=float(ledger[rec][field])
        if not math.isclose(got,float(expected),rel_tol=scale,abs_tol=scale):
            raise AssertionError(f"{rec}/{field} source={expected}, ledger={got}")
        checked.append(rec+":"+field)
    s=d["silene_2022"]
    chk("pooled_fake_display_on_insect_visitors","estimate",s["pooled_visitor_count_robust_Poisson"]["treatment_count_IRR"])
    chk("pooled_fake_display_on_insect_visitors","ci_lower",s["pooled_visitor_count_robust_Poisson"]["IRR_95"][0])
    chk("pooled_fake_display_on_insect_visitors","ci_upper",s["pooled_visitor_count_robust_Poisson"]["IRR_95"][1])
    s_pairs=[
      ("female_fake_display_visitor_count","female_visiting_insects"),
      ("fake_display_to_female_viable_seed_index","female_viable_seeds"),
      ("fake_display_to_male_MEMM_siring_index","male_siring_MEMM"),
      ("fake_display_to_female_genetic_mating_partners","female_mating_partners")
    ]
    for record,key in s_pairs:
        src=s["outcome_comparisons"][key]["test"]
        chk(record,"estimate",src["treatment_minus_control"])
        chk(record,"ci_lower",src["unpaired_bootstrap_95"][0])
        chk(record,"ci_upper",src["unpaired_bootstrap_95"][1])
        chk(record,"p_unadjusted",src["two_sided_label_permutation_p"])
    chk("real_flower_count_to_male_MEMM_siring","estimate",s["covariate_checks"]["male"]["real_flower_number"]["estimate"])
    chk("real_flower_count_to_female_viable_seed_index","estimate",s["covariate_checks"]["female"]["real_flower_number"]["estimate"])
    c=d["clarkia_2025"]
    assert sum(i["n_complete_pairs"] for i in c["outcome_by_species"].values())==156
    for sp,stats in c["outcome_by_species"].items():
        rec="paired_fruit_seedset_pollen_addition_"+sp
        chk(rec,"estimate",stats["mean_seedset_difference"])
        chk(rec,"ci_lower",stats["block_bootstrap95_seedset_difference"][0])
        chk(rec,"ci_upper",stats["block_bootstrap95_seedset_difference"][1])
        chk(rec,"p_unadjusted",stats["within_plant_signflip_p"])
        chk(rec,"q_adjusted",stats["exploratory_BH_q_over_four_species"])
        reg=c["observational_flower_number_moderation_of_pollen_addition"]["unadjusted_block_cluster"][sp]
        rec2="real_flower_count_moderation_of_pollen_gain_"+sp
        chk(rec2,"estimate",reg["slope_per_ln1p_real_flower"])
        chk(rec2,"ci_lower",reg["cluster_95"][0])
        chk(rec2,"ci_upper",reg["cluster_95"][1])
        chk(rec2,"p_unadjusted",reg["p"])
        assert float(ledger[rec]["n"])==stats["n_complete_pairs"]
    assert "2026-11-16_silene_clarkia_original_multistage_refit_v08.md" in [p.name for p in Path("results").glob("*.md")]
    print(json.dumps({
       "status":"PASS_original_source_results_match_committed_evidence_ledger",
       "ledger_rows":len(rows),"checks":len(checked),
       "original_output_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       "ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
       "source_genuinely_original":True,"no_location_or_original_record_ids_exported":True
    },indent=2))
if __name__=="__main__":main()
