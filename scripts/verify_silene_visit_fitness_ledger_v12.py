#!/usr/bin/env python3
"""Validate v0.12 individual-level Silene visitation/fitness ledger vs fresh raw fit.

Must run in SAME job after original source script. Confirms 17 new
non-independent summary records match output. Does not validate causal effects.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path

ROOT=Path("analysis_results_silene_visitors_vs_repeat_v12/silene2022_visitors_repeat_proxy_vs_genetic_fitness_v12.json")
LEDGER=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_silene_real_visits_vs_genetic_fitness_v12.md")

def main():
    data=json.loads(ROOT.read_text(encoding="utf-8"))
    rows=list(csv.DictReader(LEDGER.open(encoding="utf-8",newline="")))
    assert data["status"]=="success_original_2022_Silene_same_plant_visitation_fitness_audit"
    assert data["new_arrivals_cannot_be_identified_from_unmarked_visitors"] is True
    assert data["source_exposure_definitions"]["n_with_both_visitor_variables"]==168
    assert data["source_exposure_definitions"]["n_returned_contacts_proxy_positive"]==85
    assert data["source_exposure_definitions"]["source_attached_repeated_individual_bee_IDs"] is False
    assert REPORT.exists() and "Makino" in REPORT.read_text()
    assert len(rows)>=125
    assert len(rows)==len({r["record"] for r in rows})
    ledger={r["record"]:r for r in rows}
    count=0
    def check(name,sample,estimate,ci=None,p=None):
        nonlocal count
        assert name in ledger,name
        row=ledger[name]
        assert int(row["n"])==sample,(name,row["n"],sample)
        assert row["measured_floral_biomass"]=="False"
        assert row["primary_doi"]=="10.1098/rspb.2022.1987"
        for field,value in [("estimate",estimate),
                            ("ci_lower",ci[0] if ci else None),
                            ("ci_upper",ci[1] if ci else None),
                            ("p_unadjusted",p)]:
            if value is None:continue
            got=float(row[field])
            assert math.isclose(got,float(value),abs_tol=1e-5,rel_tol=1e-5),(name,field,got,value)
            count+=1
    M=data["comparisons_by_sex_and_endpoint"]["male_siring_MEMM_index"]
    F=data["comparisons_by_sex_and_endpoint"]["female_viable_seed_index"]
    MA=data["comparisons_by_sex_and_endpoint"]["male_genetic_mate_number"]
    FA=data["comparisons_by_sex_and_endpoint"]["female_genetic_mate_number"]
    mname={
      "v12_male_siring_total_visits_separate_beta":("plus_visits","lnVisits"),
      "v12_male_siring_visitor_number_separate_beta":("plus_visitor_number","lnVisitors"),
      "v12_male_siring_joint_visits_beta":("plus_both_visit_and_visitor_number","lnVisits"),
      "v12_male_siring_joint_visitor_number_beta":("plus_both_visit_and_visitor_number","lnVisitors")
    }
    for name,(model,term) in mname.items():
        v=M["models"][model]["terms"][term]
        check(name,84,v["beta"],v["95_HC3"],v["p"])
    for name,model,term in [
      ("v12_female_viable_seeds_total_visits_separate_beta","plus_visits","lnVisits"),
      ("v12_female_viable_seeds_observed_visitors_separate_beta","plus_visitor_number","lnVisitors")]:
        v=F["models"][model]["terms"][term]
        check(name,83,v["beta"],v["95_HC3"],v["p"])
    for name,dt,model,term,n in [
      ("v12_male_genetic_mates_observed_visits_separate_beta",MA,"plus_visits","lnVisits",83),
      ("v12_male_genetic_mates_observed_visitors_separate_beta",MA,"plus_visitor_number","lnVisitors",83),
      ("v12_female_genetic_mates_observed_visitors_separate_beta",FA,"plus_visitor_number","lnVisitors",82)
    ]:
        v=dt["models"][model]["terms"][term]
        check(name,n,v["beta"],v["95_HC3"],v["p"])
    for name,dt,model,n in [
      ("v12_male_siring_LOO_MSE_baseline",M,"true_flower_gamete_and_treatment",84),
      ("v12_male_siring_LOO_MSE_totalvisits",M,"plus_visits",84),
      ("v12_male_siring_LOO_MSE_observedvisitors",M,"plus_visitor_number",84),
      ("v12_male_siring_LOO_MSE_joint",M,"plus_both_visit_and_visitor_number",84),
      ("v12_female_viable_seeds_LOO_MSE_baseline",F,"true_flower_gamete_and_treatment",83),
      ("v12_female_viable_seeds_LOO_MSE_joint",F,"plus_both_visit_and_visitor_number",83)]:
        check(name,n,dt["models"][model]["LOO_plant_MSE"])
    check("v12_male_visitors_totalvisits_spearman",84,M["visitors_and_total_visits_spearman"])
    check("v12_female_visitors_totalvisits_spearman",83,F["visitors_and_total_visits_spearman"])
    a=M["models"]
    assert a["plus_visitor_number"]["LOO_plant_MSE"]<a["plus_both_visit_and_visitor_number"]["LOO_plant_MSE"]
    assert a["plus_visits"]["LOO_plant_MSE"]<a["plus_both_visit_and_visitor_number"]["LOO_plant_MSE"]
    assert F["models"]["true_flower_gamete_and_treatment"]["LOO_plant_MSE"]<F["models"]["plus_both_visit_and_visitor_number"]["LOO_plant_MSE"]
    assert len([k for k in ledger if k.startswith("v12_")])==17
    print(json.dumps({"status":"PASS_Silene_v12_original_source_ledger",
                      "source_plant_visitor_records":168,"evidence_rows":len(rows),
                      "new_v12_dependent_records_checked":17,"numeric_checks":count,
                      "male_visitors_total_visits_spearman":M["visitors_and_total_visits_spearman"],
                      "female_visitors_total_visits_spearman":F["visitors_and_total_visits_spearman"],
                      "source_sha256":hashlib.sha256(ROOT.read_bytes()).hexdigest(),
                      "ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
                      "cannot_identify_unique_bee_arrivals":True},indent=2))
if __name__=="__main__":main()
