#!/usr/bin/env python3
"""Verify v0.10 Sagittaria original-data report ledger against a fresh source run.

Must run AFTER sagittaria2026_original_sequential_reproduction_v09.py in SAME
GitHub Actions job. This avoids manually transcribing a source report into
claims that drift from the actual source version. Validates key estimates,
uncertainty, denominators, and environmental confounding declarations.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path

SRC=Path("analysis_results_sagittaria2026_sequential_fitness_v09/sagittaria_trifolia_2026_original_sequential_compensation_v09.json")
LEDGER=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_sagittaria_sequential_fruiting_genotype_matched_v10.md")

def main():
    src=json.loads(SRC.read_text(encoding="utf-8"))
    rows=list(csv.DictReader(LEDGER.open(encoding="utf-8",newline="")))
    keys={r["record"]:r for r in rows}
    assert len(keys)==len(rows), "duplicate evidence records"
    assert len(rows)>=98
    assert src["status"]=="success_original_figshare_sequential_fruiting_analysis"
    assert src["source_tables"]["inflorescence_rows"]==655
    assert src["source_tables"]["total_rows"]==109
    assert src["n_genotypes_total"]==47
    assert src["repro_plants_to_total_matched"]==109
    assert src["repro_plants_without_matched_total_seeds"]==4
    assert src["genotype_inconsistencies_on_matched_individuals"]==1
    assert src["treatment_site_confounded_no_factorial_pollen_or_light_test"] is True
    assert REPORT.exists()
    assert "site" in REPORT.read_text().lower()

    verified=[]
    def check(rec, n, estimate, ci=None, p=None, tol=1e-5):
        if rec not in keys:raise AssertionError("missing evidence row "+rec)
        row=keys[rec]
        if int(float(row["n"]))!=n:
            raise AssertionError(f"{rec}: n {row['n']} differs from source {n}")
        fields=[("estimate",estimate)]
        if ci is not None:fields.extend([("ci_lower",ci[0]),("ci_upper",ci[1])])
        if p is not None:fields.append(("p_unadjusted",p))
        for field,x in fields:
            got=float(row[field])
            if not math.isclose(got,float(x),abs_tol=tol,rel_tol=tol):
                raise AssertionError(f"{rec}: {field} got {got} expected {x}")
            verified.append(rec+":"+field)
        if row["measured_floral_biomass"]!="False":
            raise AssertionError("incorrect measured floral mass flag")
        if row["analysis_file"]!="results/2026-11-16_sagittaria_sequential_fruiting_genotype_matched_v10.md":
            raise AssertionError("source report is not the audited Sagittaria primary-data report")

    s=src["by_treatment"]
    for env in ("Net","Shade","Open"):
        t=s[env]
        rec="sagittaria_"+env.lower()+"_late_minus_first_fruit_fraction_matched"
        x=t["late_minus_first_fruit_fraction_same_plant"]
        check(rec, t["n_plant_pairs_with_first_and_later_female_flowers"],
              x["mean"],x["genotype_bootstrap95"])
        x=t["plant_level_total_seeds_TFS"]
        check("sagittaria_"+env.lower()+"_whole_plant_mean_seeds",
              x["n_plants"],x["mean"],x["genotype_bootstrap95"])
    t=s["Net"]["proportion_plants_initially_zero_but_later_fruit"]
    assert t["n_initial_zero"]==27 and t["n_rescued_later"]==21
    check("sagittaria_net_first_zero_later_any_fruit",27,t["prob_recovery_given_initial_zero"])
    root=src["cross_site_common_genotype_contrasts"]
    for env,field,suffix in [
        ("Net_minus_Shade","first_fruit_fraction","sagittaria_net_shade_shared_genotype_first_fruit_fraction"),
        ("Net_minus_Shade","later_fruit_fraction","sagittaria_net_shade_shared_genotype_later_fruit_fraction"),
        ("Net_minus_Shade","late_minus_first_fraction","sagittaria_net_shade_shared_genotype_recovery_difference"),
        ("Net_minus_Shade","TFS","sagittaria_net_shade_shared_genotype_total_seed_difference"),
        ("Net_minus_Shade","inflorno","sagittaria_net_shade_shared_genotype_inflorescence_number_difference"),
        ("Net_minus_Open","TFS","sagittaria_net_open_shared_genotype_total_seed_difference"),
        ("Shade_minus_Open","TFS","sagittaria_shade_open_shared_genotype_total_seed_difference")]:
        x=root[env][field]
        check(suffix,x["n_shared_genotypes"],x["mean_a_minus_b"],
              x["genotype_pair_bootstrap95"],x["within_genotype_signflip_p"])
        assert x["site_confounded_observational_contrast"] is True
    print(json.dumps({
        "status":"PASS_v0.10_original_source_ledgers",
        "source_status":src["status"],
        "stage_evidence_rows":len(rows),
        "numeric_assertions":len(verified),
        "source_json_sha256":hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
        "n_original_plant_rows":src["source_tables"]["total_rows"],
        "n_genotypes":src["n_genotypes_total"],
        "causal_pollination_effect_not_claimed":True,
        "no_record_ids_exported":True,
    },indent=2))
if __name__=="__main__":main()
