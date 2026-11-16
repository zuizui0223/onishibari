#!/usr/bin/env python3
"""Check v0.11 manuscript-grade source ledger and curve claims against fresh data.

To be invoked AFTER joffard2026_pollen_carryover_shape_cv_v11.py in the SAME
GitHub Actions job. A pass validates the nine new data-dependent evidence rows;
it does not create an independent biological replication or validate causality.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path
SRC=Path("analysis_results_pollen_carryover_shape_v11/pollen_carryover_source_heldout_bee_shapes_v11.json")
LEDGER=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_pollen_carryover_fresh_bouts_v11.md")

def main():
    d=json.loads(SRC.read_text(encoding="utf-8"))
    r=list(csv.DictReader(LEDGER.open(encoding="utf-8",newline="")))
    assert d["status"]=="success_original_held_out_bee_shape_models"
    assert d["source_QA"]["n_bees"]==25
    assert d["source_QA"]["n_original_stigmas"]==186
    assert d["source_QA"]["n_unique_bee_position_cells"]==185
    assert d["source_QA"]["n_duplicate_bee_positions"]==1
    assert d["held_out_bee_and_plant_ids_not_exported"] is True
    assert REPORT.exists() and "25" in REPORT.read_text() and "zero" in REPORT.read_text().lower()
    assert len(r)>=107
    records={row["record"]:row for row in r}
    if len(r)!=len(records):raise AssertionError("duplicate ledger record")
    checks=[]
    def check(name,est,lower=None,upper=None,p=None,rel=1e-5):
        if name not in records:raise AssertionError("missing "+name)
        row=records[name]
        if int(row["n"])!=25:raise AssertionError("must be 25 independent bees, not 185 stigmas")
        for k,v in [("estimate",est),("ci_lower",lower),("ci_upper",upper),("p_unadjusted",p)]:
            if v is None:continue
            got=float(row[k])
            if not math.isclose(got,float(v),rel_tol=rel,abs_tol=rel):
                raise AssertionError(f"{name}/{k} current={got} source={v}")
            checks.append(name+":"+k)
        assert row["measured_floral_biomass"]=="False"
        assert "joffard" in name
    z=d["ordinal_zero_grain_risk"]
    check("joffard2026_stigma_zero_grains_ordinal_GEE_OR",z["gee_logit_or_per_additional_visit"],
          *z["gee_robust_OR_95"],z["p"])
    s=d["candidates"]
    curve=s["power_log_position"]
    check("joffard2026_within_bee_logpos_curve_beta",curve["coefficient_per_basis"][0],
          *curve["cluster95_by_basis"][0])
    for suffix,j in [("first_to_second_drop_beta",0),("after_second_slope_beta",1)]:
        src=s["continuous_kink_at_2"]
        check("joffard2026_kink2_"+suffix,src["coefficient_per_basis"][j],
              *src["cluster95_by_basis"][j],src["cluster_p_by_basis"][j])
    for key,tag in [
      ("exponential_position","exponential"),
      ("power_log_position","power"),
      ("continuous_kink_at_2","kink2")]:
        check("joffard2026_bee_holdout_MSE_"+tag,s[key]["loo_bee_mean_MSE"])
    for key,tag in [
       ("exponential_position_vs_power_log_position","exponential_minus_power"),
       ("power_log_position_vs_continuous_kink_at_2","power_minus_kink2")]:
        c=d["pairwise_cv"][key]
        check("joffard2026_CV_"+tag+"_MSE",c["mean_bee_mse_A_minus_B"],
              *c["paired_bee_bootstrap_95_mse_difference"])
        assert c["paired_bee_bootstrap_95_mse_difference"][0]<0<c["paired_bee_bootstrap_95_mse_difference"][1]
    assert d["selected_by_exploratory_cv"]=="power_log_position"
    assert d["source_pollen_deposition_not_compatible_tubes_or_seeds"] is True
    print(json.dumps({
      "status":"PASS_original_source_pollen_carryover_CV_v11",
      "original_unique_bee_positions":d["source_QA"]["n_unique_bee_position_cells"],
      "independent_bees":25,"ledger_rows":len(r),
      "dependent_analysis_rows_checked":9,"numeric_checks":len(checks),
      "source_output_sha256":hashlib.sha256(SRC.read_bytes()).hexdigest(),
      "evidence_ledger_sha256":hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
      "no_unique_exhaustion_kink_claimed":True,
      "same_study_rows_not_independent_replications":True
    },indent=2))
if __name__=="__main__":main()
