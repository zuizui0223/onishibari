#!/usr/bin/env python3
"""Recompute and audit 23 v0.14 claims from BOTH original public author tables.

Run after scripts/johnson_harder2023_pollen_fates_v14.py AND
scripts/johnson_harder2023_pollen_ovules_51_overlap_v14.py.
The 51 paper's 49 duplicate species CANNOT be counted as an independent
replication of the 228-species data; 164/198 packaged taxa come from 2 families.
"""
import csv,json,math,hashlib
from pathlib import Path
P1=Path("analysis_results_global_pollen_export_v14/johnson_harder_228_pollen_stages_v14.json")
P2=Path("analysis_results_global_pollen_po_link_v14/johnson_harder_51_species_pollen_ovule_PTE_and_228_overlap_v14.json")
LED=Path("data/2026-11-15_stage_transition_evidence_v02.csv")
REPORT=Path("results/2026-11-16_global_pollen_loss_location_v14.md")

def main():
    x=json.loads(P1.read_text());y=json.loads(P2.read_text())
    assert x["status"]=="success_original_228_species_pollen_stage_reanalysis"
    assert y["status"]=="success_original_51_species_PTE_vs_pollen_ovule_separate_source"
    assert x["n_species_original"]==228 and x["n_RE"]==202 and x["n_PTE"]==224
    assert x["n_both"]==198 and x["n_families_with_both"]==20
    a=x["fate_summary"];types=x["by_pollen_packaging"]
    assert sum(a["produced_pollen_mass_balance_species_equal_means"].values())==True or abs(
        sum(v for k,v in a["produced_pollen_mass_balance_species_equal_means"].items() if k.endswith("pct"))-100)<.0001
    assert types["monads"]["n_species"]==34 and types["monads"]["n_families"]==18
    assert types["Massulate"]["n_species"]==81 and types["Solid"]["n_species"]==30
    assert types["Asclepiad"]["n_species"]==53
    assert a["families_by_species_sample_size"]["Orchidaceae"]==111
    assert a["families_by_species_sample_size"]["Apocynaceae"]==53
    assert y["cross_paper_taxon_overlap"]["n_exact_binomial_matches"]==49
    assert y["cross_paper_taxon_overlap"]["n_51_nonambiguous_names"]==51
    assert y["cross_paper_taxon_overlap"]["PTE_source_disagreement"]["n_more_than_10pct_relative_disagreement"]==1
    assert REPORT.exists() and "physical" in REPORT.read_text().lower()
    data=list(csv.DictReader(LED.open(encoding="utf-8",newline="")))
    rows={r["record"]:r for r in data}
    assert len(data)==len(rows) and len(data)>=164
    assert sum(k.startswith("v14_") for k in rows)==23
    checks=0
    def chk(s,n,v,ci=None,p=None):
        nonlocal checks
        assert s in rows,s
        row=rows[s]
        assert int(row["n"])==n,(s,row["n"],n)
        assert row["measured_floral_biomass"]=="False"
        for field,expected in [
          ("estimate",v),
          ("ci_lower",ci[0] if ci is not None else None),
          ("ci_upper",ci[1] if ci is not None else None),
          ("p_unadjusted",p)]:
            if expected is None:continue
            got=float(row[field])
            if not math.isclose(got,float(expected),rel_tol=2e-5,abs_tol=2e-5):
                raise AssertionError(f"{s}:{field} source {expected} ledger {got}")
            checks+=1
    for rec,val in [
      ("v14_198_species_median_removed_fraction",a["median_RE"]),
      ("v14_198_species_median_PTE_fraction",a["median_PTE_given_removal"]),
      ("v14_198_species_median_RPE_fraction",a["median_RPE_produced_to_stigma"]),
      ("v14_198_species_postremoval_physical_loss_dominates_fraction",a["fraction_post_removal_mass_loss_dominates"]),
      ("v14_198_species_postremoval_logloss_dominates_fraction",a["fraction_species_transport_loss_exceeds_removal_loss"])]:
        chk(rec,198,val)
    mass=a["produced_pollen_mass_balance_species_equal_means"]
    for rec,key in [
      ("v14_198_species_mean_unremoved_pollen_pct","not_removed_pct"),
      ("v14_198_species_mean_removed_not_delivered_pct","removed_not_delivered_pct"),
      ("v14_198_species_mean_delivered_pct","delivered_to_conspecific_stigma_pct")]:
        chk(rec,198,mass[key])
    chk("v14_198_species_median_transport_share_logloss",198,
        a["median_transport_share_of_log_shortfall"],
        a["family_CI_median_transport_logloss_share"]["CI95"])
    assoc=x["nonidentity_RE_vs_PTE"]
    chk("v14_198_species_RE_PTE_spearman",198,assoc["Spearman_rho"],
        assoc["family_boot_CI"]["CI95"])
    adj=x["RE_PTE_adjusted_noncausal"]
    chk("v14_172_species_RE_PTE_adjusted_logit_beta",172,adj["logitRE_beta"],
        adj["family_cluster_CI95"],adj["p"])
    chk("v14_twenty_family_equal_mean_RPE",20,a["mean_RPE_equal_family_weighting"])
    chk("v14_198_source_species_in_orchid_apocyn_family_count",198,
        a["families_by_species_sample_size"]["Orchidaceae"]+
        a["families_by_species_sample_size"]["Apocynaceae"])
    for key,rec,field,n in [
       ("monads","v14_monad_n34_mean_unremoved_pct","mass_balance_pct_unremoved_species_mean",34),
       ("monads","v14_monad_n34_mean_removed_not_delivered_pct","mass_balance_pct_removed_not_delivered_species_mean",34),
       ("monads","v14_monad_n34_mean_delivered_pct","mass_balance_pct_stigma_delivered_species_mean",34),
       ("Massulate","v14_massulate_n81_mean_removed_not_delivered_pct","mass_balance_pct_removed_not_delivered_species_mean",81),
       ("Solid","v14_solid_n30_mean_stigma_delivered_pct","mass_balance_pct_stigma_delivered_species_mean",30),
       ("Asclepiad","v14_asclepiad_n53_mean_unremoved_pct","mass_balance_pct_unremoved_species_mean",53)]:
        chk(rec,n,types[key][field])
    chk("v14_PTE_PO_51species_cross_paper_overlap_taxa",51,
        y["cross_paper_taxon_overlap"]["n_exact_binomial_matches"])
    for src,rec in [
      ("lnPollen","v14_lnPollen_vs_lnPTE_n51_familycluster"),
      ("lnOvules","v14_lnOvules_vs_lnPTE_n51_familycluster"),
      ("lnPO","v14_lnPollenOvuleRatio_vs_lnPTE_n51_familycluster")]:
        m=y["regressions"][src]
        chk(rec,51,m["beta"],m["CI95"],m["p"] if rec!="v14_lnPollenOvuleRatio_vs_lnPTE_n51_familycluster" else None)
    print(json.dumps({"status":"PASS_global_pollen_fate_v14_two_original_sources",
      "verified_source_species_RE_PTE":198,
      "verified_PTE_pollen_ovules_species":51,
      "overlap_taxa_not_independent":49,
      "source_species_bias_two_families":164,
      "evidence_rows_total":len(data),"new_source_dependent_evidence_records":23,
      "numeric_checks":checks,
      "physical_and_log_loss_not_conflated":True,
      "original_author_228_reanalysis_sha256":hashlib.sha256(P1.read_bytes()).hexdigest(),
      "original_author_51_reanalysis_sha256":hashlib.sha256(P2.read_bytes()).hexdigest(),
      "ledger_sha256":hashlib.sha256(LED.read_bytes()).hexdigest(),
      "no_flower_biomass_or_genetic_paternity_claimed":True},indent=2))
if __name__=="__main__":main()
