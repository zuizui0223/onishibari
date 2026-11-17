#!/usr/bin/env python3
"""Numerical/provenance assertions on the 2010 Daphne historical barcode archive v43.

No claims of missing original sequence data based only on ESearch zero hits.
"""
from pathlib import Path
import json
OUT=Path("analysis_results_daphne_2010_barcode_archive_provenance_v43")
FILE=OUT/"historical_daphne_2010_barcode_record_provenance_v43.json"
REPORT=Path("results/2026-11-17_daphne_2010_barcode_sequence_recovery_audit_v43.md")
def main():
    a=json.loads(FILE.read_text(encoding="utf8"))
    assert a["status"]=="GenBank_first_record_provenance_original_Daphne_2010_v43"
    assert a["barcoding_presentation_doi"]=="10.11519/jfsc.121.0.332.0"
    assert a["sequence_record_count_inspected"]==62
    assert a["records_with_any_team_first_reference_name_match"]==0
    d={r["accession"]:r for r in a["record_metadata"]}
    assert len(d)==62
    assert d["GU736360"]["organism"]=="Daphne tangutica"
    assert d["GU736359"]["organism"]=="Daphne mezereum"
    assert len([x for x in d if x.startswith("GU")])==2
    assert sum(v["source_search_label"]=="Daphne_ITS2_2010" for v in d.values())==60
    assert sum(v["source_search_label"]=="Daphne_plastid_spacer_2010" for v in d.values())==2
    assert len(a["team_matching_accessions"])==0
    tallies=a["record_first_reference_title_tally"]
    assert len(tallies)==2
    assert next(v for k,v in tallies.items() if "Back-and-forth hermaphroditism" in k)==60
    assert next(v for k,v in tallies.items() if "Stellera chamaejasme" in k)==2
    assert "Zhang" in d["GU736360"]["first_reference_authors"]
    assert "Stellera chamaejasme" in d["GU736360"]["first_reference_title"]
    assert d["GQ167547"]["organism"]=="Daphne rodriguezii"
    assert d["GQ167544"]["organism"]=="Daphne oleoides"
    for k,v in d.items():
        if v["source_search_label"]=="Daphne_ITS2_2010":
            assert "Alonso" in v["first_reference_authors"]
            assert "Hermaphroditism" in v["first_reference_title"].title()
        assert not v["team_author_match_in_first_reference"]
    for name in ("original_seven_barcode_markers_six_plastid_one_multicopy_nuclear",
                 "first_reference_metadata_not_exhaustive_source_provenance",
                 "no_original_2010_sequence_accession_claim_without_voucher_crosswalk",
                 "search_subset_not_all_historical_NCBI_records",
                 "no_focal_Onishibari_new_genetic_measurement"):
        assert a["strict_guard"][name]
    assert REPORT.exists()
    txt=REPORT.read_text(encoding="utf8")
    assert all(x in txt for x in ("62件","GU736360","GQ167547","2010","ILS","独立"))
    print(json.dumps({
     "status":"PASS_DAPHNE_HISTORICAL_BARCODE_ORIGIN_V43",
     "inspected_accessions":62,
     "original_Kawahara_project_voucher_matches_confirmed":0,
     "not_a_negative_proof_of_any_original_accession":True,
     "not_an_onishibari_original_sequence_result":True},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
