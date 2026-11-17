#!/usr/bin/env python3
"""v0.40 | Validate provenance + taxonomic concept + functional-sex evidence for Daphne.

DO NOT turn a morphological dioecy/gynodioecy claim into direct evidence that
Onishibari residual pistils produce viable offspring. Likewise, historic
'D. kamtchatica var. jezoensis' denotes current D. jezoensis, not the modern
D. kamtschatica type variety. Northern Honshu 'var. koreana' plants are not
proven to match the Korean lineage now accepted at POWO as D. kamtschatica
var. koreana. Sources, geographic concepts and assay levels are separate.

This checks a hand-reviewed reference ledger. No new field/sequence data.
"""
from __future__ import annotations
import csv,io,json,re
from pathlib import Path
from collections import Counter,defaultdict

HERE=Path(__file__).resolve().parents[1]
INPUT=HERE/"data/derived/daphne_sex_taxon_concept_provenance_v40.csv"
OUT=HERE/"analysis_results_daphne_sex_taxon_v40"
COLS={"claim_id","source_year","name_as_printed","resolved_taxon_concept",
      "geographic_scope","axis","evidence_level","observed_or_described_state",
      "paraphrased_claim","has_direct_functional_assay_for_claim","source_url",
      "interpretation_guard"}
P="Daphne pseudomezereum var. pseudomezereum"
J="Daphne jezoensis"
K="Daphne kamtschatica var. kamtschatica"
C="Daphne kamtschatica var. koreana"
H="Daphne sp. Honshu koreana sensu Nitta2020 [unresolved]"
I="Daphne kiusiana"
L="Daphne laureola"

def read():
    txt=INPUT.read_text(encoding="utf-8-sig")
    with io.StringIO(txt,newline="") as h:
        r=csv.DictReader(h)
        assert r.fieldnames and set(r.fieldnames)==COLS, ("Unexpected columns",r.fieldnames)
        d=list(r)
    assert len(d)==32 and len({v["claim_id"] for v in d})==32
    for v in d:
        assert all(str(t).strip() for t in v.values()),v["claim_id"]
        assert v["source_url"].startswith("https://"),v["claim_id"]
        v["has_direct_functional_assay_for_claim"]=v["has_direct_functional_assay_for_claim"].lower()
        assert v["has_direct_functional_assay_for_claim"] in ("true","false"),v["claim_id"]
        assert "PLACEHOLDER" not in v["source_url"] and "unknown DOI" not in v["source_url"]
    return d

def byid(d):
    return {x["claim_id"]:x for x in d}

def checks(d):
    a=byid(d)
    assert set(a)=={
      "P01","P02","P03","P04","P05","P06","P07",
      "J01","J02","J03","J04","J05","J06","J07","J08","J09",
      "K01","C01","C02","C03","C04","I01","I02","I03",
      "L01","L02","L03","M01","M02","M03","M04","M05"}
    # Focal contradiction is explicitly sourced, neither claims functional
    # fecundity of male-appearing flowers directly.
    assert a["P02"]["observed_or_described_state"]=="gynodioecy"
    assert a["P03"]["observed_or_described_state"]=="dioecy"
    assert a["P04"]["observed_or_described_state"]=="ambiguous_dioecy"
    assert a["P02"]["source_year"]=="2010"
    assert a["P03"]["source_year"]=="2020"
    assert all(a[k]["resolved_taxon_concept"]==P for k in ("P01","P02","P03","P04","P05","P06","P07"))
    assert all(a[k]["has_direct_functional_assay_for_claim"]=="false" for k in ("P01","P02","P03","P04","P05","P06"))
    assert a["P06"]["axis"]=="self_compatibility"
    assert a["P06"]["observed_or_described_state"]=="low_self_compatibility_reported_unquantified"
    assert a["P06"]["has_direct_functional_assay_for_claim"]=="false"
    assert a["M01"]["resolved_taxon_concept"]=="Daphne miyabeana"
    assert a["M02"]["observed_or_described_state"]=="self_compatible_spontaneous_selfing_low_vs_outcross"
    assert a["M02"]["has_direct_functional_assay_for_claim"]=="true"
    assert a["M03"]["axis"]=="pollination_function"
    assert a["M04"]["observed_or_described_state"]=="gynodioecy"
    assert a["M05"]["observed_or_described_state"]=="dioecy"
    assert all(a[k]["resolved_taxon_concept"]=="Daphne miyabeana" for k in ("M01","M02","M03","M04","M05"))
    # One historically named genotype must not count as independent
    # D. kamtschatica evidence for present-day D. kamtschatica.
    assert a["J01"]["name_as_printed"]=="Daphne kamtchatica var. jezoensis"
    assert a["J01"]["resolved_taxon_concept"]==J
    assert a["J03"]["resolved_taxon_concept"]==J
    assert all(a[k]["resolved_taxon_concept"]==J for k in ("J07","J08","J09"))
    assert all(a[k]["has_direct_functional_assay_for_claim"]=="true" for k in ("J07","J08","J09","L03"))
    assert a["J07"]["name_as_printed"]=="Daphne kamtchatica var. jezoensis"
    assert a["J08"]["observed_or_described_state"]=="same_flower_self_zero_other_flower_same_plant_fruit_positive"
    assert a["J09"]["observed_or_described_state"]=="H_flower_number_predicts_paternity_ALREADY_REPORTED"
    assert a["L03"]["observed_or_described_state"]=="emasculation_geitonogamous_self_fruit_33_7_vs_18_2_percent"
    assert a["K01"]["resolved_taxon_concept"]==K
    assert a["K01"]["observed_or_described_state"]=="hermaphrodite"
    assert a["K01"]["has_direct_functional_assay_for_claim"]=="false"
    assert a["I02"]["resolved_taxon_concept"]==I
    assert a["I02"]["has_direct_functional_assay_for_claim"]=="true"
    assert a["I02"]["observed_or_described_state"]=="self_compatible_high_fruiting"
    assert a["J03"]["observed_or_described_state"]=="self_compatible_weak_assurance"
    # Honor Korean phylogeny/synonyms vs old Honshu unresolved var label.
    assert a["C01"]["resolved_taxon_concept"]==H
    assert a["C01"]["observed_or_described_state"]=="unknown"
    assert a["C02"]["resolved_taxon_concept"]==C
    assert a["C03"]["resolved_taxon_concept"]==C
    assert a["C04"]["resolved_taxon_concept"]==H
    assert a["C04"]["observed_or_described_state"]=="summergreen_Japanese_koreana_not_genetically_confirmed_Korean_C"
    assert a["P07"]["observed_or_described_state"]=="P_and_J_southern_northern_races_of_one_species_unpublished"
    assert a["L01"]["source_url"]=="https://doi.org/10.2307/2657083"
    assert a["L01"]["has_direct_functional_assay_for_claim"]=="true"
    assert a["L01"]["observed_or_described_state"]=="vector_required_self_compatible"
    assert a["P04"]["source_year"]=="undated"
    assert a["P05"]["source_year"]=="undated"
    assert a["C03"]["source_year"]=="undated"
    # Verify no false line of ancestry from mere 3 tips.
    assert "dioecy"!=a["J01"]["observed_or_described_state"]
    # Every named external taxon used in a planned staged comparison is distinct.
    assert len({P,J,K,I,C,H,"Daphne miyabeana"})==7

def summary(d):
    a=byid(d)
    group=defaultdict(list)
    axis=defaultdict(list)
    for x in d:
        group[x["resolved_taxon_concept"]].append(x["claim_id"])
        axis[x["axis"]].append(x["claim_id"])
    states={k:sorted(set(x["observed_or_described_state"] for x in d if
          x["resolved_taxon_concept"]==k and x["axis"]=="sex_structure"))
          for k in sorted(group)}
    focal_direct=[x["claim_id"] for x in d if x["resolved_taxon_concept"]==P and
       x["has_direct_functional_assay_for_claim"]=="true"]
    assert focal_direct==[]
    result={
      "status":"PASS_DAPHNE_SEX_CONCEPT_PROVENANCE_V40",
      "source_count":len(d),
      "unique_assertion_ids":len(a),
      "taxon_concepts_present":sorted(group),
      "claim_counts_by_taxon":{k:len(v) for k,v in sorted(group.items())},
      "claim_ids_by_axis":{k:sorted(v) for k,v in sorted(axis.items())},
      "sex_structure_state_claims_by_taxon":states,
      "focal_onishibari_claim_discordance":{
         "morphological_and_reported_system":"UNRESOLVED_CONFLICT",
         "ffpri_2010_claim":"gynodioecy",
         "jjb_2020_table_claim":"dioecy",
         "botanical_garden_statement":"dioecy_unclear",
         "focal_direct_cross_viable_offspring_assays_in_ledger":0,
         "focal_residual_pistil_ability_to_produce_viable_seed":"UNKNOWN",
         "focal_self_compatibility":"LOW_REPORTED_2011_KAKEN_UNQUANTIFIED_NOT_CONFIRMED_PHYSIOLOGICAL_SI",
         "focal_hermaphrodite_functionality":"UNKNOWN",
         "NOT_a_reclassification_from_this_meta_analysis":True
       },
      "nonfocal_functional_contrast":{
        "D_jezoensis":"male_biased_H_can_be_self_compatible_but_low_seeding_and_inbreeding_depression",
        "D_kiusiana":"H_hand_cross_has_substantial_fruit_in_reported_populations",
        "same_morphological_gynodioecy_does_not_force_same_function":True,
        "D_miyabeana_2011_gyn_vs_2022_dioecy_claim_conflict":True,
        "historic_var_jezoensis_alias_is_modern_jezoensis":True,
        "modern_D_kamtschatica_sensu_stricto_not_empirically_selfing_assayed_here":True,
        "1989_jezoensis_emasculation_pollen_clogging_previously_published":True,
        "2008_jezoensis_same_vs_other_self_source_treatments_non_equivalent":True,
        "2019_jezoensis_H_flower_number_paternity_association_already_reported":True,
        "2001_laureola_H_geitonogamous_self_emasculation_result_with_only_five_H_plants":True
      },
      "taxonomy_guard":{
        "Honshu_var_koreana":"UNRESOLVED_OTU_NOT_AUTOMATICALLY_KOREAN_C",
        "Korean_var_koreana":"Daphne_kamtschatica_var_koreana_acceptance_by_POWO",
        "Korean_cp_and_nuclear_marker_not_same_individual_in_existing_source":True,
        "no_three_tip_linear_ancestral_order":True,
        "P_jezoensis_same_species_claim_only_unpublished_phylogeny":True,
        "Japanese_2020_koreana_montane_phenology_not_Korean_2025_taxon_assignment":True,
      },
      "next_data_gate":"Focal H-morph hand-cross and hand-self versus verified F positive control; genetically viable seed and actual dry cost, with vouchers",
      "not_an_original_field_experiment":True,
      "not_multilocus_nuclear_phylogeny_of_focal_P_C_K":True,
      "no_new_genetic_mating_or_fitness_measurements":True
    }
    return result

def main():
    d=read()
    checks(d)
    r=summary(d)
    OUT.mkdir(exist_ok=True)
    (OUT/"daphne_source_sex_claim_contradictions_v40.json").write_text(
      json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("DAPHNE_TAXON_SEX_PROVENANCE_V40",json.dumps(r,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
