#!/usr/bin/env python3
"""v0.48: Source-traced 2005 Daphne laureola five-population stage filters.

Medrano, Alonso and Herrera (2005) Heredity 94:37-43,
DOI 10.1038/sj.hdy.6800550. All numerical data below come from original
printed Tables 1 and 2 (open source www.nature.com/articles/6800550.pdf).
This is an AUDIT of published aggregates, not new raw seeds/individuals.

Separates:
 (a) frequency of *no* stigma pollen load in F vs H plants (2002);
 (b) genetic selfing of *mature seeds* in H and F (2002);
 (c) much lower adult inbreeding coefficient FIS than expected from
     the mating system in 4 of 5 sites (published qualitative result).

Without per-adult or seed individual data there is NO direct lifetime
survival ratio and no proof of causal inbreeding depression. Female t estimates
may exceed 1 (unconstrained estimators); cap at 1 in a predeclared sensitivity
instead of interpreting 1.168 as a literal biological fraction. No original
source counts are independent new publications or P focal Onishibari data.
"""
from __future__ import annotations
from pathlib import Path
import json, math

OUT=Path("analysis_results_daphne_life_stage_filters_v48")
SOURCE="https://www.nature.com/articles/6800550.pdf"
DOI="10.1038/sj.hdy.6800550"
# Table 1 percent female, NoPollenStigmas_F/H = 72 flowers per sex and site.
# Table 2 exact maternal families, genotyped seeds, t_s, t_m.
# t_m=null for sites with 1 polymorphic allozyme; use t_s for that site.
ROWS=[
 ("CV",56.1,950,14.2,40.3,0.0,19,358,.024,.024,6,116,.745,.746,"low"),
 ("RH",34.5,1235,13.7,58.3,6.9,13,249,.026,.020,2,65,1.200,1.168,"low"),
 ("NC",30.1,1615,12.8,84.7,11.1,20,400,.001,.001,None,None,None,None,"low"),
 ("CT",28.1,1640,12.8,73.6,11.1,20,372,.125,None,4,84,.972,None,"low"),
 ("CE",21.4,1575,13.7,31.9,1.4,16,299,.026,None,None,None,None,None,"similar"),
]
SCORING={"CV":57,"RH":55,"NC":133,"CT":96,"CE":145}
ADULT_GENOTYPED={"CV":43,"RH":20,"NC":29,"CT":26,"CE":20}
FIS_DESC={
 "CE":"0.64 approx, near equilibrium, authors exception",
 "NC":"higher than other three sites but notably below estimated equilibrium",
 "CT":"low, below expected equilibrium",
 "RH":"near 0.11, below expected equilibrium",
 "CV":"near 0.11, below expected equilibrium"
}
def val(x):return round(float(x),7)
def expected_F_inbreeding(female_fraction, tH, tF):
    t=female_fraction*tF+(1-female_fraction)*tH
    return (1-t)/(1+t),t

def main():
    OUT.mkdir(exist_ok=True)
    details=[]
    for (site,percent,elev,length,noF,noH,HFfam,HFseed,tHs,tHm,
         Ffam,Fseed,tFs,tFm,fit) in ROWS:
        assert 0<=percent<=100 and 0<=noF<=100 and 0<=noH<=100
        assert 0<=tHs<=1 and (tHm is None or 0<=tHm<=1)
        assert (tFm is None or tFm>=0)
        # H single-locus only when no t_m is available.
        th=tHm if tHm is not None else tHs
        # Female default complete outcross if original source has no
        # F-family estimates; whether actual D. laureola outcrossing is
        # 100% is not proven for unsampled F families.
        tfraw=tFm if tFm is not None else (tFs if tFs is not None else 1.0)
        tfclip=min(1.0,max(0.0,tfraw))
        fexp,tpop=expected_F_inbreeding(percent/100,th,tfclip)
        fexp_simple,_=expected_F_inbreeding(percent/100,0.,1.)
        print_row={
          "site":site,"sex_ratio_percent_female":percent,"adults_scored_for_sex":SCORING[site],
          "adult_plants_genotyped_years_2003_2004":ADULT_GENOTYPED[site],
          "elevation_metres":elev,"H_mean_corolla_length_mm":length,
          "stigma_without_any_pollen_F_percent":noF,
          "stigma_without_any_pollen_H_percent":noH,
          "female_minus_H_no_pollen_percentage_points":val(noF-noH),
          "mating_H_n_maternal_seed_families":HFfam,
          "mating_H_n_genotyped_seeds":HFseed,
          "H_t_s_single_locus_as_reported":tHs,
          "H_t_m_multilocus_as_reported_or_null":tHm,
          "H_selected_t_for_pop_outcross_scenario":th,
          "H_estimated_fraction_mature_genotyped_seeds_selfed":val(1-th),
          "F_n_seed_families":Ffam,
          "F_n_genotyped_seeds":Fseed,
          "F_t_s_as_reported":tFs,
          "F_t_m_as_reported":tFm,
          "F_t_used_for_biologically_clipped_calculation":val(tfclip),
          "F_t_clipped_or_assumed_full":tFm is None and tFs is None or tfraw>1.,
          "pop_mating_fraction_outcross_scenario":val(tpop),
          "pop_mating_equilibrium_Fexp_scenario":val(fexp),
          "pop_mating_equilibrium_Fexp_allHself_allFoutcross":val(fexp_simple),
          "published_adult_FIS_vs_equilibrium_classification":fit,
          "adult_FIS_original_paper_figure_approx_text_only":FIS_DESC[site],
          "not_an_adult_survival_probability":True
        }
        details.append(print_row)
    assert len(details)==5
    assert sum(r["mating_H_n_maternal_seed_families"] for r in details)==88
    assert sum(r["mating_H_n_genotyped_seeds"] for r in details)==1678
    assert sum(ADULT_GENOTYPED.values())==138
    assert sum(SCORING.values())==486
    assert sum(r["F_n_seed_families"] or 0 for r in details)==12
    assert sum(r["F_n_genotyped_seeds"] or 0 for r in details)==265
    assert all(r["female_minus_H_no_pollen_percentage_points"]>=30 for r in details)
    assert all(0<=r["pop_mating_equilibrium_Fexp_scenario"]<=1 for r in details)
    assert 1-max(r["H_estimated_fraction_mature_genotyped_seeds_selfed"] for r in details)<=.00101
    assert min(r["H_estimated_fraction_mature_genotyped_seeds_selfed"] for r in details)>=.875
    assert sum(x["published_adult_FIS_vs_equilibrium_classification"]=="low"
               for x in details)==4
    pct=sum(r["sex_ratio_percent_female"] for r in details)/len(details)
    out={
      "status":"PASS_2005_DAPHNE_LAUREOLA_SOURCE_MATING_STAGE_AUDIT_V48",
      "published_article":"Medrano Alonso Herrera (2005) Heredity 94 37-43",
      "doi":DOI,"printed_source_pdf":SOURCE,
      "source_tables":"Table 1 and Table 2, paper pages 40-41",
      "n_geographic_populations":5,
      "n_adult_plants_scored_morph_sex_not_genotyped":sum(SCORING.values()),
      "n_adult_plants_genotyped_for_FIS_years_2003_2004":sum(ADULT_GENOTYPED.values()),
      "n_H_maternal_seed_families":sum(r["mating_H_n_maternal_seed_families"] for r in details),
      "n_H_genotyped_mature_seeds":sum(r["mating_H_n_genotyped_seeds"] for r in details),
      "n_F_genotyped_mature_seeds":sum(r["F_n_genotyped_seeds"] or 0 for r in details),
      "H_mature_seed_selfing_estimate_range":[
         min(r["H_estimated_fraction_mature_genotyped_seeds_selfed"] for r in details),
         max(r["H_estimated_fraction_mature_genotyped_seeds_selfed"] for r in details)],
      "stigma_no_pollen_F_range_percent":[
         min(r["stigma_without_any_pollen_F_percent"] for r in details),
         max(r["stigma_without_any_pollen_F_percent"] for r in details)],
      "stigma_no_pollen_H_range_percent":[
         min(r["stigma_without_any_pollen_H_percent"] for r in details),
         max(r["stigma_without_any_pollen_H_percent"] for r in details)],
      "n_populations_authors_report_adult_FIS_below_mating_equilibrium":4,
      "mean_female_sex_frequency_percent_unweighted_sites":val(pct),
      "site_original_aggregates_and_derived_equilibrium_scenarios":details,
      "inference_constraints":{
        "printed_site_mating_estimates_not_raw_seed_genotypes":True,
        "2_populations_H_t_m_absent_due_to_one_variable_locus":True,
        "2_populations_F_outcrossing_unmeasured_assumed_1_for_scenario":True,
        "RH_F_reported_t_m_exceeds_1_and_is_clipped_to_1_in_scenario":True,
        "FIS_points_approximate_figure_do_not_use_as_exact_numeric_data":True,
        "reported_four_of_five_late_filter_is_authors_indirect_figure_result":True,
        "equilibrium_Fexp_assumes_stable_mating_sex_frequency_and_no_population_structure":True,
        "no_causal_late_inbreeding_depression_directly_tracked_to_adults":True,
        "no_original_sitewise_adult_genotype_data_acquired":True,
        "no_maternal_paternal_fitness_equality_claim":True,
        "same_published_study_five_populations_not_five_independent_publications":True,
        "no_direct_Onishibari_reproduction_experiment":True
      }}
    path=OUT/"daphne_laureola_2005_5pop_mature_seed_to_adult_filter_v48.json"
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_LAUREOLA_LIFE_STAGE_V48",json.dumps({
        "status":out["status"],"n_sites":5,"H_seed_families":out["n_H_maternal_seed_families"],
        "H_genotyped_seeds":out["n_H_genotyped_mature_seeds"],
        "H_selfed_mature_seed_fraction_range":out["H_mature_seed_selfing_estimate_range"],
        "F_no_pollen_stigmas_percent_range":out["stigma_no_pollen_F_range_percent"],
        "H_no_pollen_stigmas_percent_range":out["stigma_no_pollen_H_range_percent"],
        "adult_FIS_lower_than_scenario_in_four_of_five_published":True,
        "equilibrium_scenarios":[{"population":r["site"],
              "t":r["pop_mating_fraction_outcross_scenario"],
              "Fexp":r["pop_mating_equilibrium_Fexp_scenario"],
              "source_estimate_F_T_clipped_or_assumed":r["F_t_clipped_or_assumed_full"]}
             for r in details],
        "not_direct_survival_or_our_mating_experiment":True
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
