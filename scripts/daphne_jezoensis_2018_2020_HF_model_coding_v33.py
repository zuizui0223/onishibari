#!/usr/bin/env python3
"""v0.33: fatherhood-floral scaling after adding all HF-coded male-function candidates.

Retests H-only father model using H+HF adult genotype risk set. HF-coded
plants have demonstrated paternal SSR offspring, but the meaning of the HF
sex category is not causally independent of fruiting status. This is a
source-coding sensitivity, not randomized evidence that HF causes paternity.

Same source: Shibata & Kudo 2023, DOI 10.1002/ajb2.16239, Zenodo 8253921.
Output aggregate numbers only (no adult genotype ID, coordinates, or raw rows).
"""
import json
from pathlib import Path
import numpy as np,pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import (
  build,estimate,boot
)
OUT=Path("analysis_results_daphne_HF_coding_father_model_v33")
LAMBDAS=(5.,10.,20.)

def f(x):return round(float(x),7)
def fruit_counts(ad,year):
    sex=f"Sex_{year}.g";flow=f"NaturalFlowerNo_{year}.g";fr=f"NaturalFruitNo_{year}.g"
    x=ad[[sex,flow,fr]].copy()
    x["flowers"]=pd.to_numeric(x[flow],errors="coerce")
    x["fruits"]=pd.to_numeric(x[fr],errors="coerce")
    x=x[x.flowers.gt(0)&x.fruits.notna()&x.fruits.ge(0)&x.fruits.le(x.flowers)]
    out={}
    for label,codes in [("H_exact",["H"]),("HF_exact",["HF"]),
                        ("H_plus_HF",["H","HF"]),
                        ("F_exact",["F"])]:
        z=x[x[sex].isin(codes)]
        out[label]={
          "genets_with_valid_natural_flower_fruit":len(z),
          "fruit_positive_genets":int(z.fruits.gt(0).sum()),
          "fruit_positive_genet_fraction":f(z.fruits.gt(0).mean()),
          "flower_sum":int(z.flowers.sum()),
          "fruit_sum":int(z.fruits.sum()),
          "pooled_fruit_rate":f(z.fruits.sum()/z.flowers.sum())}
    assert out["H_plus_HF"]["genets_with_valid_natural_flower_fruit"]==(108 if year==2018 else 79)
    assert out["H_plus_HF"]["fruit_sum"]==(67 if year==2018 else 96)
    return out

def main():
    OUT.mkdir(exist_ok=True)
    files=data()
    ad=files["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    years=[]
    for year in (2018,2020):
        sourceH=build(year,files)
        sourceHF=build(year,files,include_HF=True)
        expected={2018:{"H":180,"HF":183,"n_HF":3},
                  2020:{"H":206,"HF":214,"n_HF":8}}[year]
        assert int(sourceH["counts"].sum())==expected["H"]
        assert int(sourceHF["counts"].sum())==expected["HF"]
        assert sourceHF["counts"].shape[1]-sourceH["counts"].shape[1]==expected["n_HF"]
        point=[]
        for lam in LAMBDAS:
            only=estimate(sourceH,lam,extra=True)
            incl=estimate(sourceHF,lam,extra=True)
            uncertainty=boot(sourceHF,lam)
            point.append({
              "lambda_m":int(lam),
              "strict_H_beta":only["beta_hat"],
              "strict_H_assigned_child_n":int(sourceH["counts"].sum()),
              "H_plus_HF_beta":incl["beta_hat"],
              "H_plus_HF_cluster_CI95":uncertainty["beta_mother_boot_ci95"],
              "H_plus_HF_mother_holdout_gain_per_child":incl["mother_leave_one_out_gain_per_child"],
              "n_mothers_H_only":len(sourceH["counts"]),
              "n_mothers_H_plus_HF":len(sourceHF["counts"]),
              "n_assigned_offspring_H_plus_HF":int(sourceHF["counts"].sum()),
              "n_eligible_H_plus_HF_fathers":sourceHF["counts"].shape[1],
              "difference_beta_not_causal_HF_effect":f(incl["beta_hat"]-only["beta_hat"])
            })
        years.append({
          "year":year,
          "paternity_source_data_counts_strict_H":sourceH["qc"],
          "paternity_source_data_counts_H_plus_HF":sourceHF["qc"],
          "adult_fruit_rate_sex_code_sensitivity":fruit_counts(ad,year),
          "flower_to_father_choice_model":point
        })
    out={"status":"source_verified_Daphne_H_plus_HF_father_choice_v33",
      "paper_doi":"10.1002/ajb2.16239",
      "source_doi":"10.5281/zenodo.8253921",
      "yearly":years,
      "limits":{
        "HF_sex_code_may_be_assigned_using_prior_fruiting_outcome":True,
        "HF_male_function_verified_for_subset_only":True,
        "comparing_H_vs_HplusHF_conditions_on_different_sire_sets":True,
        "natural_fruit_numbers_are_not_genetically_verified_viable_maternal_offspring":True,
        "assigned_95percent_parentage_not_lifetime_sire_output":True,
        "male_flower_number_not_actual_construction_dry_mass":True,
        "same_original_2023_study_no_independent_replication":True,
        "not_Daphne_pseudomezereum_experiment":True}}
    (OUT/"daphne_2018_2020_HF_inclusion_father_choice_v33.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_HF_CODE_FATHER_FIT_V33",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
