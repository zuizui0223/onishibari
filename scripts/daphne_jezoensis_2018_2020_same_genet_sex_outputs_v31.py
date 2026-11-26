#!/usr/bin/env python3
"""v0.31 | On the same adult genotype, are male and residual female functions both observed?

Original *Daphne jezoensis* Shibata & Kudo 2023 (AJB 110:e16239),
author Zenodo 10.5281/zenodo.8253921, source-pinned in v28.

Adult 'NaturalFruitNo' is fruit count, not viable selfed/outcross seed count.
Offspring with an assigned father are a sampled subset, not total male fitness.
No claim of trade-off: sex morphology and parental success are observational.
No personally identifying location/genotype/ramet source rows are output.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,stringify

OUT=Path("analysis_results_daphne_same_genet_dual_fitness_v31")
PFILES={2018:"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
        2020:"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"}
RNG=np.random.default_rng(20261116)
def fval(x):
    return round(float(x),6) if pd.notna(x) and np.isfinite(float(x)) else None

def inspect_year(d,year):
    adult=d["GenetID_geno_phenotype_1820Nopadult_211109.csv"].copy()
    p=d[PFILES[year]]
    father=p[p.Father.notna()&pd.to_numeric(p.Probability,errors="coerce").ge(.95)].copy()
    assert len(father)==(209 if year==2018 else 261)
    father["genet"]=stringify(father.Father)
    fathers=father.genet.value_counts()
    sc=f"Sex_{year}.g";fl=f"NaturalFlowerNo_{year}.g";fr=f"NaturalFruitNo_{year}.g"
    assert set([sc,fl,fr,"GenetID2"])<=set(adult)
    adult["genet"]=stringify(adult.GenetID2)
    assert not adult.genet.duplicated().any()
    adult["sex"]=adult[sc].astype(str)
    adult["flowers"]=pd.to_numeric(adult[fl],errors="coerce")
    adult["fruits"]=pd.to_numeric(adult[fr],errors="coerce")
    adult["assigned_paternal_seeds"]=adult.genet.map(fathers).fillna(0).astype(int)
    adult["fruit_valid"]=adult.fruits.notna()&adult.flowers.gt(0)&adult.fruits.ge(0)&adult.fruits.le(adult.flowers)
    adult["fruit_rate"]=np.where(adult.fruit_valid,adult.fruits/adult.flowers,np.nan)
    output={"year":year,"n_adult_genets":len(adult),
            "n_assigned_paternal_seeds":len(father),
            "morphs":{}}
    for sex in ("F","H","HF"):
        a=adult[adult.sex.eq(sex)].copy()
        b=a[a.flowers.gt(0)].copy()
        c=b[b.fruit_valid].copy()
        paternal=b[b.assigned_paternal_seeds>0]
        fruiters=c[c.fruits>0]
        dual=c[(c.fruits>0)&(c.assigned_paternal_seeds>0)]
        output["morphs"][sex]={
          "adult_genets_with_exact_label":len(a),
          "adult_genets_with_observed_positive_flower_count":len(b),
          "adult_flowering_with_valid_same_year_fruit_count":len(c),
          "adult_flowering_with_missing_or_invalid_fruit_count":len(b)-len(c),
          "n_genets_with_observed_fruit_positive":len(fruiters),
          "n_genets_with_at_least_one_high_confidence_assigned_paternal_offspring":len(paternal),
          "n_genets_with_both_recorded_fruit_and_assigned_paternal_offspring":len(dual),
          "n_fruits_total_on_valid_genets":int(c.fruits.sum()),
          "n_flowers_total_on_valid_genets":int(c.flowers.sum()),
          "mean_fruits_per_valid_flowering_genet":fval(c.fruits.mean()),
          "pooled_fruit_per_flower":fval(c.fruits.sum()/c.flowers.sum()) if len(c) else None,
          "mean_natural_fruit_rate_genet_equal":fval(c.fruit_rate.mean()),
          "n_assigned_paternal_offspring_from_flowering_genets":int(b.assigned_paternal_seeds.sum()),
          "assigned_paternal_offspring_from_dual_output_genets":int(dual.assigned_paternal_seeds.sum()),
          "n_fruit_positive_genets_without_sampled_assigned_fatherhood":int((c.fruits.gt(0)&c.assigned_paternal_seeds.eq(0)).sum()),
          "n_assigned_father_genets_without_recorded_fruits":int((c.fruits.eq(0)&c.assigned_paternal_seeds.gt(0)).sum()),
          "n_adult_flowering_father_genets_without_valid_fruit_count":int((b.assigned_paternal_seeds.gt(0)&(~b.fruit_valid)).sum()),
          "naive_rank_association_fruits_and_sire_count_in_same_genets":
             fval(spearmanr(c.fruits,c.assigned_paternal_seeds).statistic) if len(c)>5 and c.fruits.nunique()>1 else None,
          "naive_rank_association_fruit_rate_and_sire_count_in_same_genets":
             fval(spearmanr(c.fruit_rate,c.assigned_paternal_seeds).statistic) if len(c)>5 and c.fruit_rate.nunique()>1 else None,
        }
    H=adult[adult.sex.eq("H")&adult.flowers.gt(0)&adult.fruit_valid].copy()
    output["_private_H_ids_for_year_overlap"]=set(H.genet.tolist())
    return output

def main():
    OUT.mkdir(exist_ok=True)
    d=data()
    result=[inspect_year(d,y) for y in (2018,2020)]
    overlap=len(result[0].pop("_private_H_ids_for_year_overlap")&
                result[1].pop("_private_H_ids_for_year_overlap"))
    report={
      "status":"original_Daphne_jezoensis_same_genet_male_paternal_and_residual_fruits_v31",
      "original_article":"10.1002/ajb2.16239",
      "source_record":"10.5281/zenodo.8253921",
      "by_year":result,
      "exact_H_with_both_years_numeric_flower_fruit":overlap,
      "quality_and_inference":{
        "natural_fruits_not_genotyped_viable_offspring":True,
        "sire_count_is_only_assigned_sampled_offspring_not_total_fitness":True,
        "fatherhood_assigned_95percent":True,
        "fruit_no_is_not_hand_self_or_hand_outcross_fertility_assay":True,
        "HF_semantics_not_defined_by_this_table":True,
        "siring_zero_can_reflect_sampling_detection_not_sterility":True,
        "no_causal_male_female_resource_tradeoff_from_correlation":True,
        "original_research_not_new_independent_experiment":True,
        "no_focal_D_pseudomezereum_data":True}}
    assert sum(j["n_assigned_paternal_seeds"] for j in result)==470
    (OUT/"daphne_same_genet_male_siring_female_fruits_v31.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_DUAL_SEX_FITNESS_V31",json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
