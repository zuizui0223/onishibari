#!/usr/bin/env python3
"""v0.34: source-original within-genet 2018->2020 sex-code and fruit/paternity trajectories.

Author original: Shibata & Kudo (2023), AJB e16239, Zenodo 8253921.
2016 publication already documents ~15% sex phenotype switches across 6 years;
this is a different 2-year interval in one source plot, NOT an independent
replication and NOT proof of functional physiological sex change.
HF category meaning is not source-verified; retain exact code transitions.
Output aggregate only; original genet IDs/coordinates/genotypes never released.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data, stringify

OUT=Path("analysis_results_daphne_2018_2020_sex_states_v34")
PFILES={2018:"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
        2020:"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"}
def n(x):return int(x)
def foo(x):return round(float(x),6) if pd.notna(x) and np.isfinite(x) else None
def groupcode(x):
    k=str(x)
    if k in ("H","HF","F"):return k
    if k=="nan" or k in ("dead","NF","0","-"):return "unclassified_or_nonflower"
    return "mixed_or_other"
def inspect(d):
    adult=d["GenetID_geno_phenotype_1820Nopadult_211109.csv"].copy()
    assert len(adult)==492 and adult.GenetID2.nunique()==492
    adult["key"]=stringify(adult.GenetID2)
    assert not adult.key.duplicated().any()
    for year in (2018,2020):
        adult[f"sex_{year}"]=adult[f"Sex_{year}.g"].map(groupcode)
        adult[f"raw_{year}"]=adult[f"Sex_{year}.g"].fillna("NA").astype(str)
        adult[f"flo_{year}"]=pd.to_numeric(adult[f"NaturalFlowerNo_{year}.g"],errors="coerce")
        adult[f"fruit_{year}"]=pd.to_numeric(adult[f"NaturalFruitNo_{year}.g"],errors="coerce")
        adult[f"valid_{year}"]=adult[f"flo_{year}"].gt(0)&adult[f"fruit_{year}"].notna()&adult[f"fruit_{year}"].ge(0)&adult[f"fruit_{year}"].le(adult[f"flo_{year}"])
        papa=d[PFILES[year]]
        assigned=papa[papa.Father.notna()&pd.to_numeric(papa.Probability,errors="coerce").ge(.95)]
        assert len(assigned)==(209 if year==2018 else 261)
        adult[f"paternal_{year}"]=adult.key.map(stringify(assigned.Father).value_counts()).fillna(0).astype(int)
    panel=adult[adult["valid_2018"]&adult["valid_2020"]].copy()
    trans=panel.groupby(["sex_2018","sex_2020"],dropna=False).size().reset_index(name="count")
    labels=sorted(set(adult.sex_2018)|set(adult.sex_2020))
    matrix=[]
    for a in labels:
        for b in labels:
            m=panel[(panel.sex_2018==a)&(panel.sex_2020==b)]
            if len(m):
                matrix.append({"code_2018":a,"code_2020":b,"n_genets":len(m),
                 "fruited_2018":int(m.fruit_2018.gt(0).sum()),
                 "fruited_2020":int(m.fruit_2020.gt(0).sum()),
                 "paternal_2018_any":int(m.paternal_2018.gt(0).sum()),
                 "paternal_2020_any":int(m.paternal_2020.gt(0).sum()),
                 "mean_flowers_2018":foo(m.flo_2018.mean()),
                 "mean_flowers_2020":foo(m.flo_2020.mean()),
                 "mean_fruit_rate_2018":foo((m.fruit_2018/m.flo_2018).mean()),
                 "mean_fruit_rate_2020":foo((m.fruit_2020/m.flo_2020).mean())})
    hpool=panel[panel.sex_2018.isin(["H","HF"])&panel.sex_2020.isin(["H","HF"])].copy()
    exactH=panel[(panel.sex_2018=="H")&(panel.sex_2020=="H")].copy()
    switch=panel[(panel.sex_2018!=panel.sex_2020)&panel.sex_2018.isin(["H","HF","F"])&
                  panel.sex_2020.isin(["H","HF","F"])]
    out={
        "status":"original_2018_2020_genet_sex_code_longitudinal_v34",
        "original_doi":"10.1002/ajb2.16239",
        "source_doi":"10.5281/zenodo.8253921",
        "n_source_adult_genets":len(adult),
        "n_genets_with_valid_positive_flowers_and_fruits_both_years":len(panel),
        "n_both_years_male_function_codes_H_or_HF":len(hpool),
        "n_stable_exact_H":len(exactH),
        "n_exact_F_H_HF_label_switch_in_both_year_flowering_genets":len(switch),
        "all_source_sex_code_counts":{
          str(year):{str(k):int(v) for k,v in adult[f"sex_{year}"].value_counts(dropna=False).items()}
          for year in (2018,2020)},
        "flowering_both_year_transitions":matrix,
        "source_male_function_subset":{
          "sex_code_H_HF_transition_counts":{
            f"{a}->{b}":int(len(hpool[(hpool.sex_2018==a)&(hpool.sex_2020==b)]))
            for a in ("H","HF") for b in ("H","HF")},
          "year_2018_fruit_positive":int(hpool.fruit_2018.gt(0).sum()),
          "year_2020_fruit_positive":int(hpool.fruit_2020.gt(0).sum()),
          "both_years_fruit_positive":int((hpool.fruit_2018.gt(0)&hpool.fruit_2020.gt(0)).sum()),
          "2018_fruit_zero_to_2020_positive":int((hpool.fruit_2018.eq(0)&hpool.fruit_2020.gt(0)).sum()),
          "2018_fruit_positive_to_2020_zero":int((hpool.fruit_2018.gt(0)&hpool.fruit_2020.eq(0)).sum()),
          "both_years_no_fruits":int((hpool.fruit_2018.eq(0)&hpool.fruit_2020.eq(0)).sum()),
          "with_paternal_assigned_2018":int(hpool.paternal_2018.gt(0).sum()),
          "with_paternal_assigned_2020":int(hpool.paternal_2020.gt(0).sum()),
          "with_paternal_assigned_both_years":int((hpool.paternal_2018.gt(0)&hpool.paternal_2020.gt(0)).sum())},
        "stable_H_year_2018_fruit_and_year_2020_natural_flower":{
          "n_stable_H":len(exactH),
          "n_fruiting_2018":int(exactH.fruit_2018.gt(0).sum()),
          "n_nonfruiting_2018":int(exactH.fruit_2018.eq(0).sum()),
          "2020_flo_mean_if_fruited_2018":foo(exactH.loc[exactH.fruit_2018.gt(0),"flo_2020"].mean()),
          "2020_flo_mean_if_no_fruit_2018":foo(exactH.loc[exactH.fruit_2018.eq(0),"flo_2020"].mean()),
          "baseline_2018_flo_mean_if_fruited_2018":foo(exactH.loc[exactH.fruit_2018.gt(0),"flo_2018"].mean()),
          "baseline_2018_flo_mean_if_no_fruit_2018":foo(exactH.loc[exactH.fruit_2018.eq(0),"flo_2018"].mean())},
        "quality":{
          "two_timepoints_no_2019_intermediate_measurement":True,
          "HF_is_not_verified_as_exogenous_female_function_treatment":True,
          "observed_2018_to_2020_label_change_not_proof_of_biological_sex_reversal":True,
          "conditioning_on_flowering_both_years_is_selection":True,
          "nonassigned_siring_not_necessarily_zero_true_paternity":True,
          "individual_raw_identifiers_and_genotypes_not_in_report":True,
          "same_original_study_as_v30_v33":True,
          "no_focal_D_pseudomezereum_data":True}
    }
    return out
def main():
    OUT.mkdir(exist_ok=True)
    result=inspect(data())
    (OUT/"daphne_2018_2020_genet_sex_code_transition_v34.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_SEX_CODE_TRANSITION_V34",json.dumps(result,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
