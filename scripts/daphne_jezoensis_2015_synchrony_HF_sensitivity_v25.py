#!/usr/bin/env python3
"""v0.25 robust original 2015 Daphne phenological asynchrony sensitivity.

HF plants are source-classified ambiguously. Compare H-only with H+all HF
counted hypothetically as donor candidate flowers (NOT functional verified pollen).
Also compare sex-temporal centers with absolute open flower counts vs equal
independent-ramet normalized open fractions. No actual compatible mates known.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from daphne_jezoensis_2015_sex_synchrony_v24 import (
    load, DAYS, LABELS, ZIP_SHA, CSV_SHA, time_stats
)

OUT=Path("analysis_results_daphne_2015_synchrony_sensitivity_v25")
RNG=np.random.default_rng(20261116)
BOOT=10000

def near(x):
    return round(float(x),7)

def ci(v):
    return [near(x) for x in np.quantile(np.asarray(v,dtype=float),[.025,.975])]

def mean_date_distribution(mat):
    denom=np.trapezoid(mat,DAYS,axis=-1)
    cen=np.trapezoid(mat*DAYS,DAYS,axis=-1)/denom
    return cen

def summarize():
    d=load()
    d["normalized_open"]=d.open/d.total
    groups={}
    for sex,n_expected in [("F",76),("H",90),("HF",15)]:
        a=d[d.Sex==sex]
        raw=a.pivot(index="RametID",columns="day",values="open")[list(DAYS)].to_numpy(float)
        normalized=a.pivot(index="RametID",columns="day",values="normalized_open")[list(DAYS)].to_numpy(float)
        assert raw.shape==(n_expected,5)
        assert normalized.shape==(n_expected,5)
        groups[sex]={"raw":raw,"frac":normalized}
    f=groups["F"]["raw"];h=groups["H"]["raw"];hf=groups["HF"]["raw"]
    fn=groups["F"]["frac"];hn=groups["H"]["frac"];hfn=groups["HF"]["frac"]
    source={"F":f.sum(axis=0),"H":h.sum(axis=0),"HF":hf.sum(axis=0)}
    donor_without=source["H"]
    donor_with=source["H"]+source["HF"]
    fsum=source["F"]
    # Each case is a source-backed display-opportunity scenario, not pollen availability.
    point={
      "source_counts_H_only":time_stats(fsum,donor_without),
      "source_counts_H_plus_HF_hypothetical_donor":time_stats(fsum,donor_with),
      "equal_plant_normalized_F_vs_H":time_stats(fn.mean(axis=0),hn.mean(axis=0)),
      "equal_plant_normalized_F_vs_H_plus_HF_hypothetical_donor":time_stats(
        fn.mean(axis=0),np.concatenate([hn,hfn],axis=0).mean(axis=0))
    }
    bydate=[]
    for i,day in enumerate(DAYS):
        bydate.append({
          "date":LABELS[i],
          "F_open_sample_count":int(fsum[i]),
          "H_open_sample_count":int(source["H"][i]),
          "HF_ambiguous_open_sample_count":int(source["HF"][i]),
          "H_plus_HF_total_if_all_pollen_donors":int(donor_with[i]),
          "H_only_div_F":near(donor_without[i]/fsum[i]),
          "H_plus_HF_candidate_div_F":near(donor_with[i]/fsum[i]),
          "n_HF_plant_open":int((hf[:,i]>0).sum()),
          "sex_balanced_H_open_per_individual_over_F":near(
             h[:,i].mean()/f[:,i].mean()),
          "sex_balanced_H_plus_HF_mean_per_individual_over_F":near(
             np.concatenate([h[:,i],hf[:,i]]).mean()/f[:,i].mean())
        })
    fi=RNG.integers(0,len(f),size=(BOOT,len(f)))
    hi=RNG.integers(0,len(h),size=(BOOT,len(h)))
    hfi=RNG.integers(0,len(hf),size=(BOOT,len(hf)))
    bF=f[fi].sum(axis=1)
    bH=h[hi].sum(axis=1)
    bHF=hf[hfi].sum(axis=1)
    bFfrac=fn[fi].mean(axis=1)
    bHfrac=hn[hi].mean(axis=1)
    bHFfrac=hfn[hfi].mean(axis=1)
    bHall=bH+bHF
    bHfrac_all=(bHfrac*len(h)+bHFfrac*len(hf))/(len(h)+len(hf))
    boot={
       "center_diff_F_minus_H_raw":ci(mean_date_distribution(bF)-mean_date_distribution(bH)),
       "center_diff_F_minus_HplusHF_candidate_raw":ci(mean_date_distribution(bF)-mean_date_distribution(bHall)),
       "center_diff_F_minus_H_equal_plant":ci(mean_date_distribution(bFfrac)-mean_date_distribution(bHfrac)),
       "center_diff_F_minus_HplusHF_candidate_equal_plant":ci(mean_date_distribution(bFfrac)-mean_date_distribution(bHfrac_all)),
       "H_only_div_F_final_day":ci(bH[:,4]/bF[:,4]),
       "HplusHF_candidate_div_F_final_day":ci(bHall[:,4]/bF[:,4]),
       "H_only_div_F_peak_Apr27":ci(bH[:,2]/bF[:,2]),
       "HplusHF_candidate_div_F_peak_Apr27":ci(bHall[:,2]/bF[:,2])
    }
    # Verify original source and whether the primary counterfactual late donor pool
    # includes enough HF flowers to change the H-only narrative.
    assert source["F"][4]==654 and source["H"][4]==47
    assert int((hf[:,4]>0).sum())==12
    return {
       "status":"source_v25_donor_sex_code_and_plant_weight_sensitivity",
       "source":"10.5281/zenodo.8253921",
       "original_DOI":"10.1002/ajb2.16239",
       "original_zip_sha256":ZIP_SHA,"original_CSV_sha256":CSV_SHA,
       "sex_group_size":{"F":76,"H":90,"HF_ambiguous":15},
       "by_date":bydate,
       "timing_metrics":point,
       "independent_plant_bootstrap_95%":boot,
       "validation_notes":{
         "HF_is_not_verified_pollen_donor":True,
         "including_all_HF_as_donor_is_counterfactual_only":True,
         "no_neighborhood_specific_compatibility":True,
         "female_open_flowers_not_verified_receptive":True,
         "male_open_flowers_not_male_fitness_or_viable_pollen":True,
         "same_published_2015_field_census_as_v23_v24":True,
         "no_Daphne_pseudomezereum_measurements":True,
         "do_not_infer_late_season_pollen_limitation_from_flower_abundance":True
       }}
def main():
    OUT.mkdir(exist_ok=True)
    d=summarize()
    path=OUT/"daphne_2015_synchrony_HF_ambiguous_sensitivity_v25.json"
    path.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_SEX_SYNCHRONY_ROBUSTNESS_V25",json.dumps(d,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
