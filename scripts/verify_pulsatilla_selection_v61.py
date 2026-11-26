#!/usr/bin/env python3
"""Pin all v61 original source identity, missingness and uncertainty numerics."""
import argparse,json,math
from pathlib import Path
def eq(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:eq(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):eq(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-8),(path,a,b)
    else:
        assert a==b,(path,a,b)
def contrast_sensitivity(v,delta,q_all,q_intact,lost_all_mean):
    groups=v["cohort_flow_by_arm"]
    comp=v["tipping_by_delta_original_assigned"][str(delta)]
    A,I=groups["all"],groups["intact"]
    sumA=comp["typed_sum_all"]+A["known_mature_outcome_untyped_harvest_seeds_sum"]*(1-delta*q_all)
    sumI=comp["typed_sum_intact"]+I["known_mature_outcome_untyped_harvest_seeds_sum"]*(1-delta*q_intact)
    sumA+=A["no_matching_mature_outcome_in_source"]*lost_all_mean
    return ((I["original_assigned_single_flower_bisexual"] /
             A["original_assigned_single_flower_bisexual"] * sumA - sumI) /
             I["no_matching_mature_outcome_in_source"])
if __name__=="__main__":
    q=argparse.ArgumentParser()
    q.add_argument("--actual",type=Path,required=True)
    q.add_argument("--expected",type=Path,default=Path("data/2026-11-24_pulsatilla_selection_tipping_v61.json"))
    args=q.parse_args()
    a=json.loads(args.actual.read_text());b=json.loads(args.expected.read_text())
    eq(a,b)
    assert a["source_single_flower_bisexual_original_assigned"]==81
    assert a["source_genotyped_bisexual_original_2024"]==60
    assert a["source_genotyped_bisexual_confirmed_original_census_ID"]==59
    assert a["source_genotyped_bisexual_without_exact_original_census_ID"]==1
    assert [a["cohort_flow_by_arm"][arm]["genotyped_mature_outcome"] for arm in ("intact","half","all")]==[42,11,6]
    assert [a["cohort_flow_by_arm"][arm]["mature_outcome_known_but_selfing_untyped"] for arm in ("intact","half","all")]==[3,1,1]
    assert [a["cohort_flow_by_arm"][arm]["no_matching_mature_outcome_in_source"] for arm in ("intact","half","all")]==[9,3,5]
    assert a["cohort_flow_by_arm"]["all"]["known_mature_outcome_untyped_harvest_seeds_sum"]==113
    assert a["cohort_flow_by_arm"]["intact"]["known_mature_outcome_untyped_harvest_seeds_sum"]==126
    assert math.isclose(contrast_sensitivity(a,.95,0,0,0),75.37930079180079,rel_tol=1e-8)
    assert math.isclose(contrast_sensitivity(a,.95,.5,.5,0),55.19180079180079,rel_tol=1e-8)
    assert math.isclose(contrast_sensitivity(a,.95,0,0,25),137.8793007918008,rel_tol=1e-8)
    assert math.isclose(a["tipping_by_delta_original_assigned"]["0.95"]["typed_all_minus_intact"],38.47413588395731,rel_tol=1e-8)
    assert a["nested_typed_flower_and_typed_seed_bootstrap"]["0.95"]["all_minus_intact_genotype_complete_case_Beta_plant_bootstrap95"][0]>0
    print("v61 ALL SOURCE identity, treatment cohort, measured immature followup seed, complete-case and loss sensitivity constraints PASS")
