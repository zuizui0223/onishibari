#!/usr/bin/env python3
"""Complete pinned v63 source and numerical verification; no raw ID output."""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
def eq(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:eq(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):eq(x,y,f"{path}[{i}]")
    elif isinstance(b,(int,float)) and not isinstance(b,bool):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--actual",type=Path,required=True)
    parser.add_argument("--expected",type=Path,default=Path("data/2026-11-24_pulsatilla_origin_uncertainty_v63.json"))
    a=parser.parse_args()
    actual=json.loads(a.actual.read_text());saved=json.loads(a.expected.read_text())
    eq(actual,saved)
    assert actual["source_first_81_enrolled_bisexual"]==81
    assert actual["exact_enrollment_linked_typed_n"]=={"intact":42,"half":11,"all":6}
    assert actual["typed_seed_counts_by_arm"]["all"]=={"genotyped_seeds":49,"selfed_genotyped_seeds":0}
    assert actual["monte_carlo_draws"]==20000
    assert actual["reused_v61_source_sha256"]["phenotypes"]=="e49fef29cefa6e1ed419500a0e312a49524f04876cca5eb21dda14657fe11a98"
    z=actual["delta_x_missing_cap_uncertainty"]["0.95"]
    assert abs(z["0"]["v62_point_id_lower"]-3.6173834653001293)<1e-8
    assert abs(z["0"]["structural_emasculation_point_lower"]-12.563216798633462)<1e-8
    beta=z["50"]["uncertainty"]["Beta_per_flower_unknown_selfing_all_arms_seed_only"]
    structural=z["50"]["uncertainty"]["structural_emasculation_zero_within_flower_selfing_seed_only"]
    assert beta["lower_bound_samples"]["mc_percentile95"][1]<0
    assert structural["lower_bound_samples"]["mc_percentile95"][0]>0
    assert structural["joint_probability_bounds_strictly_positive"]>.999
    assert len(actual["delta_x_missing_cap_uncertainty"])==7
    for case in actual["delta_x_missing_cap_uncertainty"].values():
        assert len(case)==7
    print("v63 source SHA, exact enrollment, 0/49 full-removal typed selfing, two structural assumptions and 20k-draw full numerical gates PASS")
