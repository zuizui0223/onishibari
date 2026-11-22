#!/usr/bin/env python3
"""Verify exact source-locked aggregate v53 under deterministic RNG nulls."""
import argparse,json,math
from pathlib import Path
def cmp(a,b,path="root"):
    if isinstance(b,dict):
        assert set(a)==set(b),(path,"keys")
        for k in b:cmp(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):cmp(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-7),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",required=True,type=Path)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-22_daphne_sire_repeatability_v53.json"))
    v=p.parse_args()
    a=json.loads(v.actual.read_text())
    b=json.loads(v.expected.read_text())
    cmp(a,b)
    assert a["same_strict_H_genet_in_both_years"]==45
    assert a["repeatability"]["overlapping_sires_with_at_least_one_assigned_both_years"]==15
    assert a["repeatability"]["observed_successful_same_H_sires_2018"]==24
    assert a["repeatability"]["observed_successful_same_H_sires_2020"]==23
    print("DAPHNE V53: fixed archive, exact-year father risk sets, conditional nulls, 45 paired H genets: ALL numeric gates passed")
