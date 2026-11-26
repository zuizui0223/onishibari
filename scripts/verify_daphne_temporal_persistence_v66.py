#!/usr/bin/env python3
"""Numerical, not bitwise, replay gate for v66.

Independent source-verified reruns differ slightly in finite-difference
L-BFGS-B convergence. Preserve exact counts, seeds, categories and MC
quantiles; use explicit small tolerances for fitted quantities. In the H
mixture vs independent comparison the models are effectively identical,
so exact counts of positive score differences of order 1e-6 are unstable.
These sign counts are diagnostics, never evidence of ecological direction.
An optimizer's success flag is reported, not treated as a biological result.
"""
import argparse,json,math
from pathlib import Path

def compare(a,b,path="root",notes=None):
    if notes is None: notes=[]
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:compare(a[k],b[k],path+"."+k,notes)
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),(path,"length")
        for i,(u,v) in enumerate(zip(a,b)):compare(u,v,f"{path}[{i}]",notes)
    elif path.endswith(".optimizer_success_count"):
        assert isinstance(a,int) and isinstance(b,int) and abs(a-b)<=1,(path,a,b)
        if a!=b:notes.append((path,b,a,"convergence diagnostic"))
    elif path.startswith("root.by_sex[0].comparisons.mixture_minus_independent.") and path.endswith(".positive_genets"):
        assert isinstance(a,int) and 0<=a<=45
        if a!=b:notes.append((path,b,a,"near-zero score sign count"))
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a),(path,a)
        if path.endswith(".concentration"):
            assert math.isclose(a,b,rel_tol=1e-4,abs_tol=1e-5),(path,a,b)
        elif ".simulated_overlap." in path:
            # Replayed MC draws use the same seed; minor changes in
            # fitted beta parameters may affect occasional simulated counts.
            assert math.isclose(a,b,rel_tol=0,abs_tol=.005),(path,a,b)
        else:
            assert math.isclose(a,b,rel_tol=1e-6,abs_tol=2e-5),(path,a,b)
    else:
        assert a==b,(path,a,b)
    return notes

def validate(a):
    h,f=a["by_sex"]
    assert h["sex"]=="H" and h["paired_genets"]==45
    assert h["positive_genets"]==[6,6] and h["positive_both"]==1
    assert h["flower_totals"]==[2404,1979] and h["fruit_totals"]==[31,22]
    assert f["sex"]=="F" and f["paired_genets"]==19 and f["positive_both"]==4
    assert a["seed"]==20261126 and a["bootstrap_draws"]==5000 and a["overlap_simulations"]==10000
    for group in (h,f):
        assert all(v["optimizer_success_count"]>=group["paired_genets"]-1
                   for v in group["leave_one_genet_out"].values())
    assert h["joint_models"]["mixture"]["parameters"]["weight_persistent"]<1e-4
    c=h["comparisons"]["persistent_minus_independent"]["joint"]
    assert c["mean"]<0 and c["genet_resample95"][0]<0<c["genet_resample95"][1]
    for c in h["comparisons"]["mixture_minus_independent"].values():
        assert abs(c["mean"])<1e-4
    tr=h["new_genet_only_crossyear_transfer"]
    assert [t["unseen_genets"] for t in tr]==[26,60]
    assert tr[0]["models"]["beta_binomial"]["mean_logscore"]<tr[0]["models"]["binomial"]["mean_logscore"]
    assert tr[1]["models"]["beta_binomial"]["mean_logscore"]>tr[1]["models"]["binomial"]["mean_logscore"]

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-26_daphne_temporal_persistence_v66.json"))
    z=p.parse_args();a=json.loads(z.actual.read_text());b=json.loads(z.expected.read_text())
    validate(a);validate(b);notes=compare(a,b)
    print("Numerical replay differences in diagnostic flags/sign counts:",json.dumps(notes))
    print("v66 PASS exact source counts and tolerance-bounded numerical estimates; inferential direction and transfer limitations unchanged")
