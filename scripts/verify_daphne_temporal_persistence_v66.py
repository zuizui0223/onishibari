#!/usr/bin/env python3
"""Validate fixed v66 original-source numerical aggregate and fresh rerun."""
import argparse,json,math
from pathlib import Path
def eq(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:eq(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),(path,"length")
        for i,(u,v) in enumerate(zip(a,b)):eq(u,v,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-8),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-26_daphne_temporal_persistence_v66.json"))
    z=p.parse_args();a=json.loads(z.actual.read_text());b=json.loads(z.expected.read_text())
    eq(a,b)
    h,f=a["by_sex"]
    assert h["sex"]=="H" and h["paired_genets"]==45 and h["positive_genets"]==[6,6] and h["positive_both"]==1
    assert h["flower_totals"]==[2404,1979] and h["fruit_totals"]==[31,22]
    assert f["sex"]=="F" and f["paired_genets"]==19 and f["positive_both"]==4
    assert h["joint_models"]["mixture"]["parameters"]["weight_persistent"]<1e-4
    assert h["comparisons"]["persistent_minus_independent"]["joint"]["mean"]<0
    lo,hi=h["comparisons"]["persistent_minus_independent"]["joint"]["genet_resample95"]
    assert lo<0<hi
    tr=h["new_genet_only_crossyear_transfer"]
    assert [t["unseen_genets"] for t in tr]==[26,60]
    assert tr[0]["models"]["beta_binomial"]["mean_logscore"]<tr[0]["models"]["binomial"]["mean_logscore"]
    assert tr[1]["models"]["beta_binomial"]["mean_logscore"]>tr[1]["models"]["binomial"]["mean_logscore"]
    print("v66 PASS original paired-genet counts, exact numerical aggregate, two-year conditional forecasts and unseen-genet transfers")
