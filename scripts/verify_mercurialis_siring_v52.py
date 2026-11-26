#!/usr/bin/env python3
"""Compare original-source v0.52 derived aggregate against pinned summary."""
import argparse,json,math
from pathlib import Path
def check(a,b,path="root"):
    if isinstance(b,dict):
        assert set(a)==set(b),(path,"keys")
        for k in b:check(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(av,bv) in enumerate(zip(a,b)):check(av,bv,f"{path}[{i}]")
    elif isinstance(b,float):
        assert math.isfinite(a) and math.isclose(a,b,rel_tol=1e-5,abs_tol=1e-6),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-22_mercurialis_siring_decomposition_v52.json"))
    q=p.parse_args()
    a=json.loads(q.actual.read_text());b=json.loads(q.expected.read_text())
    check(a,b)
    assert a["source_parents"]==181 and a["complete_bisexual_size_subset_n"]==167
    assert a["formula_identity_audit"]["paired_mismatch_above_1e_minus_5"]==1
    print("V52 pinned original archive + same-genet decomposition aggregate: ALL checks passed")
