#!/usr/bin/env python3
"""Reproduce full v58 source-backed numerical source and marker gates."""
import argparse,json,math
from pathlib import Path
def rec(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:rec(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):rec(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_oenothera_locus_robustness_v58.json"))
    z=p.parse_args()
    a=json.loads(z.actual.read_text());b=json.loads(z.expected.read_text())
    rec(a,b)
    q=a["strict_mendelian_uniqueness_stability"]
    assert q["full_7_locus_unique_count_5_typed_min"]==274
    assert q["full_unique_same_sire_under_all_seven_single_locus_omissions_4_typed_min"]==108
    assert q["full_unique_ever_ambiguous_after_one_locus_drop"]==166
    assert a["wrong_mother_negative_control"]["actual_mother_offspring_with_one_or_more_mendelian_exclusions"]==1
    assert len(a["heldout_locus_predictive_evidence"])==3
    assert all(h["mother_bootstrap95"][0]>0 for h in a["heldout_locus_predictive_evidence"])
    assert a["pinned_genotype_sha256"]=="8348dd28cbc4532cea22f80676d241cea322527bec35e8ec1995edd812be30e7"
    print("v58 ORIGINAL SOURCE HASH, 7 SSR LOLO, HOLDOUT and 26-mother negative controls PASS")
