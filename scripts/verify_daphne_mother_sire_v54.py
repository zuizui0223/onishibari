#!/usr/bin/env python3
"""Check deterministic rerun of the source-pinned Daphne v54 aggregate."""
import argparse,json,math
from pathlib import Path
def compare(a,b,path="root"):
    if isinstance(b,dict):
        assert set(a)==set(b),(path,"keys")
        for k in b:compare(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(float,int)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_daphne_mother_sire_clustering_v54.json"))
    o=p.parse_args()
    a=json.loads(o.actual.read_text())
    b=json.loads(o.expected.read_text())
    compare(a,b)
    assert a["source_zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    assert [v["source_assigned_offspring"] for v in a["yearly"]]==[180,206]
    assert [v["sampled_mothers_with_eligible_H_offspring"] for v in a["yearly"]]==[39,33]
    for year in a["yearly"]:
        assert year["mother_equal_same_father_pair_fraction_observed"]>year["spatially_weighted_exact_margins_nulls"]["5"]["mother_equal_metric"]["null_95"][1]
    print("v0.54: source archive, mean/paired/null/acceptance and maternal-equal QA matched pinned aggregate")
