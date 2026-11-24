#!/usr/bin/env python3
"""Assert pinned v56 source audit, confidence thresholds, and numerical estimates."""
import argparse,json,math
from pathlib import Path
def compare(a,b,path="root"):
    if isinstance(b,dict):
        assert set(a)==set(b),(path,"keys")
        for key in b:compare(a[key],b[key],path+"."+key)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_daphne_confidence_prior_v56.json"))
    z=p.parse_args()
    a=json.loads(z.actual.read_text())
    b=json.loads(z.expected.read_text())
    compare(a,b)
    assert a["source_zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    assert [x["strict_H_eligible_father_genets"] for x in a["by_year"]]==[105,71]
    expected_n=[[110,107,106,101],[131,129,128,126]]
    for i,year in enumerate(a["by_year"]):
        assert [x["n_confident_strictH_natural_offspring"] for x in year["natural_only_high_confidence_thresholds"]]==expected_n[i]
        for t in year["natural_only_high_confidence_thresholds"]:
            assert t["mother_equal_sibling_observed"]>t["null_fixed_sire_totals_spatial_lambda5"]["null_95"][1]
            assert t["adversarial_paternity_relabelling_bound"]["minimum_reassignment_for_upper_null"]>0
    print("v56 ALL numeric and source provenance checks passed")
