#!/usr/bin/env python3
"""v55 numerical validation against source-pinned published archive aggregate."""
import argparse,json,math
from pathlib import Path
def verify(actual,saved,path="root"):
    if isinstance(saved,dict):
        assert set(actual)==set(saved),(path,"keys")
        for k in saved:verify(actual[k],saved[k],path+"."+k)
    elif isinstance(saved,list):
        assert len(actual)==len(saved),(path,"length")
        for i,(a,b) in enumerate(zip(actual,saved)):verify(a,b,f"{path}[{i}]")
    elif isinstance(saved,float):
        assert isinstance(actual,(int,float)) and math.isfinite(actual) and math.isclose(actual,saved,rel_tol=1e-8,abs_tol=1e-9),(path,actual,saved)
    else:
        assert actual==saved,(path,actual,saved)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_daphne_treatment_specific_sire_v55.json"))
    a=p.parse_args()
    actual=json.loads(a.actual.read_text());saved=json.loads(a.expected.read_text())
    verify(actual,saved)
    assert actual["sample_treatment_csv_sha256"]=="739dafbbee07de46ee3c6a10086b754f1e33b697d111a9ec8ee22d5a8a35d9b1"
    assert [x["source_all_genotyped_treatment_counts"]["Natural"] for x in actual["years"]]==[158,199]
    assert [x["strict_H_father_p095_by_treatment"]["Natural"] for x in actual["years"]]==[110,131]
    for y in actual["years"]:
        t=next(x for x in y["by_treatment"] if x["treatment"]=="Natural")
        assert t["mothers_with_2_or_more"]>=18
        for l in ["5","10","20"]:
            assert t["mother_equal_observed"]>t["spatial_fixed_margin_nulls"][l]["mother_equal_null"]["null_95"][1]
    print("DAPHNE V55: original archives, DNAID treatments, all six cohorts, spatial nulls and numerical aggregate exactly reproduced")
