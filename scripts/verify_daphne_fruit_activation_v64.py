#!/usr/bin/env python3
"""Numerical source-archive gate for H/F Daphne jezoensis fruit activation v64."""
import argparse,json,math
from pathlib import Path
def compare(a,b,where="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(where,"keys")
        for k in b:compare(a[k],b[k],where+"."+str(k))
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),(where,"list_length")
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f"{where}[{i}]")
    elif isinstance(b,(int,float)) and not isinstance(b,bool):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-9),(where,a,b)
    else:
        assert a==b,(where,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_daphne_fruit_activation_v64.json"))
    x=p.parse_args()
    observed=json.loads(x.actual.read_text())
    saved=json.loads(x.expected.read_text())
    compare(observed,saved)
    assert observed["original_zip_sha256"]=="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
    expected={
        (2018,"H"):(105,4765,53,14,23),
        (2018,"F"):(63,619,255,37,34),
        (2020,"H"):(71,2573,28,9,14),
        (2020,"F"):(57,719,205,32,28),
    }
    for d in observed["sex_year"]:
        key=(d["year"],d["morph_exact_label"])
        observed_tuple=(d["n_original_genets_valid"],d["flowers_total"],d["fruits_total"],
                d["fruit_positive_genets"],d["single_largest_fruit_genet_excluded"]["n_original_fruits_removed_by_one_genet"])
        assert observed_tuple==expected[key],(key,observed_tuple)
        assert d["null_tests"]["all_flowers_within_morph_year"]["summary"]["fruit_positive"]["null_tail_p_one_sided"]<.001
        assert d["null_tests"]["flower_count_quartile_conditioned"]["summary"]["fruit_positive"]["null_tail_p_one_sided"]<.001
        assert d["null_tests"]["spatial_10m_shift_0m_plus_flower_half"]["summary"]["fruit_positive"]["null_tail_p_one_sided"]<.006
        assert d["single_largest_fruit_genet_excluded"]["flower_quartile_conditional_null"]["fruit_positive"]["null_tail_p_one_sided"]<.003
        assert d["null_tests"]["all_flowers_within_morph_year"]["summary"]["draws"]==6000
    h=observed["strict_H_temporal_recurrence"]
    assert h["same_H_genets_with_valid_natural_fruit_both_years"]==45
    assert h["fruit_positive_H_2018"]==h["fruit_positive_H_2020"]==6
    assert h["positive_both_years"]==1
    assert h["hypergeometric_one_sided_repeat_excess_p"]>.5
    print("v64 source SHA + H/F x year fruit counts, 6000-rep spatial/flower nulls, largest genet ablation and temporal H check PASS")
