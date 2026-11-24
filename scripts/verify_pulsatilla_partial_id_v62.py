#!/usr/bin/env python3
"""Full original-source downstream numerical verification for v0.62 partial ID."""
import argparse,json,math
from pathlib import Path
def check(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:check(a[k],b[k],path+"."+str(k))
    elif isinstance(b,list):
        assert isinstance(a,list) and len(a)==len(b),(path,"len")
        for i,(x,y) in enumerate(zip(a,b)):check(x,y,f"{path}[{i}]")
    elif isinstance(b,(float,int)) and not isinstance(b,bool):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_pulsatilla_partial_identification_v62.json"))
    z=p.parse_args()
    a=json.loads(z.actual.read_text());b=json.loads(z.expected.read_text())
    check(a,b)
    assert a["published_source_census_has_pistil_ovule_counts_for_unmatched_flowers"] is False
    q=a["δ_sensitivity"]["0.95"]
    assert 21.70<q["boundaries"]["max_intact_missing_effective_seeds_for_guaranteed_positive"]<21.71
    rows=q["symmetric_missing_outcome_hypothetical_cap_scenarios"]
    assert len(rows)==13 and len(a["δ_sensitivity"])==6
    t10=next(x for x in rows if x["hypothetical_cap_effective_seeds_per_untracked_all_flower"]==10)
    t25=next(x for x in rows if x["hypothetical_cap_effective_seeds_per_untracked_all_flower"]==25)
    assert t10["sharp_treatment_mean_difference_lower"]>0
    assert t25["sharp_treatment_mean_difference_lower"]<0<t25["sharp_treatment_mean_difference_upper"]
    dose=q["dose_response_sharp_pairwise_at_shared_lost_cap"][0]
    assert dose["half_vs_intact"]["sign"]=="guaranteed_positive"
    assert dose["all_vs_half"]["sign"]=="unknown"
    print("v62 complete original-enrollment source, δxK, 3-arm sharp-bounds numerical verification PASS")
