#!/usr/bin/env python3
"""Verify pinned original Rhodes genetic source and all v59 permutations."""
import argparse,json,math
from pathlib import Path
def cmp(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:cmp(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):cmp(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(float,int)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-10),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_oenothera_child_specific_v59.json"))
    q=p.parse_args()
    a=json.loads(q.actual.read_text());b=json.loads(q.expected.read_text())
    cmp(a,b)
    assert a["source_genotypes_sha256"]=="8348dd28cbc4532cea22f80676d241cea322527bec35e8ec1995edd812be30e7"
    assert a["source_correlated_paternity_sha256"]=="c703d56287b62f7550931e46f54591fa7a872438420bf401814eea3186320eb6"
    assert (a["parent_genets"],a["offspring_genotyped"],a["independent_mothers"],a["loci"])==(323,521,26,7)
    assert a["informative_mother_offspring_locus_pairs"]==3445
    assert a["results"]["same_mother_same_treatment"]["design"]["groups"]==68
    for key,v in a["results"].items():
        for z in v["epsilon_sensitivity"]:
            assert z["child_specific_increment_over_mother_group_null_nats"]>0
            assert z["mother_bootstrap_95_increment"][0]>0
            assert z["conditional_permutation_upper_tail_p"]<=.001
    print("v0.59 SOURCE VERIFIED: family + treatment membership, 7 marker holdout, 2999 permutations, 26-mother confidence bounds")
