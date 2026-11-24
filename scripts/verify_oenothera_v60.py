#!/usr/bin/env python3
"""Source archive and complete numeric verification for Oenothera v60."""
import argparse,json,math
from pathlib import Path
def check(a,b,path="root"):
    if isinstance(b,dict):
        assert isinstance(a,dict) and set(a)==set(b),(path,"keys")
        for k in b:check(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):check(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--actual",type=Path,required=True)
    p.add_argument("--expected",type=Path,default=Path("data/2026-11-24_oenothera_paternal_locus_v60.json"))
    v=p.parse_args()
    a=json.loads(v.actual.read_text())
    b=json.loads(v.expected.read_text())
    check(a,b)
    assert a["source_genotypes_sha256"]=="8348dd28cbc4532cea22f80676d241cea322527bec35e8ec1995edd812be30e7"
    assert a["source_paternity_reference_sha256"]=="c703d56287b62f7550931e46f54591fa7a872438420bf401814eea3186320eb6"
    assert (a["parents"],a["offspring"],a["maternal_genets"],a["source_SSR_loci"])==(323,521,26,7)
    assert a["maternal_offspring_informative_locus_pairs"]==3445
    assert a["coherent_permutation_metmorphic_max_absolute_score_difference"]<1e-12
    assert a["permutation_design"]["nonmaternal_adult_candidates_permuted"]==297
    assert a["permutation_design"]["nonmaternal_prefix_groups"]==2
    for mode in ("global","within_prefix"):
        for x in a["results"][mode]["epsilon_sensitivity"]:
            assert x["child_specific_paternal_identity_coherence_increment_nats"]>0
            assert x["equal_mother_resample_95increment"][0]>0
            assert x["conditional_null_upper_p"]<=.003
            assert x["number_null_locus_permutations"]==399
    print("v60 VERIFIED: Zenodo/Dryad original-source SHA, coherent father-label invariance, 399 paternal SSR locus permutations and 26-mother numerical gates")
