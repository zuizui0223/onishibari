#!/usr/bin/env python3
"""Check full original-source reproducibility of independent Oenothera v57."""
import argparse,json,math
from pathlib import Path
def eq(a,b,path="root"):
    if isinstance(b,dict):
        assert set(a)==set(b),(path,"keys")
        for k in b:eq(a[k],b[k],path+"."+k)
    elif isinstance(b,list):
        assert len(a)==len(b),(path,"length")
        for i,(x,y) in enumerate(zip(a,b)):eq(x,y,f"{path}[{i}]")
    elif isinstance(b,float):
        assert isinstance(a,(int,float)) and math.isfinite(a) and math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9),(path,a,b)
    else:
        assert a==b,(path,a,b)
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--actual",type=Path,required=True)
    ap.add_argument("--expected",type=Path,default=Path("data/2026-11-24_oenothera_genotype_only_v57.json"))
    p=ap.parse_args()
    a=json.loads(p.actual.read_text())
    b=json.loads(p.expected.read_text())
    eq(a,b)
    q=a["Mendelian_identity"]
    assert q["parents_sampled"]==323 and q["offspring_genotyped"]==521
    assert q["mother_genotype_inconsistent_loci_total"]==1
    assert q["offspring_with_zero_compatible_father"]==207
    assert q["offspring_with_unique_zero_mismatch_father"]==274
    assert q["offspring_with_multiple_zero_mismatch_fathers"]==40
    assert a["genotype_workbook_sha256"]=="8348dd28cbc4532cea22f80676d241cea322527bec35e8ec1995edd812be30e7"
    assert a["fruit_MLTR_source_sha256"]=="c703d56287b62f7550931e46f54591fa7a872438420bf401814eea3186320eb6"
    assert all(x["genotyped_families_with_original_MLTR_match"]==60 for x in a["uncalibrated_genotype_only_weight_sensitivity"])
    print("V57 source MD5/SHA, Mendelian 0-as-missing and all sensitivity/paired mother numeric checks PASS")
