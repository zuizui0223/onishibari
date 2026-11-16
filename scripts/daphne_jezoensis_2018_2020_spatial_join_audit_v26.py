#!/usr/bin/env python3
"""v0.26: identify if original 2018/2020 Daphne maternal fruit records
can be joined without ambiguity to genotyped adult locations and local
potential pollen donors. Metadata diagnostics only, no inferential model.

Key safety: do NOT infer same individual from numeric overlap alone;
require original fruit/flower counts or sex label as validation and
inspect ambiguities before gene-flow / neighborhood claims.
"""
from __future__ import annotations
import hashlib, io, json, urllib.request, zipfile
from pathlib import Path
import pandas as pd
import numpy as np

RECORD="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
HASHES={
 "Nopporo_control.csv":"8b36a936ac7880d6fd7abb9928cd1b0677ffc75e80f1abe6a967c0e67118b07b",
 "GenetID_geno_phenotype_1820Nopadult_211109.csv":"97790d8b3467d2f9e860f97607aa60525a17ddf1b97b0645fdb1b179ce99135e"
}
OUT=Path("analysis_results_daphne_2018_2020_local_mates_v26")

def get(url, limit=2_000_000):
    req=urllib.request.Request(url,headers={"User-Agent":"published-Daphne-original-mate-match-audit-v26","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as r:b=r.read(limit+1)
    if len(b)>limit:raise RuntimeError("source read exceeds cap")
    return b
def files():
    record=json.loads(get(RECORD))
    fs=record["files"]
    entries=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    item=entries["Data for Figures_Analyses .zip"]["links"]
    archive=get(item.get("content") or item["self"])
    assert hashlib.sha256(archive).hexdigest()==ZIP_SHA
    out={}
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        for p in z.namelist():
            n=Path(p).name
            if n in HASHES:
                raw=z.read(p)
                assert hashlib.sha256(raw).hexdigest()==HASHES[n]
                out[n]=pd.read_csv(io.BytesIO(raw))
    assert set(out)==set(HASHES)
    return out
def vals(x):
    return {str(k):int(v) for k,v in x.value_counts(dropna=False).items()}
def year_profile(ctrl,ad,year):
    c=ctrl[(ctrl.Year==year)&(ctrl.Treatment=="Control")&ctrl.Sex.isin(["F","H","HF"])].copy()
    col_sex=f"Sex_{year}.g"
    col_fl=f"NaturalFlowerNo_{year}.g"
    col_fr=f"NaturalFruitNo_{year}.g"
    assert {col_sex,col_fl,col_fr,"GenetID2","RametNo","Ramet","GenetID_extraction"}<=set(ad)
    c["flow"]=pd.to_numeric(c.FlowerNo,errors="coerce")
    c["fru"]=pd.to_numeric(c.FruitNo,errors="coerce")
    ret={"year":year,"control_observed_sex_counts":vals(c.Sex),
      "control_n":len(c),"valid_fruit_positive_flower":int((c.flow.gt(0)&c.fru.ge(0)&c.fru.le(c.flow)).sum()),
      "control_non_numeric_flower":int(c.flow.isna().sum()),
      "control_non_numeric_fruit":int(c.fru.isna().sum()),
      "control_year_ramet_unique":int(c.RametID.nunique()),
      "adult_sex_counts":vals(ad[col_sex]),
      "adult_flower_nonnumeric":int(pd.to_numeric(ad[col_fl],errors="coerce").isna().sum()),
      "adult_fruit_nonnumeric":int(pd.to_numeric(ad[col_fr],errors="coerce").isna().sum()),
      "candidate_id_matches":[]}
    # Check unique ID mappings of each plausible ID column, validating both
    # naturally observed flower and fruit counts, not just same integers.
    for idcol in ["GenetID2","RametNo","Ramet","GenetID_extraction","DNAID"]:
        if idcol not in ad:continue
        dat=ad[[idcol,col_sex,col_fl,col_fr]].copy()
        left=c[["RametID","Sex","flow","fru"]].copy()
        left["id_compare"]=left.RametID.astype(str).str.strip().str.replace(r"\.0$","",regex=True)
        dat["id_compare"]=dat[idcol].astype(str).str.strip().str.replace(r"\.0$","",regex=True)
        nonblank=dat.id_compare.ne("")&dat.id_compare.ne("nan")&dat.id_compare.ne("None")
        dat=dat[nonblank]
        unique=dat.id_compare.nunique()
        dup=dat.id_compare.duplicated().sum()
        # Avoid exploding repeated adult key matches; report both uniqueness
        # and validation with unique IDs only.
        unique_ids=dat.id_compare.value_counts().loc[lambda x:x.eq(1)].index
        dat=dat[dat.id_compare.isin(unique_ids)]
        merged=left.merge(dat,on="id_compare",how="left",validate="m:1")
        hit=merged[idcol].notna()
        af=pd.to_numeric(merged[col_fl],errors="coerce")
        ar=pd.to_numeric(merged[col_fr],errors="coerce")
        nonmissing=hit&left.flow.notna()&left.fru.notna()&af.notna()&ar.notna()
        exact=nonmissing&merged.flow.eq(af)&merged.fru.eq(ar)
        sexmatch=hit&merged.Sex.astype(str).eq(merged[col_sex].astype(str))
        # Count plus alternative where fruit is measured at genus level and not
        # identical to ramet; too many equal numeric zeros can match by accident.
        ret["candidate_id_matches"].append({
          "adult_column":idcol,
          "adult_key_unique":int(unique),
          "adult_key_duplicate_rows":int(dup),
          "row_id_match":int(hit.sum()),
          "joint_valid_flower_and_fruit":int(nonmissing.sum()),
          "both_source_counts_identical":int(exact.sum()),
          "sex_label_identical":int(sexmatch.sum()),
          "all_three_identical":int((exact&sexmatch).sum()),
          "n_rows_original_female_with_ID_hit":int((hit&left.Sex.eq("F")).sum())
        })
    return ret

def main():
    OUT.mkdir(exist_ok=True)
    ff=files()
    c=ff["Nopporo_control.csv"]
    a=ff["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(c)==728 and len(a)==492
    report={
      "status":"original_Daphne_spatial_ID_eligibility_audit_v26_no_effect_fit",
      "author_article_doi":"10.1002/ajb2.16239",
      "author_source_doi":"10.5281/zenodo.8253921",
      "zip_sha256":ZIP_SHA,
      "original_CSV_sha256":HASHES,
      "source_control_rows":len(c),
      "source_adult_rows":len(a),
      "unique_genetIDs":int(a.GenetID2.nunique()),
      "n_adult_location_pair_nonmissing":int(a[["East.g","North.g"]].notna().all(axis=1).sum()),
      "adult_columns":[str(k) for k in a.columns],
      "original_control_sex_counts":vals(c.Sex),
      "original_control_treatment_counts":vals(c.Treatment),
      "year_diagnostics":[year_profile(c,a,y) for y in [2018,2020]],
      "conditions":{
        "possible_identifier_overlap_not_a_validated_genet_join":True,
        "need_1to1_ID_mapping_and_phenotype_concordance":True,
        "need_reproductive_cohort_and_sex_definition":True,
        "local_neighbors_not_pollen_donors_without_compatibility_evidence":True,
        "original_2015_pheno_source_is_different_year":True,
        "no_D_pseudomezereum_empirical_results":True}}
    (OUT/"daphne_local_mate_join_eligibility_v26.json").write_text(
      json.dumps(report,ensure_ascii=False,indent=2,default=str)+"\n")
    printable={k:v for k,v in report.items() if k not in ("adult_columns",)}
    print("DAPHNE_SPATIAL_JOIN_AUDIT_V26",json.dumps(printable,ensure_ascii=False,default=str),flush=True)
if __name__=="__main__":main()
