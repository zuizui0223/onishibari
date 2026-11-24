#!/usr/bin/env python3
"""Privacy-preserving audit of treatment coding in original 2018/2020 maternal SSR records.
Print only redacted patterns, counts, and codebook indications, never original IDs.
"""
import collections,json,re,io,zipfile
from pathlib import Path
import pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,get,API,ZIP_SHA,stringify
src=data()
adult=src["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
for year in [2018,2020]:
    nm=f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv"
    r=src[nm].copy()
    for column in ["SeedID","Offspring","DNAID"]:
        val=r[column].astype(str)
        pattern=val.map(lambda s:re.sub(r"\d+","#",s))
        # Never output specific seed identifiers, even prefixes longer than four.
        short=val.map(lambda s:re.match(r"^[A-Za-z]{1,4}",s).group(0) if re.match(r"^[A-Za-z]{1,4}",s) else "_nonalphabetic_")
        print("LABEL_CODE_AUDIT_V54",json.dumps({
           "year":year,"column":column,"rows":len(r),
           "type":str(r[column].dtype),
           "redacted_label_patterns":{str(k):int(v) for k,v in pattern.value_counts().head(20).items()},
           "nonidentifying_alpha_prefix_counts":{str(k):int(v) for k,v in short.value_counts().head(20).items()},
           "first_character_frequency":{str(k):int(v) for k,v in val.str[0].value_counts().head(20).items()}},ensure_ascii=False))
    codes=adult[["GenetID2",f"Sex_{year}.g"]].copy()
    codes["Mother_key"]=stringify(codes.GenetID2)
    r["Mother_key"]=stringify(r.Mother)
    j=r.merge(codes[["Mother_key",f"Sex_{year}.g"]],on="Mother_key",how="left",validate="m:1")
    s=r.loc[r.Father.notna()&r.Probability.ge(.95)]
    print("MOTHER_SEX_SOURCE_V54",json.dumps({
        "year":year,
        "mother_sex_by_all_seed_rows":{str(k):int(v) for k,v in j[f"Sex_{year}.g"].fillna("_unknown").value_counts().items()},
        "mother_sex_by_p095_seed_rows":{str(k):int(v) for k,v in j.loc[j.Father.notna()&j.Probability.ge(.95),f"Sex_{year}.g"].fillna("_unknown").value_counts().items()},
        "unique_source_mothers":int(r.Mother.nunique()),
        "assigned_p095":int(len(s))
        },ensure_ascii=False))
# Source scripts: print file names and count of treatment-coding terms, not statements with IDs
record=json.loads(get(API))
fs=record["files"];en=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
sourcezip=get(en["Data for Figures_Analyses .zip"]["links"].get("content") or en["Data for Figures_Analyses .zip"]["links"]["self"])
with zipfile.ZipFile(io.BytesIO(sourcezip)) as z:
    import hashlib
    sample_members=[p for p in z.namelist() if Path(p).name=="Sample_DATA.201820_seed.csv"]
    assert len(sample_members)==1
    sample_raw=z.read(sample_members[0])
    samples=pd.read_csv(io.BytesIO(sample_raw))
    assert {"DNAID","Info"}<=set(samples)
    assert not samples.DNAID.duplicated().any()
    for year in (2018,2020):
        fname=f"NF{year}_2_primer16_res1.2.3_var1_150000_220702.csv"
        original=src[fname]
        jj=original.merge(samples[["DNAID","Info"]],on="DNAID",how="left",validate="one_to_one")
        assert len(jj)==len(original)
        assign=jj.loc[jj.Father.notna()&jj.Probability.ge(.95)]
        print("TREATMENT_SOURCE_JOIN_V54",json.dumps({
            "year":year,"archive_member_sha256":hashlib.sha256(sample_raw).hexdigest(),
            "sample_table_rows":len(samples),
            "source_csv_rows":len(original),
            "unmatched_info":int(jj.Info.isna().sum()),
            "all_treatment_counts":{str(k):int(v) for k,v in jj.Info.fillna("_unmatched").value_counts().items()},
            "p095_assigned_treatment_counts":{str(k):int(v) for k,v in assign.Info.fillna("_unmatched").value_counts().items()},
        },ensure_ascii=False))
    matches=[]
    for p in z.namelist():
        if not p.lower().endswith((".r",".rmd",".txt")):continue
        raw=z.read(p)
        if len(raw)>70000:continue
        txt=raw.decode("utf-8-sig",errors="ignore")
        terms={a:len(re.findall(a,txt,re.I)) for a in ["SeedID","Offspring","Day","Night","Natural","NF2018","NF2020","Treatment"]}
        matches.append({"file":Path(p).name[:80],"matched_terms":terms})
    print("TREATMENT_CODEBOOK_SOURCE_FILE_INDEX_V54",json.dumps(matches[:40],ensure_ascii=False))
