#!/usr/bin/env python3
"""v0.28 Source-original 2018/2020 Daphne father assignments: sex of true sires.

Uses author-released 16-SSR parentage, joined to original adult genet/year
sex codes. A direct father assignment (posterior >=.95) is stronger evidence
of male reproduction than flower morphology. HF-coded adult sex is ambiguous
until this is checked. No individual parent IDs / genotypes / coordinates logged.
Not a separate experiment, and not D. pseudomezereum.
"""
from __future__ import annotations
import hashlib,io,json,urllib.request,zipfile
from pathlib import Path
import pandas as pd,numpy as np

API="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
HASHES={
"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv":"fefd5e4a31880fbf1c1e70f55b897e0ac37e78837ebd40408b7a851f94016d19",
"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv":"957c293549d54b3d3a1ef5683359091c207844450215ddee11580cc23faa1d78",
"GenetID_geno_phenotype_1820Nopadult_211109.csv":"97790d8b3467d2f9e860f97607aa60525a17ddf1b97b0645fdb1b179ce99135e"
}
OUT=Path("analysis_results_daphne_HF_father_source_v28")

def get(url,cap=2100000):
    request=urllib.request.Request(url,headers={"User-Agent":"original-Daphne-HF-sire-assignment-v28","Accept":"*/*"})
    with urllib.request.urlopen(request,timeout=80) as res:b=res.read(cap+1)
    if len(b)>cap:raise RuntimeError("source file too large")
    return b
def data():
    j=json.loads(get(API))
    fs=j["files"];es=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    a=es["Data for Figures_Analyses .zip"]["links"]
    blob=get(a.get("content") or a["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    d={}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for p in z.namelist():
            name=Path(p).name
            if name in HASHES:
                b=z.read(p)
                assert hashlib.sha256(b).hexdigest()==HASHES[name],name
                d[name]=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig")
    assert set(d)==set(HASHES)
    return d
def stringify(s):
    return s.astype(str).str.strip().str.replace(r"\.0$","",regex=True)
def cats(label):
    x=str(label)
    if x=="H":return "H_exact"
    if x=="HF":return "HF_exact"
    if x=="F":return "F_exact"
    if x=="NF":return "NF_exact"
    if "H" in x:return "male_mixed_code"
    if x=="dead":return "dead"
    return "other_unknown"
def by_year(df,adult,year):
    sex=f"Sex_{year}.g"
    fl=f"NaturalFlowerNo_{year}.g"
    assert sex in adult and fl in adult
    source_n=len(df)
    f=df[df.Father.notna()&pd.to_numeric(df.Probability,errors="coerce").ge(.95)].copy()
    f["key"]=stringify(f.Father)
    ad=adult[["GenetID2",sex,fl]].copy()
    assert not ad.GenetID2.duplicated().any()
    ad["key"]=stringify(ad.GenetID2)
    assert not ad.key.duplicated().any()
    j=f.merge(ad,on="key",how="left",validate="m:1",indicator=True)
    assert len(j)==len(f) and j["_merge"].eq("both").all()
    j["category"]=j[sex].map(cats)
    j["flowers"]=pd.to_numeric(j[fl],errors="coerce")
    summary={}
    for cat,g in j.groupby("category"):
        summary[cat]={
           "assigned_offspring":int(len(g)),
           "unique_father_genets":int(g.key.nunique()),
           "n_genet_flower_count_nonmissing":int(g[["key","flowers"]].drop_duplicates("key").flowers.notna().sum()),
           "n_unique_father_genets_with_positive_flowers":int(g.loc[g.flowers.gt(0),"key"].nunique())}
    candidates=ad.copy()
    candidates["category"]=candidates[sex].map(cats)
    candidates["flowers"]=pd.to_numeric(candidates[fl],errors="coerce")
    candidate_stats={}
    for cat,g in candidates.groupby("category"):
        candidate_stats[cat]={
          "adult_genets":int(len(g)),
          "adult_genets_with_positive_flowers":int(g.flowers.gt(0).sum()),
          "source_total_flowers_where_positive":int(g.flowers.loc[g.flowers.gt(0)].sum()),
          "with_assigned_at_95percent":int(g.key.isin(j.key).sum())}
    out={
       "year":year,
       "parentage_offspring_all":source_n,
       "parentage_assigned_95":len(f),
       "assigned_mother_genets":int(j.Mother.nunique()),
       "assigned_father_genets":int(j.key.nunique()),
       "paternal_sex_code_categories":summary,
       "year_source_adult_sex_categories":candidate_stats,
       "n_assigned_without_known_adult":int(j["_merge"].ne("both").sum()),
       "functional_HF_fathers_verified":bool(summary.get("HF_exact",{}).get("unique_father_genets",0)>0),
       "HF_sire_conclusion_is_from_95percent_genetic_parentage_not_visual_sex_alone":True
    }
    return out

def main():
    OUT.mkdir(exist_ok=True)
    d=data()
    adult=d["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(adult)==492
    infiles={2018:"NF2018_2_primer16_res1.2.3_var1_150000_220702.csv",
             2020:"NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"}
    out={"status":"author_source_parentage_potential_HF_fathers_v28",
       "source_doi":"10.5281/zenodo.8253921",
       "article_doi":"10.1002/ajb2.16239",
       "original_zip_SHA256":ZIP_SHA,"CSV_SHA256":HASHES,
       "by_year":[by_year(d[f],adult,year) for year,f in infiles.items()],
       "inferential_boundaries":{
         "only_assigned_seeds_at_95percent":True,
         "father_counts_are_not_independent_seed_replicates":True,
         "composite_sex_codes_not_assumed_functional_H_without_father_confirmation":True,
         "failure_to_assign_father_is_not_proof_of_no_fertility":True,
         "no_Daphne_pseudomezereum_experimental_data":True,
         "local_fruit_density_association_v27_not_a_causal_effect":True}}
    assert sum(x["parentage_offspring_all"] for x in out["by_year"])==584
    assert sum(x["parentage_assigned_95"] for x in out["by_year"])==470
    (OUT/"daphne_2018_2020_father_genet_sex_code_v28.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("DAPHNE_HF_TRUE_PATERNITY_V28",json.dumps(out,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
