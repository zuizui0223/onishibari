#!/usr/bin/env python3
"""Fresh original ZIP cohort and paternal sampling audit, v67.
No coordinates, ID values, raw genotype rows, or individual counts exported.
Do not expand paternal sampled offspring to population output; maternal fruit
counts and source seed sampling frames may differ. Aggregate evidence only.
"""
import io,zipfile,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from daphne_jezoensis_2018_2020_HF_fathers_v28 import get,API,ZIP_SHA,HASHES,stringify
from daphne_jezoensis_treatment_stratified_paternity_v55 import SAMPLE_SHA,SAMPLES_NAME
OUT=Path("analysis_results_daphne_v67_source_audit")
ADULT="GenetID_geno_phenotype_1820Nopadult_211109.csv"
def source():
    # Use one publicly accessible official file download; verify exact bytes.
    url="https://zenodo.org/records/8253921/files/Data%20for%20Figures_Analyses%20.zip?download=1"
    try:
        blob=get(url)
    except OSError:
        j=json.loads(get(API));fs=j["files"]
        e=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
        links=e["Data for Figures_Analyses .zip"]["links"]
        blob=get(links.get("content") or links["self"])
    assert hashlib.sha256(blob).hexdigest()==ZIP_SHA
    hashes={**HASHES,SAMPLES_NAME:SAMPLE_SHA}
    files={}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for member in z.namelist():
            name=Path(member).name
            if name not in hashes:continue
            raw=z.read(member);assert hashlib.sha256(raw).hexdigest()==hashes[name],name
            assert name not in files
            files[name]=pd.read_csv(io.BytesIO(raw),encoding="utf-8-sig")
    assert set(files)==set(hashes)
    return files

def year(files,yr):
    adult=files[ADULT].copy()
    assert len(adult)==492
    adult["key"]=stringify(adult.GenetID2);assert adult.key.is_unique
    adult["sx"]=adult[f"Sex_{yr}.g"].astype(str)
    adult["flowers"]=pd.to_numeric(adult[f"NaturalFlowerNo_{yr}.g"],errors="coerce")
    adult["fruits"]=pd.to_numeric(adult[f"NaturalFruitNo_{yr}.g"],errors="coerce")
    adult["valid"]=adult.flowers.gt(0)&adult.fruits.ge(0)&adult.fruits.le(adult.flowers)
    adult["locvalid"]=adult[["East.g","North.g"]].apply(pd.to_numeric,errors="coerce").notna().all(axis=1)
    h=adult[adult.sx.eq("H")&adult.valid].set_index("key")
    risk=adult[adult.sx.eq("H")&adult.flowers.gt(0)&adult.locvalid]
    assert set(h.index)==set(risk.key),"v55 sire risk set differs from v64 maternal cohort"
    assert (len(h),int(h.fruits.sum()),int(h.flowers.sum()))=={2018:(105,53,4765),2020:(71,28,2573)}[yr]
    sample=files[SAMPLES_NAME][["DNAID","Info"]]
    p=files[f"NF{yr}_2_primer16_res1.2.3_var1_150000_220702.csv"].merge(sample,on="DNAID",validate="one_to_one")
    assert p.Info.notna().all()
    natural=p[p.Info.eq("Natural")].copy()
    natural["mk"]=stringify(natural.Mother);natural["fk"]=stringify(natural.Father)
    natural["prob"]=pd.to_numeric(natural.Probability,errors="coerce")
    adults=adult.set_index("key")
    assert natural.mk.isin(adults.index).all()
    natural["mother_sx"]=natural.mk.map(adults.sx)
    # All maternal sample records are female-component labels in this data.
    assert natural.mother_sx.str.startswith("F").all()
    assert (len(natural),int(natural.mk.nunique()))=={2018:(158,32),2020:(199,32)}[yr]
    sensitivity=[]
    for cutoff in (.95,.975,.99,.995):
        confirmed=natural[natural.Father.notna()&natural.prob.ge(cutoff)]
        sires=confirmed[confirmed.fk.isin(h.index)]
        assert not (sires.mk==sires.fk).any()
        counts=sires.fk.value_counts()
        hh=h.copy();hh["s"]=hh.index.map(counts).fillna(0).astype(int)
        have=hh.fruits.gt(0);seen=hh.s.gt(0)
        row={"confidence_cutoff":cutoff,
           "natural_assignments_any_sire":len(confirmed),"natural_assignments_strict_H":len(sires),
           "natural_H_assignments_from_strict_F_mothers":int(sires.mother_sx.eq("F").sum()),
           "n_H_genets_with_sampled_siring":int(seen.sum()),
           "n_H_genets_with_own_fruit_and_sampled_siring":int((have&seen).sum()),
           "H_own_fruits_from_genets_without_sampled_siring":int(hh.loc[~seen,"fruits"].sum()),
           "sampled_H_siring_from_fruit_positive_genets":int(hh.loc[have,"s"].sum()),
           "sampled_H_siring_from_fruit_zero_genets":int(hh.loc[~have,"s"].sum()),
           "raw_cohort_fruit_to_assignment_ratio":float(hh.fruits.sum()/max(1,len(sires))),
           "fruit_positive_H_raw_ratio":float(hh.loc[have,"fruits"].sum()/max(1,hh.loc[have,"s"].sum()))}
        sensitivity.append(row)
    assert sensitivity[0]["natural_assignments_strict_H"]=={2018:110,2020:131}[yr]
    assert sensitivity[0]["natural_assignments_any_sire"]=={2018:132,2020:163}[yr]
    moms=natural.mk.value_counts().rename("sampled_seeds").to_frame().join(adults[["sx","fruits","flowers","valid"]])
    assert len(moms)==32
    valid=moms.valid.eq(True)
    vv=moms.loc[valid]
    return {"year":yr,"strict_H_sire_risk_and_fruit_cohort_identical":True,
       "strict_H_genets":len(h),"strict_H_natural_fruits":int(h.fruits.sum()),
       "natural_seed_records":len(natural),"natural_mothers":len(moms),
       "natural_mothers_with_strict_F_code":int(moms.sx.eq("F").sum()),
       "natural_seed_records_from_strict_F_code_mothers":int(natural.mother_sx.eq("F").sum()),
       "natural_mothers_with_valid_census_fruit_count":int(valid.sum()),
       "valid_mothers_recorded_fruit_total":int(vv.fruits.sum()),
       "valid_mothers_sampled_seed_total":int(vv.sampled_seeds.sum()),
       "valid_mothers_where_sample_seed_count_exceeds_recorded_fruits":int((vv.sampled_seeds>vv.fruits).sum()),
       "valid_mothers_where_sample_count_equals_recorded_fruits":int((vv.sampled_seeds==vv.fruits).sum()),
       "valid_mothers_where_sample_count_below_recorded_fruits":int((vv.sampled_seeds<vv.fruits).sum()),
       "no_population_paternity_expansion_from_these_counts":True,
       "confidence_sensitivity":sensitivity}
def main():
    f=source()
    return {"status":"v67_original_natural_sire_cohort_and_sampling_frame_audit",
       "source_doi":"10.1002/ajb2.16239","zip_sha256":ZIP_SHA,
       "by_year":[year(f,y) for y in (2018,2020)],
       "limits":["same source 2018/2020 study","no new paternity inference or truth","no offspring raw SSR",
          "no known external siring capture fraction or relative recruitment","female-component maternal genets only",
          "zero sampled sire count not zero true paternity","mother sample counts vs census do not establish sampling fraction",
          "no original individual data exported"]}
if __name__=="__main__":
    OUT.mkdir(exist_ok=True);x=main()
    (OUT/"daphne_v67_natural_source_audit.json").write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print("DAPHNE_V67_SOURCE_AUDIT",json.dumps(x,allow_nan=False),flush=True)
