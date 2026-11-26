#!/usr/bin/env python3
"""Original Silene latifolia 2004 display and 2007 survival reanalysis.

Original Dryad v11856, DOI 10.5061/dryad.20k62r83, all five source tables
verified against stable uncompressed SHA256. Original selective-breeding line
comparisons confounded by correlated leaf physiology; no causal floral-cost
effect is inferred. 'Total Flowers' is a SUM of repeated census open-flower
snapshots, not unique flower construction or total flower dry biomass.
No original individual ID exported.
"""
from __future__ import annotations
import hashlib, io, json, urllib.request, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

OUT=Path("analysis_results_silene_multiyear_original")
OUT.mkdir(exist_ok=True)
SOURCE={
"Flower counts 2004.csv":"1a8b0c49b040942d3f60a6495f74fb6fdaa54c159f5e2b776a188d24ce657594",
"flower counts 2005.csv":"a04fcddd7c74400adb53243e78277cce0c8e3d8f0aa93a35c12c9db45c40e1c1",
"Fruit and seed data.csv":"880d8180401f17c27e3f6ab97fc1efbceac19184c4d1ef62de6833ea1fa2f84b",
"Individual data.csv":"bd86a826eb2b0179db6fbf0a863c4d998a88e6fba911d4e615f0bea40ba055bf",
"Siring success and male flowering.csv":"544dceaf0007d0fa6bf4a1ec211ba2774f701719cb27c30606e52b9a0f16959f"
}
def load():
    url="https://datadryad.org/api/v2/versions/11856/download"
    req=urllib.request.Request(url,headers={"User-Agent":"Published-original-Silene-study-reanalysis/0.1"})
    with urllib.request.urlopen(req,timeout=70) as res:raw=res.read(350000)
    if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("not a public zip")
    dfs={};checked=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for info in z.infolist():
            name=Path(info.filename).name
            if name not in SOURCE:continue
            data=z.read(info)
            if hashlib.sha256(data).hexdigest()!=SOURCE[name]:raise ValueError("source member hash mismatch "+name)
            txt=data.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
            frame=pd.read_csv(io.StringIO(txt),sep=None,engine="python",na_values=["NA","N/A"],keep_default_na=True)
            frame.columns=[c.strip() for c in frame.columns]
            dfs[name]=frame
            checked.append(name)
    if set(checked)!=set(SOURCE):raise ValueError("missing one official source CSV")
    print("SOURCE_5_MEMBERS_SHA_OK",sorted(checked),flush=True)
    return dfs

def num(x):
    return pd.to_numeric(x,errors="coerce")

def scalar(x):
    if x is None:return None
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return round(float(x),6) if np.isfinite(x) else None
    if isinstance(x,dict):return {str(k):scalar(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [scalar(v) for v in x]
    return x

def categories(s):
    return {str(k):int(v) for k,v in s.astype("string").fillna("<NA>").value_counts(dropna=False).items()}

def coefficient(fit):
    c=fit.conf_int()
    return {str(k):{"b":scalar(fit.params[k]),"lo":scalar(c.loc[k,0]),"hi":scalar(c.loc[k,1]),"p":scalar(fit.pvalues[k])} for k in fit.params.index}

def main():
    a=load()
    ind=a["Individual data.csv"].copy()
    c04=a["Flower counts 2004.csv"].copy()
    c05=a["flower counts 2005.csv"].copy()
    fruit=a["Fruit and seed data.csv"].copy()
    pater=a["Siring success and male flowering.csv"].copy()
    if len(ind)!=208 or len(c04)!=208 or len(c05)!=128 or len(fruit)!=104 or len(pater)!=275:
        raise ValueError("unexpected original row counts")
    for frame,key in ((ind,"Plant"),(c04,"Plant"),(c05,"plant"),(fruit,"Plant"),(pater,"Plant")):
        frame["_id"]=frame[key].astype("string").str.strip()
    for d in (ind,c04,c05):
        d["Sex"]=d["Sex"].astype("string").str.lower().str.strip()
        d["Line"]=d["Line"].astype("string").str.lower().str.strip()
    if any(d["_id"].duplicated().any() for d in (ind,c04,c05,fruit)):
        raise ValueError("source dataset has duplicate plant IDs")
    merged=ind.merge(c04[["_id","Sex","Line","Total Flowers","N counts flowering"]],
        on="_id",how="left",validate="one_to_one",suffixes=("","_2004"),indicator=True)
    matches=(merged["_merge"]=="both").sum()
    bad_morph=(merged["Sex"]!=merged["Sex_2004"])|(merged["Line"]!=merged["Line_2004"])
    if bool(bad_morph.fillna(False).any()):raise ValueError("sex line mismatch source datasets")
    summary={
      "status":"success_original_rows_loaded",
      "study":"Delph and Herlihy Evolution online 2011 / print 2012",
      "doi":"10.1111/j.1558-5646.2011.01510.x",
      "dryad":"10.5061/dryad.20k62r83",
      "source_version":11856,
      "n_individuals":len(ind),
      "n_2004_plant_snapshot_rows":len(c04),
      "n_2005_plant_snapshot_rows":len(c05),
      "n_female_fruit_seed_source_rows":len(fruit),
      "n_paternal_line_weekly_rows":len(pater),
      "match_to_2004_counts":int(matches),
      "match_to_2005_counts":int(c05["_id"].isin(ind["_id"]).sum()),
      "fruit_source_IDs_matching_initial":int(fruit["_id"].isin(ind["_id"]).sum()),
      "sex_counts":categories(ind.Sex),
      "line_counts":categories(ind.Line),
      "repeated_2004_open_flower_snapshot_sum_not_unique_blooms":True,
      "paternity_assigned_to_selection_lines_not_individual_fathers":True,
      "selection_line_and_leaf_physiology_confounding":True,
      "original_phenotype_columns":list(ind.columns.drop("_id")),
      "alive_status_raw_codes":{}
    }
    ycols=[x for x in ind.columns if x.startswith("Alive ")]
    for key in ycols:summary["alive_status_raw_codes"][key]=categories(ind[key])
    merged["display2004"]=num(merged["Total Flowers"])
    merged["num_census_positive2004"]=num(merged["N counts flowering"])
    summary["display_2004_missing"]=int(merged.display2004.isna().sum())
    summary["display_2004_by_sex_line"]={}
    for (sex,line),g in merged.groupby(["Sex","Line"],dropna=True):
        summary["display_2004_by_sex_line"][str(sex)+"_"+str(line)]={
            "n":len(g),"median_open_census_sum":scalar(g.display2004.median()),
            "mean_open_census_sum":scalar(g.display2004.mean()),
            "n_total_flower_values_positive":int((g.display2004>0).sum())}
    raw_alive=ind[ycols[-1]]
    str_alive=raw_alive.astype("string").str.strip().str.lower()
    mapping={"1":1,"0":0,"true":1,"false":0,"yes":1,"no":0,"y":1,"n":0,"alive":1,"dead":0}
    merged["alive_end"]=str_alive.map(mapping)
    valid=merged.dropna(subset=["alive_end","display2004","Sex","Line"]).copy()
    summary["last_survival_endpoint"]=ycols[-1]
    summary["survival_n_binary_usable"]=len(valid)
    summary["survival_missing_end"]=int(merged.alive_end.isna().sum())
    summary["end_missing_by_early_2005_alive_status"]=categories(merged.loc[merged.alive_end.isna(),"Alive 4/2005"])
    summary["survival_by_sex_line"]=(
        {str(sex)+"_"+str(line):{"n_with_recorded_end":len(g),
                                 "alive_end":int((g.alive_end==1).sum()),
                                 "alive_fraction":scalar(g.alive_end.mean())}
         for (sex,line),g in valid.groupby(["Sex","Line"])})
    summary["survival_dead_count"]=int((valid.alive_end==0).sum())
    summary["survival_alive_count"]=int((valid.alive_end==1).sum())
    if len(valid)>=160 and valid.alive_end.nunique()==2:
        valid["log1p_display2004"]=np.log1p(valid.display2004.clip(lower=0))
        mod=smf.glm("alive_end ~ log1p_display2004 + C(Sex)*C(Line)",data=valid,
                    family=sm.families.Binomial()).fit(cov_type="HC3")
        summary["survival_observational_model"]={
            "n_independent_genets":len(valid),"formula":"logit alive_2007 ~ ln(1+sum_open_flower_snapshots_2004) + sex*artificial_selection_line",
            "coefficients":coefficient(mod),"overdispersion":scalar(mod.pearson_chi2/mod.df_resid),
            "not_randomized_flower_cost":True,"leaf_physiology_not_pre_treatment_controlled":True}
    else:
        summary["survival_observational_model"]={"status":"not_fitted","reason":"insufficient binary coded outcome or 2004-linked observations"}
    # Survivorship and missing-censoring sensitivity. Binary future status
    # is Y/N in original source; missing is never silently recoded as dead.
    summary["survival_temporal_endpoints"]={}
    for endpoint in ycols:
        v=merged[endpoint].astype("string").str.strip().str.lower().map(mapping)
        dat=merged.copy()
        dat["alive_bin"]=v
        dat=dat.dropna(subset=["alive_bin","display2004","Sex","Line"]).copy()
        if dat.alive_bin.nunique()<2 or len(dat)<100:continue
        dat["log1p_display2004"]=np.log1p(dat.display2004.clip(lower=0))
        m=smf.glm("alive_bin ~ log1p_display2004 + C(Sex)*C(Line)",
                  data=dat,family=sm.families.Binomial()).fit(cov_type="HC3")
        summary["survival_temporal_endpoints"][endpoint]={
            "n":len(dat),"n_alive":int(dat.alive_bin.sum()),
            "display_beta":coefficient(m).get("log1p_display2004"),
            "timepoint_not_inferred_to_be_lifetime_fitness":True}
    summary["final_survival_sensitivity_missing_extremes"]={}
    for missing_status in (0,1):
        dat=merged.dropna(subset=["display2004","Sex","Line"]).copy()
        dat["alive_bin"]=dat["alive_end"].fillna(missing_status)
        dat["log1p_display2004"]=np.log1p(dat.display2004.clip(lower=0))
        m=smf.glm("alive_bin ~ log1p_display2004 + C(Sex)*C(Line)",
            data=dat,family=sm.families.Binomial()).fit(cov_type="HC3")
        summary["final_survival_sensitivity_missing_extremes"][f"all_16_missing_as_{missing_status}"]={
            "n":len(dat),"n_alive":int(dat.alive_bin.sum()),
            "display_beta":coefficient(m).get("log1p_display2004"),
            "hypothetical_extreme_only_not_a_truth_claim":True}
    summary["sex_specific_final_survival"]={}
    for sex in ("f","m"):
        dd=valid[valid.Sex==sex].copy()
        dd["log1p_display2004"]=np.log1p(dd.display2004.clip(lower=0))
        m=smf.glm("alive_end ~ log1p_display2004 + C(Line)",
                   data=dd,family=sm.families.Binomial()).fit(cov_type="HC3")
        summary["sex_specific_final_survival"][sex]={
           "n":len(dd),"n_alive":int(dd.alive_end.sum()),
           "display_beta":coefficient(m).get("log1p_display2004"),
           "few_surviving_male_events_and_low_power":sex=="m"}
    # Female year not inferred from file title without original protocol.
    fruit["Number of Fruits"]=num(fruit["Number of Fruits"])
    fruit["Mass of all seeds (mg)"]=num(fruit["Mass of all seeds (mg)"])
    female=fruit.merge(ind[["_id","Sex","Line"]],on="_id",how="left",validate="one_to_one")
    summary["female_fruit_table_sex_counts"]=categories(female.Sex)
    summary["female_fruit_source_rows_excluded_as_male"]=int((female.Sex!="f").sum())
    female=female[female.Sex=="f"].copy()
    summary["female_seed_mass_by_line"]={}
    for line,g in female.groupby("Line"):
        summary["female_seed_mass_by_line"][str(line)]={
           "n":len(g),"n_nonmissing_seed_mass":int(g["Mass of all seeds (mg)"].notna().sum()),
           "median_seed_mass_mg":scalar(g["Mass of all seeds (mg)"].median()),
           "median_fruits":scalar(g["Number of Fruits"].median())}
    # Separate current seed-mass association with a 2004 display *snapshot*.
    # Fruit/seed source year is not assumed; not a mechanistic allocation test.
    fjoined=female.merge(c04[["_id","Total Flowers"]],on="_id",how="inner",validate="one_to_one")
    fjoined["seed_mass"]=num(fjoined["Mass of all seeds (mg)"])
    fjoined["display"]=num(fjoined["Total Flowers"])
    fjoined=fjoined.dropna(subset=["seed_mass","display","Line"]).copy()
    if (fjoined.seed_mass<0).any() or (fjoined.display<0).any():
        raise ValueError("negative original source seed mass or flower snapshot")
    if len(fjoined)>=80:
        fjoined["ln_seeds"]=np.log1p(fjoined.seed_mass)
        fjoined["ln_display"]=np.log1p(fjoined.display)
        m=smf.ols("ln_seeds ~ ln_display + C(Line)",data=fjoined).fit(cov_type="HC3")
        summary["female_seed_mass_vs_2004_flower_snapshot"]={
           "n_independent_genets":len(fjoined),
           "model":"ln1p seed mass ~ ln1p 2004 repeated open-flower snapshots + selection line",
           "display_beta":coefficient(m).get("ln_display"),
           "seed_sample_year_not_independently_verified":True,
           "trait_covariance_not_randomized":True}
    summary["reproduction_of_female_plants_not_lifetime_offspring"]=True
    outfile=OUT/"silene_original_multiyear_survival_source_audit.json"
    outfile.write_text(json.dumps(scalar(summary),ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(scalar(summary),ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":
    main()
