#!/usr/bin/env python3
"""Primary-row multi-stage reproductive audit (Brassica and Veratrum).

Source data retrieved exclusively through publicly documented Dryad version ZIP
download. Published 2025/2026 study values are kept distinct from new
exploratory re-estimation. Exports have aggregated counts and model terms,
never original plant/flower IDs, coordinates, or observational rows.

Brassica: plant (n=86), hand-pollinated flower fate (~2201) and ovule fate
(~9421) are DIFFERENT levels. Experiments 1 & 2 are distinct, cannot multiply
rates into a whole-plant fitness estimate.
Veratrum: processed pair records compare natural vs supplemented within three
manipulation extents. Whole-plant comparison ranks distinct plants by size;
cannot call it within-genet paired experimental treatment.
"""
from __future__ import annotations
import hashlib, io, json, math, urllib.request, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

OUT=Path("analysis_results_primary_stage_refit")
OUT.mkdir(exist_ok=True)
DATA={
 "Brassica":{"id":434166,"doi":"10.5061/dryad.hx3ffbgtq",
 "sha":"6df15f56f904318dcf236c157003ce37cccade495bb9bc6fac108edcbd6a175e"},
 "Veratrum":{"id":405186,"doi":"10.5061/dryad.9s4mw6mw2",
 "sha":"3a69b5c3e742ebb139769175280bc6f9bb27ec90f18165523d4c052bff32da1c"}
}
GEN=np.random.default_rng(20261115)
N_BOOT=4000

def obtain(name):
    spec=DATA[name]
    url=f"https://datadryad.org/api/v2/versions/{spec['id']}/download"
    req=urllib.request.Request(url,headers={"User-Agent":"Published-data-reproduction/0.1","Accept":"application/zip"})
    with urllib.request.urlopen(req,timeout=70) as r:
        raw=r.read(850000)
    # Dryad may generate ZIP with variable envelope timestamps. Instead
    # authenticate each table's uncompressed original bytes against the pinned
    # SHA256 from our previously inspected official original ZIP inventory.
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("not a published Dryad ZIP")
    expected={
      "Brassica":{
         "1_Table_1_-_Vegetative_mass.csv":"eee389756b1b713dabe3bb9d244a9a1096066ac8df63840496c12794370c932d",
         "3_Table_2_-_Fruit_set_flower_level_4_raceme.csv":"a6fa8c969aa9e76e12081979a8a1c431c2ba183a8bcb363ac9adda165c61a756",
         "6_Table_3_-_Ovule_fates_flower_level.csv":"4dc4718c5d2233fc95ebddbe799653e88bc60358d4b601cb9f618268aaa2d45d",
         "7_Table_3_-_Ovule_fates_ovule_level.csv":"2ccafcb2c3ab73ec2c882c60f4047cbb34315287fcdc296388773dafaeba076a"
      },
      "Veratrum":{"Fruit_Seed_Data.csv":"d020dfa7fbdbe05ab822c64004b9de83b9840e07d02d04918678974f34848b7a"}
    }
    result={}
    verified=set()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for item in z.infolist():
            basename=Path(item.filename).name
            if not basename.lower().endswith(".csv"):continue
            if item.file_size>400000:raise ValueError("unexpected oversized data member")
            blob=z.read(item)
            if basename in expected[name]:
                if hashlib.sha256(blob).hexdigest()!=expected[name][basename]:
                    raise ValueError("published original CSV member SHA mismatch "+basename)
                verified.add(basename)
            # Legacy Mac CR line endings are normalized rather than reinterpreted.
            txt=blob.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
            df=pd.read_csv(io.StringIO(txt),sep=None,engine="python",na_values=["NA","NaN"],keep_default_na=True)
            df.columns=[str(x).strip() for x in df.columns]
            result[basename]=df
    if verified!=set(expected[name]):
        raise ValueError("missing required published original files "+str(set(expected[name])-verified))
    print("SOURCE_MEMBERS_SHA256_VERIFIED",name,len(verified),"ZIP-envelope SHA mismatch is acceptable",flush=True)
    return result

def numeric(df,keys):
    for c in keys:df[c]=pd.to_numeric(df[c],errors="coerce")
    return df

def native(x):
    if x is None:return None
    if isinstance(x,(float,np.floating)):
        return round(float(x),6) if math.isfinite(x) else None
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,dict):return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(y) for y in x]
    return x

def coeff(model,terms=None):
    p=model.params; ci=model.conf_int()
    return {str(k):{"beta":native(p[k]),"ci_low":native(ci.loc[k,0]),"ci_high":native(ci.loc[k,1]),
            "p":native(model.pvalues[k])} for k in (terms or p.index)}

def bootstrap_plant_frac(df,rate_col,group_col):
    data={}
    for key,x in df.groupby(group_col):
        arr=x[rate_col].to_numpy(float)
        b=GEN.choice(arr,size=(N_BOOT,len(arr)),replace=True).mean(axis=1)
        data[str(key)]={"n_plants":len(arr),"plant_mean":native(arr.mean()),
                         "bootstrap_95":[native(v) for v in np.quantile(b,[.025,.975])]}
    if len(data)==2:
        labels=sorted(data)
        g1=df.loc[df[group_col]==int(labels[0]),rate_col].to_numpy(float)
        g2=df.loc[df[group_col]==int(labels[1]),rate_col].to_numpy(float)
        b1=GEN.choice(g1,size=(N_BOOT,len(g1)),replace=True).mean(axis=1)
        b2=GEN.choice(g2,size=(N_BOOT,len(g2)),replace=True).mean(axis=1)
        data["between_treatments"]={"75_minus_25":native(g2.mean()-g1.mean()),
                  "cluster_bootstrap_95":[native(v) for v in np.quantile(b2-b1,[.025,.975])]}
    return data

def brassica(files):
    df=files["1_Table_1_-_Vegetative_mass.csv"].copy()
    required=["Treatment","Plant","Raceme","nflowers","nfruit","intfruit","Totveg"]
    if any(c not in df for c in required):raise ValueError("missing Brassica raceme columns")
    numeric(df,required)
    original_raceme_rows=len(df)
    missing_by_column={k:int(v) for k,v in df[required].isna().sum().items()}
    print("BRASSICA_RACEME_MISSING_COUNTS",missing_by_column,flush=True)
    # Do not silently turn missing mass/fruit counts into zero. Repeated total
    # plant vegetative mass may be missing in one raceme but present in another.
    if df.groupby("Plant").Treatment.nunique().max()!=1:raise ValueError("plant treatment conflict")
    if df.groupby("Plant").Totveg.nunique().max()>1:raise ValueError("mass inconsistent within plant")
    df["Totveg"]=df.groupby("Plant")["Totveg"].transform(lambda v:v.dropna().iloc[0] if v.notna().any() else float("nan"))
    core=["Treatment","Plant","nflowers","nfruit","Totveg"]
    invalid_plants=set(df.loc[df[core].isna().any(axis=1),"Plant"].dropna())
    df=df[~df.Plant.isin(invalid_plants)].copy()
    if len(df)==0 or df.Plant.nunique()<65:raise ValueError("Too much missing original raceme data")
    plant=df.groupby(["Plant","Treatment"],as_index=False).agg(
        nflowers=("nflowers","sum"),nfruit=("nfruit","sum"),intfruit=("intfruit","sum"),
        Totveg=("Totveg","first"),n_racemes=("Raceme","size"))
    if len(plant)<65:raise ValueError(f"too few complete plants {len(plant)}")
    if (plant.nfruit>plant.nflowers).any():raise ValueError("fruit counts exceed flower numbers")
    if (plant.Totveg<=0).any() or (plant.nflowers<=0).any():raise ValueError("nonpositive mass/flowers")
    plant["fruit_rate"]=plant.nfruit/plant.nflowers
    plant["log_vegetative_mass"]=np.log(plant.Totveg)
    summary={"plant_total":len(plant),"raceme_rows":len(df),
             "source_raceme_rows":original_raceme_rows,
             "missing_by_raceme_column":missing_by_column,
             "excluded_plants_with_missing_required_fields":len(invalid_plants),
             "treatment_n_plants":native(plant.Treatment.value_counts().sort_index().to_dict()),
             "all_fruits":native(plant.nfruit.sum()),"all_flowers":native(plant.nflowers.sum()),
             "by_treatment":{}}
    for t,group in plant.groupby("Treatment"):
        summary["by_treatment"][str(int(t))]={
            "plants":len(group),"total_flowers":native(group.nflowers.sum()),
            "total_fruits":native(group.nfruit.sum()),"pooled_fraction":native(group.nfruit.sum()/group.nflowers.sum()),
            "mean_plant_fruit_fraction":native(group.fruit_rate.mean()),
            "mean_plant_flowers":native(group.nflowers.mean()),
            "mean_plant_fruits":native(group.nfruit.mean()),
            "median_vegetative_mass":native(group.Totveg.median())}
    # Same per-genet independent units as Experiment 1: plant.
    # Binomial response, adjusted for assigned pollen supplementation and mass.
    X=pd.get_dummies(plant.Treatment.astype(int).astype(str),prefix="pollination",drop_first=True)
    X["lnVeg"]=plant.log_vegetative_mass
    X=sm.add_constant(X.astype(float))
    y=np.column_stack((plant.nfruit.to_numpy(float),(plant.nflowers-plant.nfruit).to_numpy(float)))
    fit=sm.GLM(y,X,family=sm.families.Binomial()).fit(cov_type="HC3")
    summary["plant_binomial_pollen_and_mass"]={
        "terms":coeff(fit),"n_independent_plant_units":len(plant),
        "mass_not_randomized":True,"overdispersion_pearson":native(fit.pearson_chi2/fit.df_resid)}
    if "lnVeg" in fit.params.index:
        b=float(fit.params["lnVeg"])
        lo,hi=fit.conf_int().loc["lnVeg"]
        summary["plant_binomial_mass_OR_per_doubling"]={"OR":native(math.exp(b*math.log(2))),
          "ci":[native(math.exp(lo*math.log(2))),native(math.exp(hi*math.log(2)))],
          "not_causal_vegetative_mass_effect":True}
    # Flower records include only intentionally hand-supplemented flowers in Exp 1.
    f=files["3_Table_2_-_Fruit_set_flower_level_4_raceme.csv"].copy()
    rcols=["Treatment","Plant","Totveg","proppoll","Raceme","Flw","fday","Frt"]
    numeric(f,rcols)
    f=f.dropna(subset=rcols).copy()
    f=f[f.Totveg>0]
    f["logVeg"]=np.log(f.Totveg)
    if f.Frt.isin([0,1]).all() is False:raise ValueError("nonbinomial floral outcome")
    if len(f)<2000 or f.Plant.nunique()!=86:raise ValueError("unexpected flower inventory")
    mods=smf.glm("Frt ~ fday*C(Treatment) + logVeg",data=f,
             family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":f.Plant})
    summary["hand_pollen_flower_timing"]={"n_flower_records":len(f),"n_distinct_plant_clusters":int(f.Plant.nunique()),
             "model":"logit Fruit ~ fday * treatment_fraction + log_vegetative_mass (cluster SE by plant)",
             "effect_terms":coeff(mods),
             "by_treatment":{},"only_hand_pollinated_flowers":True,
             "anthesis_day_not_randomly_assigned":True}
    for t,g in f.groupby("Treatment"):
        common=[]
        for p,gg in g.groupby("Plant"):
            median=gg.fday.median()
            before=gg[gg.fday<=median]
            after=gg[gg.fday>median]
            if len(after)>0 and len(before)>0:
                common.append((before.Frt.mean(),after.Frt.mean()))
        paired=np.asarray(common)
        summary["hand_pollen_flower_timing"]["by_treatment"][str(int(t))]={
            "n_records":len(g),"n_plants":int(g.Plant.nunique()),
            "overall_hand_pollen_fruit_fraction":native(g.Frt.mean()),
            "day_10pct":native(g.fday.quantile(.1)),
            "day_90pct":native(g.fday.quantile(.9)),
            "n_plant_early_late_pairs":len(common),
            "plant_mean_early_fruit_rate":native(paired[:,0].mean()) if len(common) else None,
            "plant_mean_late_fruit_rate":native(paired[:,1].mean()) if len(common) else None,
            "plant_mean_late_minus_early":native((paired[:,1]-paired[:,0]).mean()) if len(common) else None}
    # Exp2: conditional on sampled fruits and plant genets.
    o=files["7_Table_3_-_Ovule_fates_ovule_level.csv"].copy()
    cols=["Treatment","Plant","fruitN","rPosition","n_totovule","ovule_basal","fertilized","develop","abort"]
    numeric(o,cols)
    o=o.dropna(subset=cols)
    if len(o)!=9421 or o.Plant.nunique()!=46:raise ValueError("unexpected ovule source grain")
    if not (o.fertilized==o.develop+o.abort).all():
        raise ValueError("ovule fate contradiction - check input source")
    if (o.develop>o.fertilized).any():raise ValueError("development before fertilization")
    og=o.groupby(["Plant","Treatment"],as_index=False).agg(ovules=("fertilized","size"),
                      n_fert=("fertilized","sum"),n_dev=("develop","sum"),n_abort=("abort","sum"))
    og["fert_rate"]=og.n_fert/og.ovules
    og["dev_conditional"]=og.n_dev/og.n_fert
    summary["sampled_developing_fruits_ovule_fate"]={
      "n_ovule_records":len(o),"n_plants":o.Plant.nunique(),
      "fertilized_pooled":native(o.fertilized.sum()/len(o)),
      "dev_given_fertilized_pooled":native(o.develop.sum()/o.fertilized.sum()),
      "n_fert":native(o.fertilized.sum()),"n_developed_seeds":native(o.develop.sum()),
      "n_unfert":native((o.fertilized==0).sum()),
      "n_abort_after_fertilization":native(o.abort.sum()),
      "per_plant_fertilization_rate_by_treatment":bootstrap_plant_frac(og,"fert_rate","Treatment"),
      "per_plant_development_conditional_by_treatment":bootstrap_plant_frac(og,"dev_conditional","Treatment"),
      "conditional_on_sampled_developing_fruits_not_all_flowers":True,
      "plant_resampling_not_ovule_independence":True}
    sf=files["6_Table_3_-_Ovule_fates_flower_level.csv"].copy()
    numeric(sf,["Treatment","Plant","n_totovule","n.unfert","n.abort","n.seed"])
    if len(sf)!=719:raise ValueError("expected 719 source fruit rows")
    consistent=(sf["n_totovule"]==sf["n.unfert"]+sf["n.abort"]+sf["n.seed"])
    summary["sampled_developing_fruits_ovule_fate"]["n_fruit_rows"]=len(sf)
    summary["sampled_developing_fruits_ovule_fate"]["n_fruits_with_valid_ovule_partition"]=native(consistent.sum())
    summary["sampled_developing_fruits_ovule_fate"]["all_ovule_partition_exact"]=bool(consistent.all())
    return summary

def log_effect(hand,natural):
    mask=np.isfinite(hand)&np.isfinite(natural)&(hand>0)&(natural>0)
    return np.log(hand[mask]/natural[mask]),int(mask.sum()),int((~mask).sum())

def veratrum(files):
    df=files["Fruit_Seed_Data.csv"].copy()
    columns=["Individual.hand","Individual.natural","Total.Herm.hand","Total.Fruit.hand","Fruiting.hand",
             "Total.Herm.natural","Total.Fruit.natural","Fruiting.natural","Level","Seeds.hand","Seeds.natural"]
    if not set(columns).issubset(df):raise ValueError("unexpected Veratrum paired-table columns")
    numeric(df,[k for k in columns if k not in ("Individual.hand","Individual.natural","Level")])
    report={"source_paired_records":len(df),"raw_individual_records_expected":129,
            "by_level":{},"must_not_pseudocount_zero_seed_pairs":True,
            "whole_plant_comparison_is_size_rank_matched_different_individuals":True}
    if len(df)!=97:raise ValueError("unexpected original preprocessed source row count")
    levels=["Individual","Inflorescence","Flower"]
    obs={}
    for group in levels:
        d=df[df.Level==group].copy()
        ha=d["Seeds.hand"].to_numpy(float);na=d["Seeds.natural"].to_numpy(float)
        logratio,n_pos,n_excl=log_effect(ha,na)
        if n_pos<8:raise ValueError("Too few positive seed pairs "+group)
        resample=GEN.choice(logratio,size=(N_BOOT,len(logratio)),replace=True).mean(axis=1)
        first=d["Fruiting.hand"].to_numpy(float);second=d["Fruiting.natural"].to_numpy(float)
        flr,n_fr,n_fr_excl=log_effect(first,second)
        frboot=GEN.choice(flr,size=(N_BOOT,len(flr)),replace=True).mean(axis=1) if len(flr) else None
        expect_same=(d["Individual.hand"].astype(str)==d["Individual.natural"].astype(str))
        report["by_level"][group]={
             "n_pair_records":len(d),"n_seed_logratio_pairs":n_pos,
             "excluded_nonpositive_or_missing_seed_pairs":n_excl,
             "mean_seed_per_fruit_hand":native(np.nanmean(ha)),
             "mean_seed_per_fruit_natural":native(np.nanmean(na)),
             "log_ratio_of_means":native(math.log(np.nanmean(ha)/np.nanmean(na))),
             "mean_of_pairwise_seed_logratios":native(logratio.mean()),
             "paired_seed_logratio_bootstrap_95":[native(v) for v in np.quantile(resample,[.025,.975])],
             "published_mean_seed_log_PL_for_context":{"Individual":.129,"Inflorescence":.453,"Flower":.787}[group],
             "mean_of_pairwise_fruit_logratios":native(flr.mean()) if len(flr) else None,
             "fruit_logratio_positive_pair_count":n_fr,
             "n_excluded_nonpositive_fruit_logratio":n_fr_excl,
             "fruit_logratio_bootstrap_95":[native(v) for v in np.quantile(frboot,[.025,.975])] if frboot is not None else None,
             "same_id_both_exposure_cols":native(expect_same.sum()),
             "same_id_rate":native(expect_same.mean()),
             "size_rank_matching_is_not_genet_paired_treatment":group=="Individual"}
        obs[group]=logratio
    flower=obs["Flower"];whole=obs["Individual"]
    b1=GEN.choice(flower,size=(N_BOOT,len(flower)),replace=True).mean(axis=1)
    b2=GEN.choice(whole,size=(N_BOOT,len(whole)),replace=True).mean(axis=1)
    report["unadjusted_between_level_difference_flower_minus_whole"]={
        "difference_in_mean_logratio":native(flower.mean()-whole.mean()),
        "between_subsample_bootstrap_95":[native(v) for v in np.quantile(b1-b2,[.025,.975])],
        "not_a_paired_within_individual_causal_contrast":True}
    return report

def main():
    report={"status":"failed","provenance":DATA,
        "interpretation":"new exploratory reanalysis of ORIGINAL public published individual/flower/ovule records",
        "no_original_individual_ids_exported":True,
        "distinct_studies_not_combined_into_meta_effect":True,
        "new_results_not_author_reanalysis_certified":True}
    brass=obtain("Brassica")
    vera=obtain("Veratrum")
    report["brassica"]=brassica(brass)
    report["veratrum"]=veratrum(vera)
    report["status"]="success_original_raw_data_reanalysis"
    outfile=OUT/"primary_experimental_stage_audit_20261115.json"
    outfile.write_text(json.dumps(native(report),indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(native(report),ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":
    main()
