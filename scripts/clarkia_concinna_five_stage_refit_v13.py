#!/usr/bin/env python3
"""Original Clarkia concinna five-study floral-investment conversion audit.

Miller & Kay, "Large floral displays increase self-pollination but do not affect
visitation rates in Clarkia concinna", Madrono 2025, DOI 10.3120/0024-9637-72.1.a1;
author original Dryad 2024 doi:10.5061/dryad.8pk0p2nz2.

Source provenance from prior successful audit: uncompressed primary CSV SHA256
verified for five separate experiments; never publish original observation IDs,
site coordinates, treatment plant ID or raw source tables.

Study 1: consecutive real flower probes by visiting insect taxon, observational
Study 2: manipulated true male flower number, self versus neighbor pollen-dye
        deposition on one focal stigma, source dye rank is ordinal SURROGATE.
Study 3: manipulated flower number, one focal female-phase seed and seed weight,
        many missing fruits; NOT whole-plant lifetime seed or actual dry flower mass.
Study 4: neighborhood patch size vs 20-flower/20-minute visit counts observational
Study 5: factorial manipulated real floral signal and available male pollen
        (1, 6 emasculated, 6 intact), patch size (small vs large),
        focal emasculated stigma pollen grains / deposition time.

All stage comparisons are separately labeled observational vs manipulated.
Patch size is partly natural, partly thinned, NOT randomly assigned between
geographic sites. No functional self-incompatibility or actual paternity measured.
"""
from __future__ import annotations
import hashlib,io,json,math,re,urllib.parse,urllib.request,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr,kruskal

OUT=Path("analysis_results_clarkia_concinna_stages_v13");OUT.mkdir(exist_ok=True)
DOI="10.5061/dryad.8pk0p2nz2"
FILENAME="Floral_display_and_geitonogamy_data_study_{}.csv"
SHA={
  1:"90baf97cdfe3645fd6b99e95c2838e1231b84539ccc5190ac626e9f3ef3af5fa",
  2:"f95990d23f06f5a7c92929b247679ce2a2e2bd8c1aad6ed3545a883723ebde09",
  3:"b8143f9998a284b60733cfbb6354a290aff5b3bfe35b7fe0f50eb14f0ea33e38",
  4:"b78f62f3ffb314db95bc2be8cac24bcadadb42d9883a607c188e86b37bfcaae4",
  5:"7aea28b35e02e4bbd050ab9fa5b11bcb2f8e595d4c3fb8c76d3e9d5c178a698b"}
SEED=20261116;RNG=np.random.default_rng(SEED);N_BOOT=3500

def scalar(v):
    if isinstance(v,(bool,np.bool_)):return bool(v)
    if isinstance(v,(np.integer,int)):return int(v)
    if isinstance(v,(np.floating,float)):return round(float(v),6) if np.isfinite(v) else None
    if isinstance(v,(list,tuple,np.ndarray)):return [scalar(x) for x in v]
    if isinstance(v,dict):return {str(k):scalar(x) for k,x in v.items()}
    return v
def fetch_dataset():
    base="https://datadryad.org"
    doi_url=base+"/api/v2/datasets/"+urllib.parse.quote("doi:"+DOI,safe="")
    req=urllib.request.Request(doi_url+"/versions",headers={"User-Agent":"Clarkia-conversion-primary-reanalysis-v13","Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=65) as r:vers=json.load(r)
    vals=[]
    for a in vers.get("_embedded",{}).get("stash:versions",[]):
        link=a.get("_links",{}).get("self",{}).get("href","")
        m=re.search(r"/versions/(\d+)",link)
        if m:vals.append((int(a.get("versionNumber",0)),int(m.group(1))))
    if not vals:raise ValueError("no official public Dryad version")
    version=max(vals)[1]
    req=urllib.request.Request(base+f"/api/v2/versions/{version}/download",
        headers={"User-Agent":"Clarkia-conversion-primary-reanalysis-v13","Accept":"application/zip"})
    with urllib.request.urlopen(req,timeout=75) as r:raw=r.read(850000)
    if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("published source ZIP absent")
    frames={}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for k in SHA:
            ff=[info for info in z.infolist() if Path(info.filename).name==FILENAME.format(k)]
            if len(ff)!=1:raise ValueError("missing original study "+str(k))
            b=z.read(ff[0])
            if hashlib.sha256(b).hexdigest()!=SHA[k]:raise ValueError("original source SHA mismatch "+str(k))
            a=pd.read_csv(io.BytesIO(b))
            a.columns=[str(x).strip() for x in a.columns]
            frames[k]=a
    return version,frames
def fmtfit(m):
    ci=m.conf_int()
    return {str(k):{"beta":scalar(m.params[k]),
              "CI_95":scalar(ci.loc[k].to_numpy()),"p":scalar(m.pvalues[k])}
            for k in m.params.index if k=="ln_display" or "treatment" in k or "patch" in k or k=="ln_patch" or k=="st_dia" or k=="ln_male"}
def model_contrast(m,hypothesis):
    res=m.t_test(hypothesis)
    z=float(np.squeeze(res.effect));sd=float(np.squeeze(res.sd))
    return {"estimate_log1p":scalar(z),
      "CI_95":scalar([z-1.96*sd,z+1.96*sd]),
      "p":scalar(float(np.squeeze(res.pvalue))),
      "ratio_of_1plus_grains":scalar(math.exp(z)),
      "CI_95_ratio":scalar([math.exp(z-1.96*sd),math.exp(z+1.96*sd)])}

def descriptive_group(d,groupby,response):
    out={}
    for lab,g in d.groupby(groupby,dropna=False):
        x=pd.to_numeric(g[response],errors="coerce")
        d2=x.dropna()
        out[str(lab)]={"n_source_rows":len(g),"n_measured":int(d2.size),
          "n_missing":int(x.isna().sum()),"mean":scalar(d2.mean()),"median":scalar(d2.median()),
          "n_zero":int((d2==0).sum())}
    return out

def study1(d):
    x=d.copy()
    for k in ("fl visits","fl total","ind"):x[k]=pd.to_numeric(x[k],errors="coerce")
    x=x.dropna(subset=["fl visits","fl total","ind","genus","date"]).copy()
    if len(x)!=1626 or x["fl visits"].min()<1 or x["fl total"].min()<1:
        raise ValueError("source study1 malformed")
    x["ln_display"]=np.log(x["fl total"])
    x["additional"]=np.maximum(0,x["fl visits"]-1)
    x["tag"]=x.genus.astype(str)+"::"+x["ind"].astype(str)
    x["genus"]=x.genus.astype(str)
    x["date"]=x.date.astype(str)
    # Post-source pattern: correlations are not plant-fitness or unique genet.
    base=smf.glm("additional ~ ln_display + C(genus) + C(date)",data=x,
            family=sm.families.NegativeBinomial(alpha=1)).fit(
            cov_type="cluster",cov_kwds={"groups":x.tag})
    co=base.conf_int().loc["ln_display"]
    taxa=x.genus.value_counts()
    major=taxa[taxa>=50].index.tolist()
    z=x[x.genus.isin(major)].copy()
    m=smf.glm("additional ~ ln_display*C(genus) + C(date)",data=z,
        family=sm.families.NegativeBinomial(alpha=1)).fit(
        cov_type="cluster",cov_kwds={"groups":z.tag})
    return {"source_visiting_bouts":len(x),"unique_insect_tag_groups_genus_number":int(x.tag.nunique()),
      "n_visits_exceed_number_open_flowers":int((x["fl visits"]>x["fl total"]).sum()),
      "n_source_one_open_flower":int((x["fl total"]==1).sum()),
      "n_with_multiple_within_plant_probes":int((x["fl visits"]>1).sum()),
      "n_qualifying_major_genera":len(major),
      "major_taxon_record_counts":{str(k):int(v) for k,v in taxa[taxa>=50].items()},
      "overall_NB_log_display_to_additional_probes":{
        "beta":scalar(base.params["ln_display"]), "CI_95":scalar(co.to_numpy()),
        "p":scalar(base.pvalues["ln_display"]),
        "IRR_per_display_doubling":scalar(2**base.params["ln_display"]),
        "IRR_95":scalar(2**co.to_numpy())},
      "genus_moderation_exploratory_model_terms":fmtfit(m),
      "within_bout_intraplant_visits_not_genet_identity_genetic_selfing":True,
      "same_individual_insect_tags_may_recur_across_bouts":True,
      "experimental_increase_in_flower_count_not_assigned_here":True}
def study2(d):
    x=d.copy()
    for k in ("# flrs on focal plant","target dye rank","surr dye rank"):
        x[k]=pd.to_numeric(x[k],errors="coerce")
    x=x.dropna(subset=["# flrs on focal plant","target dye rank","surr dye rank"]).copy()
    x=x[x["# flrs on focal plant"].isin([1,2,4,6])]
    x["ln_male"]=np.log(x["# flrs on focal plant"])
    x["self_any"]=(x["target dye rank"]>0).astype(int)
    x["neighbor_any"]=(x["surr dye rank"]>0).astype(int)
    x["date"]=x["date"].astype(str)
    if len(x)!=136:raise ValueError("unexpected source pigment n=136")
    # Robust HC3 inference; source manipulates male flowers 1,2,4,6.
    mself=smf.glm("self_any ~ ln_male",data=x,family=sm.families.Binomial()).fit(cov_type="HC3")
    mnear=smf.glm("neighbor_any ~ ln_male",data=x,family=sm.families.Binomial()).fit(cov_type="HC3")
    return {"source_manipulated_plant_stigmas":len(x),
       "source_male_flower_groups":{str(int(k)):{**descriptive_group(g,[],"target dye rank")} if False else {
          "n":len(g),"self_any_n":int(g.self_any.sum()),"self_any_fraction":scalar(g.self_any.mean()),
          "neighbor_any_fraction":scalar(g.neighbor_any.mean()),
          "self_rank_mean":scalar(g["target dye rank"].mean()),
          "neighbor_rank_mean":scalar(g["surr dye rank"].mean())}
          for k,g in x.groupby("# flrs on focal plant")},
       "self_pigment_arrival_vs_male_display":fmtfit(mself),
       "neighbor_pigment_arrival_vs_male_display":fmtfit(mnear),
       "self_rank_spearman":scalar(spearmanr(x["# flrs on focal plant"],x["target dye rank"]).statistic),
       "neighbor_rank_spearman":scalar(spearmanr(x["# flrs on focal plant"],x["surr dye rank"]).statistic),
       "pigment_dye_is_not_viable_self_pollen_or_offspring":True,
       "one_focal_female_stigma_per_plant":True}
def study3(d):
    x=d.copy()
    for k in ("# flrs","st dia","fl l","fl sta","fl sty","fl 50cm","stem pos",
              "# fruits (fertilized)","seed #","unfert seeds","seed weight"):
        x[k]=pd.to_numeric(x[k],errors="coerce")
    if len(x)!=130:raise ValueError("unexpected study3 original")
    x["male_remaining"]=x["# flrs"]-1
    x["ln_male"]=np.log1p(x.male_remaining)
    x["has_seed_data"]=x["seed #"].notna().astype(int)
    x["date"]=x["date established"].astype(str)
    clean=x.dropna(subset=["seed #","unfert seeds","stem pos","st dia"]).copy()
    if (clean["seed #"]<0).any() or (clean["unfert seeds"]<0).any():
        raise ValueError("negative source seeds")
    clean["seed_total"]=clean["seed #"]+clean["unfert seeds"]
    clean=clean[clean.seed_total>0].copy()
    if len(clean)<90:raise ValueError("insufficient matched source outcomes")
    # Use flower as experimental unit; single focal female flower sampled
    # per independently treated plant. Seed set is one flower, not total genet.
    y=np.column_stack((clean["seed #"].to_numpy(float),clean["unfert seeds"].to_numpy(float)))
    xx=sm.add_constant(pd.DataFrame({
       "ln_male":clean.ln_male.to_numpy(float),
       "st_dia":clean["st dia"].to_numpy(float),
       "stem_pos":clean["stem pos"].to_numpy(float),
       "day_may30":(clean["date"]=="140530").astype(float).to_numpy(float)}))
    m=sm.GLM(y,xx,family=sm.families.Binomial()).fit(cov_type="HC3")
    # Strong source missingness audit: all 130 assigned, measured fruit subset
    mod=smf.glm("has_seed_data ~ ln_male + Q('st dia')",data=x,
        family=sm.families.Binomial()).fit(cov_type="HC3")
    by={}
    for treatment,g in x.groupby("male_remaining"):
        d2=g.dropna(subset=["seed #","unfert seeds"]).copy()
        rr=d2["seed #"]/(d2["seed #"]+d2["unfert seeds"])
        by[str(int(treatment))]={"n_original_assigned_plants":len(g),
           "n_observed_fruit_seed":int(g["seed #"].notna().sum()),
           "n_complete_seed_ovule_pair":len(d2),
           "mean_mature_seed_number":scalar(g["seed #"].mean()),
           "mean_seed_weight_original_units":scalar(g["seed weight"].mean()),
           "mean_sampled_fruit_seed_fraction":scalar(rr.replace([np.inf,-np.inf],np.nan).mean())}
    return {"assigned_plant_count":130,"focal_fruit_seed_outcome_n":int(x["seed #"].notna().sum()),
       "after_complete_focal_outcome_n":len(clean),
       "male_flower_treatment_0_1_2_4_6":by,
       "plant_covariate_adjusted_male_display_to_focal_fruit_seedset_beta":{
         "beta":scalar(m.params["ln_male"]), "CI_95":scalar(m.conf_int().loc["ln_male"].to_numpy()),
         "p":scalar(m.pvalues["ln_male"]),
         "not_whole_plant_fitness":True},
       "outcome_missingness_related_to_assigned_male_count":fmtfit(mod),
       "sex_function":"female",
       "focal_seed_count_not_total_plant_seed_count":True,
       "source_seed_weight_units_uncertain":True}
def study4(d):
    x=d.copy()
    x["Total Visits"]=pd.to_numeric(x["Total Visits"],errors="coerce")
    x["patch size"]=pd.to_numeric(x["patch size"],errors="coerce")
    x=x.dropna(subset=["Total Visits","patch size","date"]).copy()
    x=x[x["patch size"]>0]
    x["ln_patch"]=np.log(x["patch size"])
    x["date"]=x.date.astype(str)
    fit=smf.glm("Q('Total Visits') ~ ln_patch + C(date)",data=x,
                 family=sm.families.NegativeBinomial(alpha=1)).fit(cov_type="HC3")
    return {"n_20minute_20flower_sessions":len(x),
       "n_zero_all_visits":int((x["Total Visits"]==0).sum()),
       "ln_patch_to_total_20min_visits_NB":fmtfit(fit),
       "patch_environment_not_randomized":True}
def study5(d):
    x=d.copy()
    req=["display treatment","patch size","grains","t to dep","t to wilt","st dia","date","location"]
    if not set(req)<=set(x.columns) or len(x)!=204:raise ValueError("unexpected original study5 source")
    for col in ("grains","t to dep","t to wilt","st dia"):
        x[col]=pd.to_numeric(x[col],errors="coerce")
    x=x.dropna(subset=req).copy()
    if len(x)!=204:raise ValueError("204 original study5 rows expected")
    x=x.rename(columns={"display treatment":"treatment","patch size":"patch","st dia":"st_dia","t to dep":"tdep","t to wilt":"twilt"})
    x["treatment"]=pd.Categorical(x.treatment.astype(str),categories=["1","6e","6"])
    x["patch"]=pd.Categorical(x.patch.astype(str),categories=["large","small"])
    if x.treatment.isna().any() or x.patch.isna().any():raise ValueError("unrecognized treatments")
    x["date"]=x.date.astype(str)
    x["location"]=x.location.astype(str).str.lower().str.strip()
    x["y"]=np.log1p(x.grains)
    x["fraction_lifespan_before_pollen"]=x.tdep/x.twilt
    x["fraction_lifespan_after_pollen"]=1-x["fraction_lifespan_before_pollen"]
    x["deposition_after_wilt"]=x.tdep>x.twilt
    allgroups={}
    for (tr,patch),g in x.groupby(["treatment","patch"],observed=True):
        allgroups[str(tr)+"_"+str(patch)]={
          "n":len(g),"mean_grains":scalar(g.grains.mean()),
          "mean_log1p_grains":scalar(g.y.mean()),
          "median_grains":scalar(g.grains.median()),
          "mean_days_until_first_pollen":scalar(g.tdep.mean()),
          "mean_days_until_wilt":scalar(g.twilt.mean()),
          "n_zero_stigma_pollen":int((g.grains==0).sum()),
          "mean_fraction_flower_life_before_first_pollen":scalar(g.fraction_lifespan_before_pollen.mean()),
          "n_deposition_after_wilt":int(g.deposition_after_wilt.sum())}
    if set(x.groupby(["treatment","patch"],observed=True).size())!={34}:
        raise ValueError("source factorial arms not 34 replicates each")
    formula="y ~ C(treatment, Treatment(reference='1')) * C(patch, Treatment(reference='large')) + st_dia + C(date)"
    model=smf.ols(formula,data=x).fit(cov_type="HC3")
    base6e="C(treatment, Treatment(reference='1'))[T.6e]"
    base6="C(treatment, Treatment(reference='1'))[T.6]"
    small="C(patch, Treatment(reference='large'))[T.small]"
    int6e=base6e+":"+small
    int6=base6+":"+small
    all_contrasts={
     "signal_only_emasculated6_vs_single_large_patch":model_contrast(model,base6e+" = 0"),
     "pollen_supplying6_vs_emasculated6_large_patch":model_contrast(model,base6+" - "+base6e+" = 0"),
     "signal_only_emasculated6_vs_single_small_patch":model_contrast(model,base6e+" + "+int6e+" = 0"),
     "pollen_supplying6_vs_emasculated6_small_patch":model_contrast(model,base6+" - "+base6e+" + "+int6+" - "+int6e+" = 0"),
     "site_small_vs_large_at_single_flower":model_contrast(model,small+" = 0"),
     "signal_only_patch_interaction":model_contrast(model,int6e+" = 0"),
     "pollen_supply_patch_interaction":model_contrast(model,int6+" - "+int6e+" = 0")
    }
    # Distinct time to first pollen measurement includes one source censored
    # zero stigma. Simple models exclude that one and do not purport Cox PH.
    timed=x[x.grains>0].copy()
    timed["ln_time_to_deposition"]=np.log(timed.tdep)
    timed["opportunity_life_frac"]=np.clip(1-timed.tdep/timed.twilt,0,1)
    timem=smf.ols("ln_time_to_deposition ~ C(patch) + C(treatment) + st_dia + C(date)",
                   data=timed).fit(cov_type="HC3")
    window=smf.ols("opportunity_life_frac ~ C(patch) + C(treatment) + st_dia + C(date)",
                    data=timed).fit(cov_type="HC3")
    return {"source_focal_stigmas":len(x),"factorial_replicates_per_arm":34,
       "six_arms":allgroups,
       "direct_three_arm_contrasts_on_log1p_pollen_grains":all_contrasts,
       "model_and_covariate_terms":fmtfit(model),
       "first_pollen_time_positive_event_n":len(timed),
       "first_pollen_zero_pollen_censored_n":int((x.grains==0).sum()),
       "ln_time_to_first_pollen_simple_OLS_not_survival_model":fmtfit(timem),
       "fraction_of_flower_life_remaining_after_first_pollen_simple_OLS":fmtfit(window),
       "mean_remaining_fraction_pollen_event_only":{
         str(k):scalar(g.opportunity_life_frac.mean()) for k,g in timed.groupby("patch",observed=True)},
       "pollen_grains_are_not_compatible_pollen_tubes":True,
       "n_site_locations":int(x.location.nunique()),
       "patch_treatment_not_fully_randomized":True,
       "stem_diameter_not_floral_construction_biomass":True,
       "censoring_not_modeled_in_OLS_time_analysis":True}

def main():
    version,frame=fetch_dataset()
    report={"status":"failed","original_source_doi":DOI,"published_version_id":version,
       "source_sha256":SHA,"original_source_counts":{str(k):len(z) for k,z in frame.items()},
       "denominator":"different independent experiments; do not multiply stage coefficients as observed fitness",
       "source_raw_individual_ids_not_exported":True}
    report["consecutive_probe_movement"]=study1(frame[1])
    report["self_vs_neighbor_pigment"]=study2(frame[2])
    report["focal_maternal_seed"]=study3(frame[3])
    report["patch_20minute_visitation"]=study4(frame[4])
    report["display_vs_pollen_presence_patch_factorial"]=study5(frame[5])
    report["status"]="success_original_five_study_stage_evidence_audit"
    (OUT/"clarkia_concinna_original_five_stage_gates_v13.json").write_text(
        json.dumps(scalar(report),indent=2,ensure_ascii=False)+"\n")
    print("CLARKIA_FIVE_STAGE_RESULT",json.dumps(scalar(report),ensure_ascii=False)[:37000],flush=True)
if __name__=="__main__":main()
