#!/usr/bin/env python3
"""Original Joffard 2026 Silene dioica pollen carryover + paternity refit.

Two INDEPENDENT experimental datasets with different units:
A) 186 collected stigmas after sequenced bumblebee probes; 25 independent
   pollinator IDs, repeated positions 1..10. Grain counts are on stigmas,
   NOT confirmed compatible pollen tubes, seeds or future recruits.
B) 1150 confidently assigned offspring in 55 fruits nested in 21 maternal
   plants (10 visual-display-manipulated, 11 controls). Treatments assigned
   to MOTHERS. Paternity profile/diversity is NOT total paternal offspring
   success or fecundity.

Official anonymous Dataverse file API. Exact SHA256 of two public source files
pinned from independent v0.9 inventory. Report aggregates only, not IDs.
No inaccessible genotyped seeds are imputed as failures.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import hashlib,io,json,math,urllib.request
from itertools import combinations
from pathlib import Path
import numpy as np,pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

OUT=Path("analysis_results_joffard2026_original_refit_v09");OUT.mkdir(exist_ok=True)
ROOT="https://entrepot.recherche.data.gouv.fr/api/access/datafile/"
DATA={
  "paternity":(765422,"0ae809e19d4b19f60e80ed6e81c2985ecdeb0882b7b857c449b567b325462a50"),
  "pollen":(765423,"8291aea33a6c2894edccc55dd911f52a5c9c014049c20fb04ce7397f7e6dee68")
}
SEED=20261116
RNG=np.random.default_rng(SEED)
BOOT=3500
PERM=6000

def load(name):
    fileid,sha=DATA[name]
    req=urllib.request.Request(ROOT+str(fileid),headers={"User-Agent":"Original-Silene-2026-pollen-paternity-stage-audit/0.1","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=60) as r:b=r.read(150000)
    if hashlib.sha256(b).hexdigest()!=sha:raise ValueError("source archive checksum mismatch: "+name)
    x=pd.read_csv(io.BytesIO(b),sep="\t",dtype="string")
    x.columns=[c.strip() for c in x.columns]
    return x

def rnd(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return round(float(x),6) if np.isfinite(x) else None
    if isinstance(x,dict):return {str(k):rnd(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [rnd(z) for z in x]
    return x

def pearson_FE_slope(xy,group):
    z=xy.copy()
    z["xc"]=z["pos"]-z.groupby(group)["pos"].transform("mean")
    z["yc"]=z["log_grains"]-z.groupby(group)["log_grains"].transform("mean")
    return float((z.xc*z.yc).sum()/(z.xc*z.xc).sum())

def pollen():
    x=load("pollen")
    required=["pollinator_ID","position_of_the_flower_in_the_sequence","female_plant_ID","number_of_pollen_grains"]
    if list(x.columns)!=required:raise ValueError("unexpected pollen schema "+str(list(x.columns)))
    # Patsy/statsmodels categorical formula terms require ordinary Python
    # object strings, not pandas StringDtype, with some runner versions.
    x["pollinator_ID"]=x["pollinator_ID"].astype(str)
    x["female_plant_ID"]=x["female_plant_ID"].astype(str)
    x["pos"]=pd.to_numeric(x[required[1]],errors="coerce")
    x["grains"]=pd.to_numeric(x["number_of_pollen_grains"],errors="coerce")
    if x[required].isna().any().any() or x.pos.isna().any() or x.grains.isna().any():raise ValueError("missing pollen source")
    if (x.grains<0).any() or ((x.pos<1)|(x.pos>10)).any():raise ValueError("out-of-range pollen data")
    x["log_grains"]=np.log1p(x.grains.astype(float))
    if len(x)!=186 or x.pollinator_ID.nunique()!=25:raise ValueError("source rows/group mismatch")
    count_per_poll=x.groupby("pollinator_ID").size()
    repeated_pollinator_position=x.duplicated(["pollinator_ID","pos"]).sum()
    if repeated_pollinator_position>0:print("WARN_duplicate_position_per_pollinator",repeated_pollinator_position,flush=True)
    perpos={}
    for p,g in x.groupby("pos"):
        gg=g.grains.to_numpy(float)
        perpos[str(int(p))]={"n_stigmas":len(g),"mean_grains":gg.mean(),"median_grains":np.median(gg),
                             "n_zero":int(np.sum(gg==0)),"fraction_zero":np.mean(gg==0),
                             "p90_grains":np.quantile(gg,.9)}
    fit=smf.ols("log_grains ~ pos + C(pollinator_ID)",data=x).fit(cov_type="cluster",cov_kwds={"groups":x.pollinator_ID})
    b=float(fit.params["pos"]);ci=fit.conf_int().loc["pos"]
    model_poisson=smf.glm("grains ~ pos + C(pollinator_ID)",data=x,
        family=sm.families.Poisson()).fit(cov_type="cluster",cov_kwds={"groups":x.pollinator_ID})
    bp=float(model_poisson.params["pos"]);cip=model_poisson.conf_int().loc["pos"]
    # Within-pollinator fixed-effect re-estimation with independent bee-level
    # resampling: avoid calling 186 stigmas independent biological replicates.
    xy=x[["pollinator_ID","pos","log_grains"]].copy()
    xy["pos"]=xy["pos"].astype(float)
    groups=[]
    for _,g in xy.groupby("pollinator_ID"):
        xx=g.pos.to_numpy(float); yy=g.log_grains.to_numpy(float)
        xx=xx-xx.mean(); yy=yy-yy.mean()
        groups.append((xx,yy))
    terms=np.array([(float(np.dot(x,y)),float(np.dot(x,x))) for x,y in groups])
    idx=RNG.integers(0,len(groups),size=(BOOT,len(groups)))
    boots=terms[idx].sum(axis=1)
    bci=np.quantile(boots[:,0]/boots[:,1],[.025,.975])
    # shuffle sequence order within bee, and not across bees.
    obs=sum(np.dot(aa,bb) for aa,bb in groups)/sum(np.dot(aa,aa) for aa,bb in groups)
    null=[]
    for i in range(PERM):
        numerator=sum(float(np.dot(aa,RNG.permutation(bb))) for aa,bb in groups)
        denom=sum(float(np.dot(aa,aa)) for aa,bb in groups)
        null.append(numerator/denom)
    pperm=(sum(abs(v)>=abs(obs)-1e-12 for v in null)+1)/(PERM+1)
    # Does the result depend on one extreme stigma?
    q99=x.grains.quantile(.99)
    trim=x[x.grains<=q99].copy()
    ft=smf.ols("log_grains ~ pos + C(pollinator_ID)",data=trim).fit(
         cov_type="cluster",cov_kwds={"groups":trim.pollinator_ID})
    # Leave-one-bee-out and exclude-first-visit sensitivity; no claim
    # that the 1st and 10th visit involve identical bee cohorts.
    loo=[]
    for bee in x.pollinator_ID.unique():
        subset=xy[xy.pollinator_ID!=bee].copy()
        loo.append(pearson_FE_slope(subset,"pollinator_ID"))
    withoutfirst=x[x.pos>=2].copy()
    fit_withoutfirst=smf.ols("log_grains ~ pos + C(pollinator_ID)",
        data=withoutfirst).fit(cov_type="cluster",
                               cov_kwds={"groups":withoutfirst.pollinator_ID})
    # Original archive contains one double same-bee / same-sequence-position
    # record; average that one cell rather than inflate its influence.
    collapsed=x.groupby(["pollinator_ID","pos"],as_index=False).agg(
        log_grains=("log_grains","mean"),grains=("grains","mean"))
    fit_collapsed=smf.ols("log_grains ~ pos + C(pollinator_ID)",
       data=collapsed).fit(cov_type="cluster",
           cov_kwds={"groups":collapsed.pollinator_ID})
    return {
       "source":"Joffard et al 2026 AJB doi:10.1002/ajb2.70276; data doi:10.57745/CF1VNW",
       "n_stigmas":len(x),"n_independent_bumblebees":x.pollinator_ID.nunique(),
       "n_distinct_female_recipient_plants":x.female_plant_ID.nunique(),
       "positions":perpos,"n_zero_grain_stigmas":int((x.grains==0).sum()),
       "grain_sum":float(x.grains.sum()),
       "bee_num_stigmas_per_sequence":{"min":int(count_per_poll.min()),"median":float(count_per_poll.median()),"max":int(count_per_poll.max())},
       "duplicate_pollinator_position_rows":int(repeated_pollinator_position),
       "bee_cluster_influence_sensitivity":{
         "leave_one_bee_out_beta_range":[float(min(loo)),float(max(loo))],
         "drop_first_visit_n":len(withoutfirst),
         "drop_first_visit_beta":float(fit_withoutfirst.params["pos"]),
         "drop_first_visit_95":[float(v) for v in fit_withoutfirst.conf_int().loc["pos"]],
         "collapse_duplicate_bee_position_n":len(collapsed),
         "collapse_duplicate_beta":float(fit_collapsed.params["pos"]),
         "collapse_duplicate_95":[float(v) for v in fit_collapsed.conf_int().loc["pos"]]},

       "bee_FE_log1p_grain_per_additional_position":{
         "beta":b,"cluster95":[float(ci.iloc[0]),float(ci.iloc[1])],"p":float(fit.pvalues["pos"]),
         "cluster_bootstrap95":bci,"within_bee_shuffle_p":pperm,
         "multiplicative_ratio_of_1plus_grains_per_next_position":math.exp(b)},
       "bee_FE_Poisson_grain_count_sensitivity":{
         "log_IRR":bp,"IRR":math.exp(bp),"cluster95_IRR":[math.exp(float(cip.iloc[0])),math.exp(float(cip.iloc[1]))],
         "p":float(model_poisson.pvalues["pos"]),"pearson_overdispersion":float(model_poisson.pearson_chi2/model_poisson.df_resid)},
       "top_1pct_grains_trimmed_log1p_FE":{
         "n_retained":len(trim),"p99_threshold":float(q99),
         "beta":float(ft.params["pos"]),
         "95":[float(v) for v in ft.conf_int().loc["pos"]],
         "p":float(ft.pvalues["pos"])},
       "original_dataverse_row_ids_not_exported":True,
       "stigma_conspecific_deposition_not_equivalent_to_compatible_tubes_or_seeds":True,
       "visits_are_sequenced_within_bee_not_display_manipulation":True
    }

def permutation_bootstrap_difference(frame,col,treatment_col="treatment",nperm=PERM):
    arr=[frame[frame[treatment_col]==k][col].dropna().to_numpy(float) for k in ("Control","Manipulated")]
    if not all(len(a)>=7 for a in arr):return {"status":"too_small_for_inference"}
    a,b=arr
    delta=float(b.mean()-a.mean())
    b0=RNG.choice(a,size=(BOOT,len(a)),replace=True).mean(axis=1)
    b1=RNG.choice(b,size=(BOOT,len(b)),replace=True).mean(axis=1)
    ci=np.quantile(b1-b0,[.025,.975])
    combined=np.r_[a,b];n0=len(a);over=0
    for _ in range(nperm):
        z=RNG.permutation(combined)
        d=z[n0:].mean()-z[:n0].mean()
        over+=abs(d)>=abs(delta)-1e-12
    return {"n_control":len(a),"n_manipulated":len(b),
      "control_mean":float(a.mean()),"manipulated_mean":float(b.mean()),
      "manipulated_minus_control":delta,
      "bootstrap95":ci,
      "randomization_permutation_two_sided_p":(over+1)/(nperm+1)}

def expected_sires_n(father_freq,m=12):
    n=int(sum(father_freq))
    if n<m:return None
    denom=math.comb(n,m)
    return float(sum(1.0 - (math.comb(n-int(count),m)/denom if n-int(count)>=m else 0.0) for count in father_freq))

def paternity():
    x=load("paternity")
    required=["Offspring_ID","Treatment","Mother_ID","Fruit_ID","Candidate_father_ID"]
    if list(x.columns)!=required:raise ValueError("unexpected paternity columns")
    if x[required].isna().any().any():raise ValueError("missing assigned offspring IDs")
    if len(x)!=1150 or x.Mother_ID.nunique()!=21:raise ValueError("unexpected paternity grain")
    if not set(x.Treatment.unique())=={"Control","Manipulated"}:raise ValueError("unknown treatment")
    if x.Offspring_ID.duplicated().any():raise ValueError("duplicate offspring in archive")
    if x.groupby("Mother_ID").Treatment.nunique().max()!=1:raise ValueError("mother assigned different treatments")
    if x.groupby("Fruit_ID").Mother_ID.nunique().max()!=1:raise ValueError("Fruit_ID reused across mothers unexpectedly")
    fruits=[]
    for (mother,fruit,treat),g in x.groupby(["Mother_ID","Fruit_ID","Treatment"]):
        n=len(g)
        freqs=g.Candidate_father_ID.value_counts()
        p=(freqs/n).to_numpy(float)
        simpson=float((1.0-(p*p).sum())*n/(n-1)) if n>1 else None
        fruits.append({"mother":mother,"fruit":fruit,"treatment":treat,"assigned_seeds":n,
                       "sire_count":len(freqs),"unbiased_simpson":simpson,
                       "effective_sire_number":float(1/(p*p).sum()),
                       "rarefied_expected_sires12":expected_sires_n(freqs.to_numpy(),12),
                       "freqs":{str(k):int(v) for k,v in freqs.items()}})
    mothers=[]
    for (mother,treat),g in pd.DataFrame(fruits).groupby(["mother","treatment"]):
        profiles=list(g.freqs)
        within=[]
        for i,j in combinations(range(len(profiles)),2):
            a=profiles[i];b=profiles[j]
            na=sum(a.values()); nb=sum(b.values())
            overlap=sum(min(a.get(s,0)/na,b.get(s,0)/nb) for s in set(a)|set(b))
            within.append(float(1-overlap))
        mothers.append({"mother":mother,"treatment":treat,"n_fruits":len(g),
          "n_assigned_seeds":int(g.assigned_seeds.sum()),
          "mean_assigned_seeds_per_fruit":float(g.assigned_seeds.mean()),
          "mean_sire_count":float(g.sire_count.mean()),
          "mean_unbiased_simpson":float(g.unbiased_simpson.mean()),
          "mean_effective_sire_number":float(g.effective_sire_number.mean()),
          "mean_rarefied_sires12":float(g.rarefied_expected_sires12.mean()),
          "n_fruits_rarefaction_eligible":int(g.rarefied_expected_sires12.notna().sum()),
          "within_mother_mean_fruit_sire_dissimilarity":np.mean(within) if within else None,
          "n_fruit_pairs":len(within)})
    m=pd.DataFrame(mothers)
    if m.mother.nunique()!=21:raise ValueError("21 mother mismatch")
    counts=m.treatment.value_counts().to_dict()
    if counts.get("Control")!=11 or counts.get("Manipulated")!=10:raise ValueError("mother group totals unexpected "+str(counts))
    metrics=("mean_sire_count","mean_unbiased_simpson","mean_effective_sire_number",
             "mean_rarefied_sires12","within_mother_mean_fruit_sire_dissimilarity")
    tests={}
    for metric in metrics:
        tests[metric]=permutation_bootstrap_difference(m,metric)
    return {"source":"Joffard et al 2026 AJB doi:10.1002/ajb2.70276; archived data doi:10.57745/CF1VNW",
       "n_assigned_offspring":len(x),"n_independent_mothers":len(m),"n_fruits":len(fruits),
       "n_assigned_candidate_fathers":x.Candidate_father_ID.nunique(),
       "mother_treatments":{str(k):int(v) for k,v in counts.items()},
       "assigned_seeds_by_treatment":{str(k):int(v) for k,v in x.Treatment.value_counts().items()},
       "fruit_sampling_n":{"min":min(f["assigned_seeds"] for f in fruits),
                           "median":float(np.median([f["assigned_seeds"] for f in fruits])),
                           "max":max(f["assigned_seeds"] for f in fruits)},
       "n_fruits_with_12_or_more_assigned":sum(f["assigned_seeds"]>=12 for f in fruits),
       "n_mothers_with_2plus_fruits":int((m.n_fruits>=2).sum()),
       "treatment_effects_mother_as_unit":tests,
       "diversity_scope":"Assigned sire composition and number, not absolute lifetime paternal fitness",
       "original_paper_1320_genotyped_not_all_assigned":True,
       "unassigned_or_unobserved_seeds_not_in_denominator":True,
       "rarefaction12_exact_hypergeometric_expectation_per_fruit":True,
       "original_dataverse_ids_not_exported":True}

def main():
    result={"status":"failed","reproducibility":"official public Dataverse file API; pinned uncompressed file SHA256",
       "data_DOI":"10.57745/CF1VNW","paper_DOI":"10.1002/ajb2.70276","seed":SEED,
       "boot_n":BOOT,"permutation_n":PERM,"source_ids_not_exported":True}
    result["pollen_visitation_sequence"]=pollen()
    result["paternal_diversity_display_experiment"]=paternity()
    result["status"]="success_original_pollen_paternity_two_experiment_reanalysis"
    (OUT/"joffard2026_stage_pollen_vs_genetic_compatibility_v09.json").write_text(json.dumps(rnd(result),indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(rnd(result),indent=2,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
