#!/usr/bin/env python3
"""Two independently published experiments, original raw-data stage analysis.

1. Silene dioica: genuinely flowers vs fake visual signal, 2022 Barbot et al.
   43 each male/female × control/visual treatment, plus 20 hand-pollinated female
   controls. Individual trial assignments parsed from source ID PREFIX FC/FT/MC/MT/FS.
   Male fitness: published MEMM posterior median 'MedMRS', NOT 1721 unweighted
   father counts; female fitness: viable seed index 'prodTOT'.
2. Four Clarkia species: James & Geber 2025, 430 original FRUIT rows,
   175 candidate within-plant control + hand-supplement pairs, 80 controls-only.
   Experimental manipulation affects one flower, NOT whole plant.
   Plant pair proxy: (block,plot,focal_sp), only size-2 unique-key pairs used,
   with identical n_flower; verify and report other pairing/missingness.

Public Dryad published-version ZIPs only; pin original member SHA256.
No original plant/offspring/site IDs or raw observation rows exported.
Statistics are exploratory; same-genus Clarkia species not independent clades.
"""
from __future__ import annotations
from collections import Counter
import hashlib,io,json,math,urllib.request,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

OUT=Path("analysis_results_silene_clarkia_stage_v08");OUT.mkdir(exist_ok=True)
N_BOOT=3500
N_PERM=5000
SEED=20261116
RNG=np.random.default_rng(SEED)
DATA={
"silene":{"version":251371,"sha":{
 "dataID_FloDis.txt":"1d245bb2dd6222bdaa05b0e23521eafd286cfa533b3826f8a5d8a65f715d3f9a",
 "paternity_FloDis.txt":"2b440e3fad9cfd2801cba56e9de9941d13f4ee936cb9a79938d22610e62050b9",
 "allFruits_FloDis.txt":"9165c2c3c80d234b6395e91147d8f20a85fe210d0a76f28767e11b8790178bb2"}},
"clarkia":{"version":332601,"sha":{
 "clarkia_data_2017_kingsnake.csv":"a6b82a14761f1de282148195d9a2d56f6ba355812419253b8a2cb24014fa0568"}}
}
def get(key):
    spec=DATA[key]
    url=f"https://datadryad.org/api/v2/versions/{spec['version']}/download"
    req=urllib.request.Request(url,headers={"User-Agent":"Reproduction-original-data-stage-refit/0.1"})
    with urllib.request.urlopen(req,timeout=90) as r:b=r.read(750001)
    if len(b)>750000 or not zipfile.is_zipfile(io.BytesIO(b)):raise ValueError("source archive cap/format")
    out={}
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        for info in z.infolist():
            fn=Path(info.filename).name
            if fn not in spec["sha"]:continue
            raw=z.read(info)
            if hashlib.sha256(raw).hexdigest()!=spec["sha"][fn]:
                raise ValueError("uncompressed original source file SHA mismatch: "+fn)
            text=raw.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
            df=pd.read_csv(io.StringIO(text),sep=None,engine="python",keep_default_na=True)
            df.columns=[str(x).strip() for x in df.columns]
            out[fn]=df
    if set(out)!=set(spec["sha"]):raise ValueError("missing original source files "+str(set(spec["sha"])-set(out)))
    return out

def numeric(df,cs):
    for c in cs:df[c]=pd.to_numeric(df[c],errors="coerce")
    return df

def val(x):
    if isinstance(x,(np.integer,int)):return int(x)
    if isinstance(x,(np.floating,float)):return round(float(x),6) if np.isfinite(x) else None
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(list,tuple)):return [val(y) for y in x]
    if isinstance(x,dict):return {str(k):val(v) for k,v in x.items()}
    return x

def means_by_treatment(df,col,group):
    ans={}
    for n,z in df.groupby(group):
        x=pd.to_numeric(z[col],errors="coerce").dropna().to_numpy(float)
        ans[str(n)]={"n_nonmissing":len(x),"n_missing":len(z)-len(x),
                     "mean":val(x.mean()) if len(x) else None,
                     "median":val(np.median(x)) if len(x) else None,
                     "min":val(x.min()) if len(x) else None,"max":val(x.max()) if len(x) else None}
    return ans

def difference_test(a,b,log1p=False):
    """a=controls b=treatment, permutation respects individual randomization."""
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a[np.isfinite(a)];b=b[np.isfinite(b)]
    if log1p:
        if np.min(np.r_[a,b]) < 0:raise ValueError("negative outcome cannot log1p")
        a=np.log1p(a);b=np.log1p(b)
    if len(a)<8 or len(b)<8:return {"status":"insufficient","n_control":len(a),"n_treated":len(b)}
    obs=b.mean()-a.mean()
    r1=RNG.integers(0,len(a),(N_BOOT,len(a)));r2=RNG.integers(0,len(b),(N_BOOT,len(b)))
    draws=b[r2].mean(axis=1)-a[r1].mean(axis=1)
    pool=np.r_[a,b];na=len(a)
    count=0
    for _ in range(N_PERM):
        permutation=RNG.permutation(len(pool))
        d=pool[permutation[na:]].mean()-pool[permutation[:na]].mean()
        count+=abs(d)>=abs(obs)-1e-12
    return {"n_control":len(a),"n_treated":len(b),"treatment_minus_control":val(obs),
        "unpaired_bootstrap_95":[val(z) for z in np.quantile(draws,[.025,.975])],
        "two_sided_label_permutation_p":val((count+1)/(N_PERM+1)),
        "scale":"log1p_original_units" if log1p else "original_units",
        "experimental_groups_not_independent_species":True}

def silene():
    files=get("silene")
    df=files["dataID_FloDis.txt"].copy()
    df["group"]=df.ID.astype("string").str.extract(r"^([A-Za-z]+)",expand=False)
    counts=df.group.value_counts().to_dict()
    if {str(x):int(c) for x,c in counts.items()}!={"FC":43,"FT":43,"MC":43,"MT":43,"FS":20}:raise ValueError("unexpected study groups")
    cols=["prodTOT","mean_nbFlo","nbVisTot","nbVis","mean_nbOvTOT","MS_EchqueFem8",
          "MedMRS","MedFec","QttTot","TauxGerm","nb_fruit_OK","nb_fruit_MAX"]
    numeric(df,cols)
    results={"source":"Barbot et al 2022, doi:10.1098/rspb.2022.1987","source_dryad":"10.5061/dryad.0rxwdbs4d",
         "raw_n":len(df),"group_n":{str(k):int(v) for k,v in counts.items()},
         "treatment_mapping":"FC female control, FT female fake-added flowers, MC male control, MT male fake-added; FS female hand pollen (distinct from fake flower)",
         "source_metrics":"prodTOT viable maternal seed production; MedMRS male MEMM adjusted siring median; MS_EchqueFem8 mating success from 8 per mother resampling",
         "outcome_comparisons":{},"covariate_checks":{}}
    pairs={
      "female_visiting_insects":("FC","FT","nbVis",False),
      "female_total_visits":("FC","FT","nbVisTot",False),
      "male_visiting_insects":("MC","MT","nbVis",False),
      "male_total_visits":("MC","MT","nbVisTot",False),
      "female_viable_seeds":("FC","FT","prodTOT",True),
      "female_mating_partners":("FC","FT","MS_EchqueFem8",False),
      "male_siring_MEMM":("MC","MT","MedMRS",True),
      "male_mating_partners":("MC","MT","MS_EchqueFem8",False),
      "female_true_flower_count":("FC","FT","mean_nbFlo",False),
      "male_true_flower_count":("MC","MT","mean_nbFlo",False),
    }
    for key,(c,t,column,log) in pairs.items():
        a=df.loc[df.group==c,column].dropna()
        b=df.loc[df.group==t,column].dropna()
        results["outcome_comparisons"][key]={"test":difference_test(a,b,log1p=log),
          "mean_and_missing":means_by_treatment(df[df.group.isin([c,t])],column,"group"),
          "nominal_intervention":"dummy visual signal increased (not production of real flowers or gametes)"}
    for sex,c,t,outcome,gametes in (
        ("female","FC","FT","prodTOT","mean_nbOvTOT"),
        ("male","MC","MT","MedMRS","QttTot")):
        data=df[df.group.isin([c,t])].copy()
        data=data.dropna(subset=[outcome,"mean_nbFlo",gametes])
        data=data[(data[outcome]>=0)&(data.mean_nbFlo>=0)&(data[gametes]>=0)].copy()
        data["treated"]=data.group.eq(t).astype(int)
        data["logY"]=np.log1p(data[outcome])
        data["logF"]=np.log1p(data["mean_nbFlo"])
        data["logG"]=np.log1p(data[gametes])
        model=smf.ols("logY ~ treated + logF + logG",data=data).fit(cov_type="HC3")
        cis=model.conf_int()
        results["covariate_checks"][sex]={"n":len(data),"outcome":outcome,
         "treated":{"estimate":val(model.params["treated"]),
             "ci":[val(z) for z in cis.loc["treated"]],
             "p":val(model.pvalues["treated"])},
         "real_flower_number":{"estimate":val(model.params["logF"]),
              "ci":[val(z) for z in cis.loc["logF"]],
              "p":val(model.pvalues["logF"])},
         "gametes_per_flower":{"estimate":val(model.params["logG"]),
              "ci":[val(z) for z in cis.loc["logG"]],
              "p":val(model.pvalues["logG"])},
         "observational_real_flower_gamete_covariates_not_randomized":True}
    male=df[df.group.isin(["MC","MT"])].copy().sort_values("mean_nbFlo",ascending=False)
    # Authors eliminated one MT extreme; only do deterministic prior-specified comparison.
    mt=male[male.group=="MT"].sort_values("mean_nbFlo",ascending=False)
    if mt.shape[0]==43:
        results["male_max_display_treatment_sensitivity"]={
          "largest_mean_true_flower_count":val(mt.mean_nbFlo.iloc[0]),
          "second_largest_mean_true_flower_count":val(mt.mean_nbFlo.iloc[1]),
          "max_vs_second_difference":val(mt.mean_nbFlo.iloc[0]-mt.mean_nbFlo.iloc[1]),
          "after_excluding_single_largest_MT":{
           "male_MEMM_log1p":difference_test(
              male[male.group=="MC"].MedMRS,
              mt.iloc[1:].MedMRS,log1p=True)}
        }
    p=files["paternity_FloDis.txt"].copy()
    moms=set(df.loc[df.group.isin(["FC","FT","FS"]),"ID"].astype(str))
    fathers=set(df.loc[df.group.isin(["MC","MT"]),"ID"].astype(str))
    pref=p["Candidate_father_ID"].astype("string").str.extract(r"^([A-Za-z]+)",expand=False)
    results["raw_paternity_scope"]={
      "n_assigned_offspring_rows":len(p),
      "n_distinct_mothers":int(p["Mother_ID"].nunique()),
      "n_distinct_candidate_fathers":int(p["Candidate_father_ID"].nunique()),
      "matched_mother_rows":int(p.Mother_ID.astype(str).isin(moms).sum()),
      "matched_father_rows":int(p.Candidate_father_ID.astype(str).isin(fathers).sum()),
      "candidate_father_prefix_counts":{str(k):int(v) for k,v in pref.value_counts().items()},
      "note":"1721 assigned offspring != all 2030 offspring genotyped; do not treat raw assignment counts as unbiased male total genetic fitness; MEMM source MedMRS handles sampling."}
    return results

def bh_adjust(pvals):
    r=np.asarray(pvals,float);n=len(r)
    order=np.argsort(r);q=np.empty(n,float);prev=1.0
    for k in reversed(range(n)):
        idx=order[k];prev=min(prev,r[idx]*n/(k+1));q[idx]=prev
    return q

def clarkia():
    d=get("clarkia")["clarkia_data_2017_kingsnake.csv"].copy()
    required=["block","plot","focal_sp","bg_sp","treatment","n_flower",
            "n_ovule","n_ovule_assigned","n_seed","herbivory"]
    if not set(required)<=set(d.columns):raise ValueError("missing source col")
    for c in ("n_flower","n_ovule","n_ovule_assigned","n_seed","herbivory","bg_Cwt","n_forb","n_grass"):
        d[c]=pd.to_numeric(d[c],errors="coerce")
    if len(d)!=430:raise ValueError("unexpected Clarkia original fruit rows")
    if not set(d.treatment)=={"cont_seed","supp_seed"}:raise ValueError("unknown treatment")
    missing=d[["n_seed","n_ovule","n_flower"]].isna().sum().to_dict()
    if (d.n_seed>d.n_ovule).any():raise ValueError("seed exceeds ovule source row")
    grouped=d.groupby(["block","plot","focal_sp"],dropna=False)
    pairs=[]
    stats=Counter()
    for key,g in grouped:
        if len(g)==2 and set(g.treatment)=={"cont_seed","supp_seed"}:
            stats["source_two_treatment_plant_pairs"]+=1
            if g.n_flower.nunique(dropna=True)>1:stats["flower_count_mismatch"]+=1;continue
            if g[["n_seed","n_ovule","n_flower"]].isna().any().any():
                stats["incomplete_outcome_pairs"]+=1;continue
            if (g.n_ovule<=0).any():stats["nonpositive_ovule_pairs"]+=1;continue
            g=g.set_index("treatment")
            con=g.loc["cont_seed"];supp=g.loc["supp_seed"]
            pairs.append({
                "block":str(key[0]),"plot":str(key[1]),"species":str(key[2]),
                "n_flower":float(con.n_flower),
                "supp_nseed":float(supp.n_seed),"control_nseed":float(con.n_seed),
                "supp_seedset":float(supp.n_seed/supp.n_ovule),
                "control_seedset":float(con.n_seed/con.n_ovule),
                "supp_ovules":float(supp.n_ovule),"control_ovules":float(con.n_ovule),
                "any_assigned_ovule":bool((g.n_ovule_assigned==1).any()),
                "any_herbivory":bool((g.herbivory==1).any())})
        elif len(g)==1 and g.treatment.iloc[0]=="cont_seed":stats["controls_only"]+=1
        else:stats["unrecognized_group"]+=1
    p=pd.DataFrame(pairs)
    if len(p)<100:raise ValueError("too few complete original paired fruit samples: "+str(len(p)))
    p["delta_seedset"]=p.supp_seedset-p.control_seedset
    p["delta_seedcount"]=p.supp_nseed-p.control_nseed
    p["delta_ovulecount"]=p.supp_ovules-p.control_ovules
    p["delta_logseed"]=np.log1p(p.supp_nseed)-np.log1p(p.control_nseed)
    report={"source":"James and Geber Ecology Monographs 2025 doi:10.1002/ecm.1649",
      "source_dryad":"10.5061/dryad.k3j9kd56r",
      "original_fruit_rows":len(d),"source_missing_fields":{str(k):int(v) for k,v in missing.items()},
      "source_pair_filter":dict(stats),"n_complete_pairs":len(p),
      "n_independent_blocks":int(p.block.nunique()),
      "outcome_definition":"filled seeds divided by counted or source-assigned ovules, per individual fruit",
      "source_design":"one control and one pollen-supplemented flower within a focal plant; NOT whole-plant pollen intervention",
      "outcome_by_species":{}}
    pvals=[]
    spnames={"C":"Clarkia cylindrica","S":"Clarkia speciosa","U":"Clarkia unguiculata","X":"Clarkia xantiana"}
    for sp,g in p.groupby("species"):
        dseed=g.delta_seedset.to_numpy(float)
        n=len(g);n_blocks=g.block.nunique()
        # Full pair sign-flip test, one-sided not prespecified, use two-sided.
        perms=RNG.choice((-1.,1.),size=(N_PERM,n),replace=True)
        null=(perms*dseed).mean(axis=1)
        obs=dseed.mean()
        pval=(np.sum(np.abs(null)>=abs(obs)-1e-12)+1)/(N_PERM+1)
        # Cluster block bootstrap; blocks resampled, not individual fruit pairs.
        blocks=list(g.block.unique())
        byblock=[g[g.block==b]["delta_seedset"].to_numpy(float) for b in blocks]
        boot=np.empty(N_BOOT)
        for i in range(N_BOOT):
            choice=RNG.integers(0,len(blocks),size=len(blocks))
            selected=np.concatenate([byblock[k] for k in choice])
            boot[i]=selected.mean()
        info={
            "name":spnames.get(sp,sp),"n_complete_pairs":n,
            "n_distinct_blocks":n_blocks,
            "mean_control_seedset":val(g.control_seedset.mean()),
            "mean_supp_seedset":val(g.supp_seedset.mean()),
            "mean_seedset_difference":val(obs),
            "block_bootstrap95_seedset_difference":[val(x) for x in np.quantile(boot,[.025,.975])],
            "within_plant_signflip_p":val(pval),
            "mean_control_fruit_seed_count":val(g.control_nseed.mean()),
            "mean_supp_fruit_seed_count":val(g.supp_nseed.mean()),
            "mean_pair_seed_count_difference":val(g.delta_seedcount.mean()),
            "mean_pair_ovule_count_difference":val(g.delta_ovulecount.mean()),
            "n_assigned_ovule_pairs":int(g.any_assigned_ovule.sum()),
            "n_herbivory_pairs":int(g.any_herbivory.sum()),
            "median_real_flowers_on_plant":val(g.n_flower.median()),
            "not_measured_whole_plant_total_seed_production":True,
            "sensitivity_exclude_assigned_ovule":{
              "n_pairs":int((~g.any_assigned_ovule).sum()),
              "mean_difference":val(g.loc[~g.any_assigned_ovule,"delta_seedset"].mean())},
            "sensitivity_exclude_herbivory":{
              "n_pairs":int((~g.any_herbivory).sum()),
              "mean_difference":val(g.loc[~g.any_herbivory,"delta_seedset"].mean())}
        }
        report["outcome_by_species"][sp]=info
        pvals.append((sp,pval))
    q=bh_adjust([v for _,v in pvals])
    for i,(sp,_) in enumerate(pvals):
        report["outcome_by_species"][sp]["exploratory_BH_q_over_four_species"]=val(q[i])
    report["boundary_note"]="Fruit-level pollen supplementation estimates cannot be scaled to plant lifetime offspring without total seed counts on unsupplemented and entirely supplemented plants."
    report["species_not_independent_phylogenetic_replicates"]=True
    return report

def main():
    r={"status":"failed","methods":"original public version ZIP SHA-verified; individual and block IDs are never exported",
        "no_plant_id_or_offspring_id_or_plot_id_exported":True,
        "limits":["No meta regression pooling Silene and Clarkia outcomes: different intervention and fitness denominators",
                  "2026 separate Silene paternal diversity paper remains independently published research, not these 2022 observations"],
        "seed":SEED,"bootstrap_reps":N_BOOT,"permutation_reps":N_PERM}
    r["silene_2022"]=silene()
    r["clarkia_2025"]=clarkia()
    r["status"]="success_original_raw_two_experiment_reanalysis"
    (OUT/"silene_clarkia_stage_separation_v08.json").write_text(json.dumps(val(r),ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(val(r),ensure_ascii=False,indent=2),flush=True)
if __name__=="__main__":main()
