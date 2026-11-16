#!/usr/bin/env python3
"""Source-verified 2026 Silene dioica floral display -> pollen carryover -> paternity.

Joffard et al. 2026, AJB DOI 10.1002/ajb2.70276, Dataverse
DOI 10.57745/CF1VNW, published original public tabular files.
Two separate experiments: pollen-stigma transfer from sequential visits by
pollen-laden bumblebees, and maternity-fruit paternity in control/manipulated
flower displays. NEVER treat these as the same individual causal mediation chain.
No pollinator/plant/offspring/father IDs exported, no direct siring fitness.
"""
from __future__ import annotations
import csv,hashlib,io,json,math,urllib.request
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

OUT=Path("analysis_results_joffard2026_pollen_paternity_v09");OUT.mkdir(exist_ok=True)
ROOT="https://entrepot.recherche.data.gouv.fr/api/access/datafile/"
FILES={
 "paternity_assignment_results.tab":(765422,"0ae809e19d4b19f60e80ed6e81c2985ecdeb0882b7b857c449b567b325462a50"),
 "pollen_counts_on_stigmas.tab":(765423,"8291aea33a6c2894edccc55dd911f52a5c9c014049c20fb04ce7397f7e6dee68")
}
SEED=20261116;N_BOOT=5000;N_PERM=6000
GEN=np.random.default_rng(SEED)
def load():
    d={}
    for name,(fid,sha) in FILES.items():
        req=urllib.request.Request(ROOT+str(fid),headers={"User-Agent":"Peer-reviewed-pollen-transfer-audit-v09/0.1","Accept":"*/*"})
        with urllib.request.urlopen(req,timeout=80) as r:b=r.read(200000)
        if hashlib.sha256(b).hexdigest()!=sha:raise ValueError("immutable public source SHA mismatch "+name)
        txt=b.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
        frame=pd.read_csv(io.StringIO(txt),sep="\t",dtype="string")
        frame.columns=[str(x).strip() for x in frame.columns]
        d[name]=frame
    return d
def roundval(x):
    if isinstance(x,(np.integer,int)):return int(x)
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(np.floating,float)):return float(round(x,6)) if np.isfinite(x) else None
    if isinstance(x,(list,tuple)):return [roundval(v) for v in x]
    if isinstance(x,dict):return {str(k):roundval(v) for k,v in x.items()}
    return x

def effect(a,b,labels=("Control","Manipulated")):
    a=np.asarray(a,float);b=np.asarray(b,float)
    a=a[np.isfinite(a)];b=b[np.isfinite(b)]
    if len(a)<4 or len(b)<4:return {"status":"too_small","n_a":len(a),"n_b":len(b)}
    observed=b.mean()-a.mean()
    boot=np.mean(b[GEN.integers(0,len(b),(N_BOOT,len(b)))],axis=1)-np.mean(a[GEN.integers(0,len(a),(N_BOOT,len(a)))],axis=1)
    x=np.r_[a,b];n=len(a);pcount=0
    for _ in range(N_PERM):
        perm=GEN.permutation(len(x))
        d=x[perm[n:]].mean()-x[perm[:n]].mean()
        if abs(d)>=abs(observed)-1e-12:pcount+=1
    return {"groups":labels,"n_a":len(a),"n_b":len(b),"mean_a":roundval(a.mean()),"mean_b":roundval(b.mean()),
      "difference_b_minus_a":roundval(observed),"bootstrap_95":[roundval(y) for y in np.quantile(boot,[.025,.975])],
      "permutation_two_sided_p":roundval((pcount+1)/(N_PERM+1))}

def pollen(df):
    required=["pollinator_ID","position_of_the_flower_in_the_sequence","female_plant_ID","number_of_pollen_grains"]
    if not set(required)<=set(df):raise ValueError("unexpected pollen source schema")
    dd=df.copy()
    for c in ["position_of_the_flower_in_the_sequence","number_of_pollen_grains"]:
        dd[c]=pd.to_numeric(dd[c],errors="coerce")
    dd=dd.dropna(subset=required)
    if not dd.number_of_pollen_grains.ge(0).all() or not dd.position_of_the_flower_in_the_sequence.between(1,10).all():
        raise ValueError("pollen count/rank outside source domain")
    dd["bee"]=dd.pollinator_ID.astype(str)
    dd["pos"]=dd.position_of_the_flower_in_the_sequence.astype(int)
    dd["grains"]=dd.number_of_pollen_grains.astype(int)
    dd["lnpos"]=np.log(dd.pos)
    if len(dd)!=186:raise ValueError("unexpected original pollen transfer record count")
    n_bees=dd.bee.nunique()
    # Source only has 25 experimental pollinators; do not treat flowers as 186 independent bees.
    regressions={}
    for family,tag in ((sm.families.Poisson(),"poisson"),
                       (sm.families.NegativeBinomial(alpha=1),"negative_binomial_alpha1")):
        m=smf.glm("grains ~ lnpos",data=dd,family=family).fit(cov_type="cluster",cov_kwds={"groups":dd.bee})
        b=float(m.params["lnpos"]);lo,hi=[float(x) for x in m.conf_int().loc["lnpos"]]
        regressions[tag]={"n_independent_pollinators":n_bees,"n_visited_stigmas":len(dd),
          "coef_ln_position":roundval(b),"cluster95_ln_position":[roundval(lo),roundval(hi)],
          "p":roundval(m.pvalues["lnpos"]),
          "IRR_per_doubling_position":roundval(math.exp(b*math.log(2))),
          "cluster95_IRR_per_doubling":[roundval(math.exp(lo*math.log(2))),roundval(math.exp(hi*math.log(2)))],
          "overdispersion_pearson":roundval(m.pearson_chi2/m.df_resid)}
    # Additional fixed pollinator ID exposure control; cluster-robust by bee (n=25).
    m=smf.ols("np.log1p(grains) ~ lnpos + C(bee)",data=dd).fit(cov_type="cluster",cov_kwds={"groups":dd.bee})
    regressions["within_bee_fixed_effect_log1p"]={
      "coef_ln_position":roundval(m.params["lnpos"]),
      "cluster95":[roundval(x) for x in m.conf_int().loc["lnpos"]],
      "p":roundval(m.pvalues["lnpos"]),
      "n_pollinator_effects":n_bees}
    descriptive=[]
    for pos,g in dd.groupby("pos"):
        descriptive.append({"position":int(pos),"stigmas":len(g),"n_bees":g.bee.nunique(),
          "mean_pollen_grains":roundval(g.grains.mean()),
          "median_pollen_grains":roundval(g.grains.median()),
          "zero_fraction":roundval((g.grains==0).mean()),
          "total_grains":roundval(g.grains.sum())})
    # Early 1–3 and late 8–10 WITHIN the SAME pollinator. Bees without both
    # windows are excluded and their missing observations are never zero filled.
    early=dd[dd.pos<=3].groupby("bee").grains.mean()
    late=dd[dd.pos>=8].groupby("bee").grains.mean()
    common=early.index.intersection(late.index)
    if len(common)<7:raise ValueError("too few pollinators with early and late stigma samples")
    observed=(np.log1p(late.loc[common])-np.log1p(early.loc[common])).to_numpy(float)
    bootstrap=observed[GEN.integers(0,len(observed),(N_BOOT,len(observed)))].mean(axis=1)
    signs=GEN.choice([-1,1],size=(N_PERM,len(observed)))
    permmeans=np.mean(observed*signs,axis=1)
    two=(1+int((np.abs(permmeans)>=abs(observed.mean())).sum()))/(N_PERM+1)
    return {
      "experiment":"Controlled bumblebee sequence over ten female stigmas, not the natural female display manipulation",
      "n_stigmas":len(dd),"n_pollinators":n_bees,"n_female_plants_in_assay":dd.female_plant_ID.nunique(),
      "pollinator_position_sampling_counts":descriptive,"regressions":regressions,
      "early_vs_late_within_bee":{
         "n_bees_with_early_1to3_and_late_8to10":len(common),
         "mean_log1p_late_minus_early":roundval(observed.mean()),
         "bee_cluster_bootstrap95":[roundval(x) for x in np.quantile(bootstrap,[.025,.975])],
         "within_bee_signflip_p":roundval(two)},
      "positions_missing_not_imputed":True,
      "no_compatibility_or_fertilization_of_deposited_grains_tested":True}

def paternity(df):
    expected=["Offspring_ID","Treatment","Mother_ID","Fruit_ID","Candidate_father_ID"]
    if not set(expected)<=set(df):raise ValueError("unexpected original parentage file")
    if len(df)!=1150:raise ValueError("expected assigned 1150 offspring (not genotyped 1320)")
    x=df.copy()
    if x.Offspring_ID.duplicated().any():raise ValueError("duplicate assigned offspring")
    if not set(x.Treatment)=={"Control","Manipulated"}:raise ValueError("unexpected treatment labels")
    if x.groupby("Mother_ID").Treatment.nunique().max()!=1:raise ValueError("treatment conflict within mother")
    # Father distributions within fruits. 1150 is the assigned subset, NOT all genotyped seeds.
    fruit=[]
    for (mother,fruitid),g in x.groupby(["Mother_ID","Fruit_ID"]):
        n=len(g);cnt=g.Candidate_father_ID.value_counts()
        v=cnt.to_numpy(float)/n
        simpson=1-float(np.sum(v**2))
        eff=1/float(np.sum(v**2))
        fruit.append({"mother":str(mother),"fruit":str(fruitid),"treatment":str(g.Treatment.iloc[0]),
                      "n_offspring_assigned":n,"n_sires_observed":len(cnt),
                      "simpson_diversity":simpson,"simpson_effective_number_of_fathers":eff})
    fruit_df=pd.DataFrame(fruit)
    bymother=fruit_df.groupby(["mother","treatment"],as_index=False).agg(
       fruits=("fruit","size"),mean_sire_richness=("n_sires_observed","mean"),
       mean_effective_fathers=("simpson_effective_number_of_fathers","mean"),
       mean_simpson=("simpson_diversity","mean"),assigned_seeds=("n_offspring_assigned","sum"))
    assert len(bymother)==21 and fruit_df.shape[0]==55,(len(bymother),len(fruit_df))
    groups=bymother.treatment.value_counts().to_dict()
    if groups.get("Control")!=11 or groups.get("Manipulated")!=10:raise ValueError("unexpected original group sizes")
    results={"experiment":"Natural pollination female real flower-display increase; NOT controlled bee pollen carryover assay",
      "n_original_offspring_assigned":len(x),
      "n_original_1320_genotyped_seeds_not_all_assigned":1320,
      "n_unassigned_or_not_in_public_file":1320-len(x),
      "n_fruit_rows":len(fruit_df),"n_distinct_mothers":len(bymother),
      "mothers_by_treatment":{str(k):int(v) for k,v in groups.items()},
      "assigned_seeds_by_treatment":{str(k):int(v) for k,v in x.Treatment.value_counts().items()},
      "paternal_metrics":{}}
    for key in ["mean_sire_richness","mean_effective_fathers","mean_simpson","assigned_seeds","fruits"]:
        a=bymother.loc[bymother.treatment=="Control",key]
        b=bymother.loc[bymother.treatment=="Manipulated",key]
        results["paternal_metrics"][key]=effect(a,b)
    # For fruits from same mother, cosine similarity of father frequency vectors.
    # Avoid treating all pairwise fruit combinations as independent; summarize
    # to one statistic per mother, then permute/boot at *mother* level.
    founder={}
    import itertools
    for (mother,treatment),g in x.groupby(["Mother_ID","Treatment"]):
        sets={fruitid:c.Candidate_father_ID.value_counts(normalize=True)
              for fruitid,c in g.groupby("Fruit_ID")}
        sims=[]
        for a,b in itertools.combinations(sets,2):
            sx=sets[a];sy=sets[b];common=sx.index.union(sy.index)
            u=sx.reindex(common,fill_value=0).to_numpy(float)
            v=sy.reindex(common,fill_value=0).to_numpy(float)
            similarity=1-0.5*float(np.sum(np.abs(u-v)))
            sims.append(similarity)
        if sims:founder[str(mother)]={"treatment":str(treatment),"mean_bray_curtis_similarity":np.mean(sims),"fruit_pairs":len(sims)}
    mtable=pd.DataFrame(list(founder.values()))
    results["mean_among_fruit_paternal_profile_similarity"]=effect(
       mtable.loc[mtable.treatment=="Control","mean_bray_curtis_similarity"],
       mtable.loc[mtable.treatment=="Manipulated","mean_bray_curtis_similarity"])
    results["mothers_with_2plus_fruits"]=len(mtable)
    results["assigned_only_diversity_may_bias_inference"]=True
    results["statistical_unit_is_mother_not_each_seed_or_fruit"]=True
    results["paternal_diversity_is_not_total_male_offspring_siring"]=True
    return results

def main():
    raw=load()
    report={"status":"failed","study":"Joffard et al 2026 doi:10.1002/ajb2.70276",
      "source":"Official public Dataverse doi:10.57745/CF1VNW",
      "input_files_sha256_verified":{k:v[1] for k,v in FILES.items()},
      "permutation_seed":SEED,"bootstrap_reps":N_BOOT,"permutation_reps":N_PERM,
      "privacy":"No pollinator, mother, father, flower, fruit or offspring IDs in output"}
    report["pollen_carryover"]=pollen(raw["pollen_counts_on_stigmas.tab"])
    report["paternity"]=paternity(raw["paternity_assignment_results.tab"])
    report["status"]="success_original_dataverse_independent_experiments_reanalysis"
    (OUT/"joffard2026_pollen_sequence_and_paternity_v09.json").write_text(json.dumps(roundval(report),indent=2,ensure_ascii=False)+"\n")
    print("ORIGINAL_STUDY_REANALYSIS",json.dumps(roundval(report),ensure_ascii=False)[:32000],flush=True)
if __name__=="__main__":main()
