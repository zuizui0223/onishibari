#!/usr/bin/env python3
"""Johnson & Harder 2023 original 228-species pollen fates, secondarily reanalysed.
Source: Royal Society Proc B DOI 10.1098/rspb.2023.1148, Figshare
DOI 10.6084/m9.figshare.24146916.v1, author file 42364722.
RE=removed/produced; PTE=conspecific stigma receipt/removed;
RPE=RE*PTE. All are stage metrics, not genetically verified paternal fitness.
The log-stage shortfall decomposition is mathematical accounting, NOT a
new empirical causal evolutionary mechanism or optimal flower biomass.
"""
import hashlib,io,json,math,urllib.request
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.formula.api as smf
from scipy.stats import spearmanr
OUT=Path("analysis_results_global_pollen_export_v14");OUT.mkdir(exist_ok=True)
URL="https://ndownloader.figshare.com/files/42364722"
SHA="e5d981d0a7b3804a7d2309b2df6e1418b86a2d9e52d850175fe69ab5a4144a9c"
RNG=np.random.default_rng(20261116);BOOT=2200

def val(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return round(float(x),7) if np.isfinite(x) else None
    if isinstance(x,(list,tuple,np.ndarray)):return [val(y) for y in x]
    if isinstance(x,dict):return {str(k):val(v) for k,v in x.items()}
    return x
def fetch():
    req=urllib.request.Request(URL,headers={"User-Agent":"Published-228-species-pollen-stage-analysis/0.14","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=70) as r:b=r.read(500000)
    if hashlib.sha256(b).hexdigest()!=SHA:raise ValueError("Original Figshare species table SHA256 mismatch")
    excel=pd.ExcelFile(io.BytesIO(b))
    if not {"Species medians","Legend"}<=set(excel.sheet_names):raise ValueError("original sheets missing")
    d=pd.read_excel(excel,"Species medians")
    if len(d)!=228:raise ValueError("expected 228 species")
    need={"removpcnt","PTE","Poltype","groom","PollinatorClass","Nectar","Region","Family","familyGBIF"}
    if not need<=set(d):raise ValueError("original fields missing")
    for c in ["removpcnt","PTE"]:d[c]=pd.to_numeric(d[c],errors="coerce")
    if d.removpcnt.notna().sum()!=202 or d.PTE.notna().sum()!=224:
        raise ValueError("source RE/PTE counts changed")
    for c in ["removpcnt","PTE"]:
        if not d[c].dropna().between(.000001,100).all():raise ValueError(c+" outside percent range")
    d["family"]=d.familyGBIF.fillna(d.Family).astype(str)
    d["RE"]=d.removpcnt/100
    d["PTEfrac"]=d.PTE/100
    v=d.dropna(subset=["RE","PTEfrac"]).copy()
    v["RPE"]=v.RE*v.PTEfrac
    v["logloss_removal"]=-np.log(v.RE)
    v["logloss_transfer"]=-np.log(v.PTEfrac)
    v["logloss_total"]=-np.log(v.RPE)
    v["transfer_share"]=v.logloss_transfer/v.logloss_total
    # Actual grain mass-balance per PRODUCED pollen, distinct from logarithmic
    # efficiency shortfall. Stage losses must sum exactly to 1 with delivered.
    v["unremoved"]=1-v.RE
    v["removed_but_not_stigma"]=v.RE*(1-v.PTEfrac)
    v["delivered"]=v.RPE
    if not np.allclose(v.unremoved+v.removed_but_not_stigma+v.delivered,1.,atol=1e-12):
        raise ValueError("source pollen fate stage mass balance invalid")
    v["post_removal_lost_majority_of_all_losses"]=(
        v.removed_but_not_stigma>v.unremoved)

    if not np.allclose(v.logloss_removal+v.logloss_transfer,v.logloss_total,atol=1e-11):
        raise ValueError("stage-loss identity invalid")
    if len(v)<190:raise ValueError("too few complete species")
    return d,v

def rho(x,y):
    if len(x)<8 or np.unique(x).size<3 or np.unique(y).size<3:return np.nan
    return float(spearmanr(x,y).statistic)
def familyboot(v,stat):
    groups=[a for _,a in v.groupby("family",dropna=False)]
    if len(groups)<2:return {"n_families":len(groups),"CI95":None}
    results=[]
    for _ in range(BOOT):
        ids=RNG.integers(0,len(groups),size=len(groups))
        sample=pd.concat([groups[j] for j in ids],ignore_index=True)
        value=float(stat(sample))
        if np.isfinite(value):results.append(value)
    return {"n_families":len(groups),"n_valid":len(results),
       "CI95":val(np.quantile(results,[.025,.975])) if len(results)>.95*BOOT else None}
def main():
    d,v=fetch()
    report={"status":"success_original_228_species_pollen_stage_reanalysis",
      "original_paper":"10.1098/rspb.2023.1148",
      "original_figshare":"10.6084/m9.figshare.24146916.v1",
      "source_sha256":SHA,
      "n_species_original":len(d),"n_RE":int(d.RE.notna().sum()),
      "n_PTE":int(d.PTEfrac.notna().sum()),
      "n_both":len(v),"n_neither":int((d.RE.isna()&d.PTEfrac.isna()).sum()),
      "n_families_with_both":int(v.family.nunique()),
      "provenance":"species medians, not 228 direct field manipulations",
      "not_all_species_measured_both_stages":True,
      "pollen_removed_may_include_wind_and_rain":True,
      "conspecific_stigma_includes_possible_self_pollen":True,
      "cannot_estimate_fertilization_or_genetic_siring":True,
      "fate_summary":{},
      "by_pollen_packaging":{},"by_grooming":{},"packaging_by_grooming":{}}
    t=report["fate_summary"]
    t.update({
      "produced_pollen_mass_balance_species_equal_means":{
        "not_removed_pct":val(100*v.unremoved.mean()),
        "removed_not_delivered_pct":val(100*v.removed_but_not_stigma.mean()),
        "delivered_to_conspecific_stigma_pct":val(100*v.delivered.mean()),
        "sums_to_100":True},
      "n_post_removal_loss_dominates_in_produced_grain_mass_balance":
        int(v.post_removal_lost_majority_of_all_losses.sum()),
      "fraction_post_removal_mass_loss_dominates":
        val(v.post_removal_lost_majority_of_all_losses.mean()),
      "fraction_family_equal_post_removal_mass_loss_dominates":
        val(v.groupby("family").post_removal_lost_majority_of_all_losses.mean().mean()),
      "mean_RPE_equal_family_weighting":
        val(v.groupby("family").RPE.mean().mean()),
      "families_by_species_sample_size":{
        str(k):int(n) for k,n in v.family.value_counts().items()},
      "mass_balance_fraction_post_removal_family_boot":
        familyboot(v,lambda z:float(z.post_removal_lost_majority_of_all_losses.mean())),
      "mass_balance_percent_lost_after_removal_family_boot":
        familyboot(v,lambda z:float(100*z.removed_but_not_stigma.mean())),
      "median_RE":val(v.RE.median()),"median_PTE_given_removal":val(v.PTEfrac.median()),
      "median_RPE_produced_to_stigma":val(v.RPE.median()),
      "mean_RPE":val(v.RPE.mean()),"geomean_RPE":val(np.exp(np.log(v.RPE).mean())),
      "median_transport_share_of_log_shortfall":val(v.transfer_share.median()),
      "n_species_transport_loss_exceeds_removal_loss":int((v.logloss_transfer>v.logloss_removal).sum()),
      "fraction_species_transport_loss_exceeds_removal_loss":val((v.logloss_transfer>v.logloss_removal).mean()),
      "family_CI_median_RPE":familyboot(v,lambda d:d.RPE.median()),
      "family_CI_median_transport_logloss_share":familyboot(v,lambda d:d.transfer_share.median())})
    for kind,x in v.groupby("Poltype"):
        report["by_pollen_packaging"][str(kind)]={
            "n_species":len(x),"n_families":int(x.family.nunique()),
            "median_RE":val(x.RE.median()),"median_PTE":val(x.PTEfrac.median()),
            "median_RPE":val(x.RPE.median()),"median_transport_share":val(x.transfer_share.median()),
            "mass_balance_pct_unremoved_species_mean":val(100*x.unremoved.mean()),
            "mass_balance_pct_removed_not_delivered_species_mean":val(100*x.removed_but_not_stigma.mean()),
            "mass_balance_pct_stigma_delivered_species_mean":val(100*x.delivered.mean()),
            "mass_balance_fraction_with_larger_post_removal_loss":
                val(x.post_removal_lost_majority_of_all_losses.mean()),
            "n_transport_stage_dominant":int((x.logloss_transfer>x.logloss_removal).sum()),
            "family_boot_CI_median_RPE":familyboot(x,lambda z:z.RPE.median())}
    for kind,x in v.groupby("groom"):
        report["by_grooming"][str(kind)]={
            "n_species":len(x),"n_families":int(x.family.nunique()),
            "median_RE":val(x.RE.median()),"median_PTE":val(x.PTEfrac.median()),
            "median_RPE":val(x.RPE.median())}
    for (kind,groom),x in v.groupby(["Poltype","groom"]):
        report["packaging_by_grooming"][str(kind)+"_"+str(groom)]={
            "n":len(x),"n_families":int(x.family.nunique()),
            "median_RE":val(x.RE.median()),"median_PTE":val(x.PTEfrac.median()),
            "median_RPE":val(x.RPE.median())}
    report["nonidentity_RE_vs_PTE"]={
        "Spearman_rho":val(rho(v.RE,v.PTEfrac)),
        "family_boot_CI":familyboot(v,lambda z:rho(z.RE,z.PTEfrac)),
        "within_packaging":{str(k):{"n":len(z),"rho":val(rho(z.RE,z.PTEfrac))}
            for k,z in v.groupby("Poltype")},
        "not_causal_tradeoff":True,
        "RPE_correlation_with_PTE_is_algebraic_and_not_an_independent_test":True}
    c=v.dropna(subset=["groom","Poltype","Region"]).copy()
    c["logitRE"]=np.log(c.RE.clip(1e-5,1-1e-5)/(1-c.RE.clip(1e-5,1-1e-5)))
    c["logitPTE"]=np.log(c.PTEfrac.clip(1e-5,1-1e-5)/(1-c.PTEfrac.clip(1e-5,1-1e-5)))
    if len(c)>=130 and c.family.nunique()>=10:
        m=smf.ols("logitPTE~logitRE+C(Poltype)+C(groom)+C(Region)",data=c).fit(
          cov_type="cluster",cov_kwds={"groups":c.family})
        ci=m.conf_int().loc["logitRE"]
        report["RE_PTE_adjusted_noncausal"]={
           "n":len(c),"n_families":int(c.family.nunique()),
           "logitRE_beta":val(m.params["logitRE"]),
           "family_cluster_CI95":val(ci.to_numpy()),
           "p":val(m.pvalues["logitRE"]),
           "adjusted":["pollen_packaging","groom","region"],
           "family_clustering_does_not_equal_phylogenetic_control":True}
    (OUT/"johnson_harder_228_pollen_stages_v14.json").write_text(
       json.dumps(val(report),ensure_ascii=False,indent=2)+"\n")
    print("GLOBAL_228_STAGE_RESULTS",json.dumps(val(report),ensure_ascii=False)[:32000],flush=True)
if __name__=="__main__":main()
