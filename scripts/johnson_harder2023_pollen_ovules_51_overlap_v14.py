#!/usr/bin/env python3
"""Cross-paper overlap and pollen/ovule differential response to PTE (v0.14).
Data are ORIGINAL authors' 2023 Dryad 51 species DOI 10.5061/dryad.w6m905qtq,
not an independent source replication of the same researchers' Figshare 228.
If genus+species match, avoid treating the two published rows as independent.
Comparisons are noncausal and do not directly observe pollen export fathering.
"""
import hashlib,io,json,urllib.request,zipfile
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.formula.api as smf
from johnson_harder2023_pollen_fates_v14 import fetch,val,RNG
OUT=Path("analysis_results_global_pollen_po_link_v14");OUT.mkdir(exist_ok=True)
URL="https://datadryad.org/api/v2/versions/250874/download"
ZSHA="c37a24d4d0d791e5657402be51e620c3a27992092677884a6e22ed2466c794ad"
CSHA="46296a24715a49de964df87d327977b25f6058480cb7e82c9126f2d7facadb05"
BOOT=2800
def namekey(g,s):
    if pd.isna(g) or pd.isna(s):return None
    g=str(g).strip().lower().split()[0]
    x=str(s).strip().lower().split()
    if not g or not x:return None
    epi=x[1] if x[0]==g and len(x)>1 else x[0]
    if not epi.isalpha() or epi in ("sp","ssp","cf","aff"):return None
    return g+" "+epi
def dryad():
    req=urllib.request.Request(URL,headers={"Accept":"*/*","User-Agent":"Original-PTE-PO-cross-comparison-v14"})
    with urllib.request.urlopen(req,timeout=70) as r:b=r.read(800000)
    with zipfile.ZipFile(io.BytesIO(b)) as z:inner=z.read("PTE_PO_data.zip")
    if hashlib.sha256(inner).hexdigest()!=ZSHA:raise ValueError("published inner source changed")
    with zipfile.ZipFile(io.BytesIO(inner)) as z:
        matches=[i for i in z.infolist() if Path(i.filename).name=="PTE_PO_data_matched.csv"]
        if len(matches)!=1:raise ValueError("expected one published original species CSV in nested archive")
        csvbytes=z.read(matches[0])
    if hashlib.sha256(csvbytes).hexdigest()!=CSHA:raise ValueError("original 51-species source changed")
    d=pd.read_csv(io.BytesIO(csvbytes))
    if len(d)!=51 or not {"genusGBIF","speciesGBIF","familyGBIF","PTE","Pollen","Ovules","PO"}<=set(d):
        raise ValueError("source field mismatch")
    for k in ("PTE","Pollen","Ovules","PO"):d[k]=pd.to_numeric(d[k],errors="coerce")
    if d[["PTE","Pollen","Ovules"]].isna().any().any() or (d[["PTE","Pollen","Ovules"]]<=0).any().any():
        raise ValueError("original species missing source fates or counts")
    d["family"]=d.familyGBIF.astype(str)
    d["namekey"]=[namekey(g,s) for g,s in zip(d.genusGBIF,d.speciesGBIF)]
    return d
def familyboot(d,col):
    fam=[q for _,q in d.groupby("family")]
    if len(fam)<2:return None
    values=[]
    for i in range(BOOT):
        idx=RNG.integers(0,len(fam),len(fam))
        x=pd.concat([fam[j] for j in idx],ignore_index=True)
        if x.lnPTE.std()<1e-10:continue
        try:r=smf.ols(f"{col} ~ lnPTE",data=x).fit().params["lnPTE"]
        except Exception:continue
        if np.isfinite(r):values.append(float(r))
    return {"n_family_groups":len(fam),"n_valid":len(values),
            "CI95":val(np.quantile(values,[.025,.975])) if len(values)>BOOT*.95 else None}
def main():
    a= dryad()
    b,_=fetch()
    print("TAXON_SOURCE_FORMAT_DIAGNOSTIC",json.dumps({
       "source228":b[["genusGBIF","speciesGBIF"]].head(5).astype(str).to_dict(orient="records"),
       "source51":a[["genusGBIF","speciesGBIF"]].head(5).astype(str).to_dict(orient="records"),
       "source228_dtypes":{k:str(v) for k,v in b[["genusGBIF","speciesGBIF"]].dtypes.items()},
       "source51_dtypes":{k:str(v) for k,v in a[["genusGBIF","speciesGBIF"]].dtypes.items()}},ensure_ascii=False),flush=True)
    a["lnPollen"]=np.log(a.Pollen);a["lnOvules"]=np.log(a.Ovules)
    a["lnPTE"]=np.log(a.PTE);a["lnPO"]=a.lnPollen-a.lnOvules
    fam=a.family.nunique()
    results={"status":"success_original_51_species_PTE_vs_pollen_ovule_separate_source",
       "author_source":"10.5061/dryad.w6m905qtq","source_SHA256":CSHA,
       "main_author_paper":"10.1002/ajb2.16177",
       "n_species":len(a),"n_families":int(fam),
       "median_PTE_percent":val(a.PTE.median()),
       "median_pollen_per_flower":val(a.Pollen.median()),
       "median_ovules_per_flower":val(a.Ovules.median()),
       "do_not_treat_same_author_papers_as_independent_evidence":True,
       "regressions":{},"cross_paper_taxon_overlap":{}}
    for col in ("lnPollen","lnOvules","lnPO"):
        fit=smf.ols(f"{col} ~ lnPTE",data=a).fit(
           cov_type="cluster" if fam>=10 else "HC3",
           cov_kwds={"groups":a.family} if fam>=10 else None)
        ci=fit.conf_int().loc["lnPTE"]
        los=[]
        for f in a.family.unique():
            z=a[a.family!=f]
            if len(z)>=20:los.append(float(smf.ols(f"{col}~lnPTE",data=z).fit().params["lnPTE"]))
        results["regressions"][col]={
           "beta":val(fit.params["lnPTE"]),"CI95":val(ci.to_numpy()),
           "p":val(fit.pvalues["lnPTE"]),
           "leave_one_family_out_beta_range":val([min(los),max(los)]) if los else None,
           "family_bootstrap":familyboot(a,col),
           "family_cluster_SE_only_if_n_families_at_least10":fam>=10}
    results["pollen_minus_ovule_slope_difference"]={
       "beta":results["regressions"]["lnPO"]["beta"],
       "note":"exactly ln(Pollen/Ovules) slope; not an independent new metric or lifetime fitness"}
    b["namekey"]=[namekey(g,s) for g,s in zip(b.genusGBIF,b.speciesGBIF)]
    orig=b.dropna(subset=["namekey"]).copy()
    origu=orig[orig.namekey.map(orig.namekey.value_counts())==1]
    qa=a.dropna(subset=["namekey"]).copy()
    qa=qa[qa.namekey.map(qa.namekey.value_counts())==1]
    join=origu[["namekey","PTE","Poltype"]].merge(
      qa[["namekey","PTE","Pollen","Ovules"]],on="namekey",suffixes=("_228","_51"))
    overlap={"n_exact_binomial_matches":len(join),
       "n_228_nonambiguous_names":len(origu),
       "n_51_nonambiguous_names":len(qa),
       "n_228_missing_or_ambiguous_names":len(b)-len(origu),
       "potential_duplicate_species_from_same_authors":True}
    if len(join)>0:
        q=np.log(join.PTE_228/join.PTE_51)
        overlap["PTE_source_disagreement"]={
          "median_abs_log_ratio":val(np.median(np.abs(q))),
          "median_signed_log_ratio":val(np.median(q)),
          "n_more_than_10pct_relative_disagreement":int((np.abs(np.exp(q)-1)>.1).sum()),
          "different_populations_or_aggregation_possible":True}
        overlap["matches_by_pollen_packaging"]={str(k):int(v) for k,v in join.Poltype.value_counts().items()}
    results["cross_paper_taxon_overlap"]=overlap
    results["no_original_species_rows_or_coordinates_exported"]=True
    (OUT/"johnson_harder_51_species_pollen_ovule_PTE_and_228_overlap_v14.json").write_text(
      json.dumps(val(results),indent=2,ensure_ascii=False)+"\n")
    print("GLOBAL_PO_51_CROSS_SOURCE",json.dumps(val(results),ensure_ascii=False)[:34000],flush=True)
if __name__=="__main__":main()
