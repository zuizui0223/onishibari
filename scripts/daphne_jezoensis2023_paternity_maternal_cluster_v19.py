#!/usr/bin/env python3
"""v0.19: source-original, family-clustered pollen dispersal distances in Daphne jezoensis.
Shibata & Kudo 2023, AJB DOI 10.1002/ajb2.16239, Zenodo 8253921.

Reconstruct only high-confidence assigned maternal/paternal genotypes as in author
Figures_Analyses.R; recalculate distances from adult source East.g/North.g without
outputting locations or genotype IDs. Nested seeds are NOT independent mothers.
Any year×mother comparisons use matched pollen recipient families and avoid
extrapolating to unassigned seeds. Published experiment, not independent replication.
"""
import json,io,hashlib,urllib.request,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path("analysis_results_daphne_paternity_family_v19")
API="https://zenodo.org/api/records/8253921"
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
SOURCE_SHA={
    "NF2018_2_primer16_res1.2.3_var1_150000_220702.csv":"fefd5e4a31880fbf1c1e70f55b897e0ac37e78837ebd40408b7a851f94016d19",
    "NF2020_2_primer16_res1.2.3_var1_150000_220702.csv":"957c293549d54b3d3a1ef5683359091c207844450215ddee11580cc23faa1d78",
    "Sample_DATA.201820_seed.csv":"739dafbbee07de46ee3c6a10086b754f1e33b697d111a9ec8ee22d5a8a35d9b1",
    "GenetID_geno_phenotype_1820Nopadult_211109.csv":"97790d8b3467d2f9e860f97607aa60525a17ddf1b97b0645fdb1b179ce99135e"
}
RNG=np.random.default_rng(20261116)

def get(url,limit=2000000):
    req=urllib.request.Request(url,headers={"User-Agent":"published-Daphne-parentage-secondary-audit-0.19","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as r:b=r.read(limit+1)
    if len(b)>limit:raise RuntimeError("source exceeds byte cap")
    return b
def load():
    metadata=json.loads(get(API,2000000))
    fs=metadata.get("files")
    entries=fs.get("entries",{}) if isinstance(fs,dict) else {f.get("key"):f for f in fs}
    zipurl=entries["Data for Figures_Analyses .zip"]["links"].get("content") or entries["Data for Figures_Analyses .zip"]["links"]["self"]
    archive=get(zipurl)
    assert hashlib.sha256(archive).hexdigest()==ZIP_SHA
    results={}
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        for n in z.namelist():
            k=Path(n).name
            if k in SOURCE_SHA:
                b=z.read(n)
                if hashlib.sha256(b).hexdigest()!=SOURCE_SHA[k]:
                    raise ValueError("source SHA mismatch: "+k)
                results[k]=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig")
    assert set(results)==set(SOURCE_SHA)
    return results

def validate_and_join(files):
    fn18="NF2018_2_primer16_res1.2.3_var1_150000_220702.csv"
    fn20="NF2020_2_primer16_res1.2.3_var1_150000_220702.csv"
    d18=files[fn18].copy();d18["Year"]=2018
    d20=files[fn20].copy();d20["Year"]=2020
    assert len(d18)==257 and len(d20)==327
    allseed=pd.concat([d18,d20],ignore_index=True)
    assert len(allseed)==584 and {"DNAID","Mother","Father","Probability"}<=set(allseed)
    samples=files["Sample_DATA.201820_seed.csv"]
    assert {"DNAID","Info"}<=set(samples)
    cross=samples[["DNAID","Info"]].dropna(subset=["DNAID"])
    conflicted=cross.groupby("DNAID")["Info"].nunique(dropna=False)
    if (conflicted>1).any():
        raise ValueError("source DNAID -> Info is non-unique with different source treatments")
    cross=cross.drop_duplicates("DNAID")
    d=allseed.merge(cross,on="DNAID",how="left",validate="m:1")
    d["Probability"]=pd.to_numeric(d["Probability"],errors="coerce")
    assigned=d[d.Father.notna() & d.Probability.ge(.95)].copy()
    # Do not coerce unknown fathers to new 'absence of pollen' outcomes.
    d["assigned_95"]=d.Father.notna() & d.Probability.ge(.95)
    return d,assigned,{"n_all_seeds":len(d),"n_assigned_95":len(assigned),
        "n_samples_source":len(samples),"n_unmatched_source_treatment":int(d.Info.isna().sum()),
        "n_dnaid_duplicate_rows":int(samples.DNAID.duplicated().sum()),
        "treatment_labels":{str(k):int(v) for k,v in d.Info.value_counts(dropna=False).items()},
        "confident_treatment_labels":{str(k):int(v) for k,v in assigned.Info.value_counts(dropna=False).items()}}

def distances(assigned,files):
    adult=files["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    cols=["GenetID2","East.g","North.g"]
    assert set(cols)<=set(adult)
    g=adult[cols].dropna(subset=["GenetID2"]).copy()
    for c in ["East.g","North.g"]:g[c]=pd.to_numeric(g[c],errors="coerce")
    g=g.dropna(subset=["East.g","North.g"])
    spread=g.groupby("GenetID2").agg(
        minx=("East.g","min"),maxx=("East.g","max"),
        miny=("North.g","min"),maxy=("North.g","max"))
    max_disagreement=np.sqrt((spread.maxx-spread.minx)**2+(spread.maxy-spread.miny)**2)
    diag={"n_adult_rows":len(adult),"n_coordinate_rows":len(g),
          "n_unique_coordinate_genets":g.GenetID2.nunique(),
          "n_genets_with_1m_coordinate_discrepancy":int((max_disagreement>1).sum()),
          "n_genets_with_any_coordinate_discrepancy":int((max_disagreement>1e-8).sum())}
    # Use first coordinate like authors' distinct-per-genet route, but refuse
    # a disagreement at the >= 1 metre threshold.
    if diag["n_genets_with_1m_coordinate_discrepancy"]:
        raise ValueError("same source genotype has conflicting locations: cannot assign distance unambiguously")
    unique=g.drop_duplicates("GenetID2")
    recipient=unique.rename(columns={"GenetID2":"Mother","East.g":"mx","North.g":"my"})
    donor=unique.rename(columns={"GenetID2":"Father","East.g":"fx","North.g":"fy"})
    d=assigned.merge(recipient,on="Mother",how="left",validate="m:1")
    d=d.merge(donor,on="Father",how="left",validate="m:1")
    have=d[["mx","my","fx","fy"]].notna().all(axis=1)
    diag["n_95_assigned_with_both_coordinates"]=int(have.sum())
    diag["n_95_assigned_without_both_coordinates"]=int((~have).sum())
    d=d[have].copy()
    d["distance_m"]=np.hypot(d.mx-d.fx,d.my-d.fy)
    d["within_10m"]=d.distance_m.lt(10)
    assert d.distance_m.notna().all()
    return d,diag

def mother_cluster_boot(x):
    fam=[df for _,df in x.groupby("Mother")]
    n=len(fam)
    if n<3:return None
    seedfracs=np.array([v.within_10m.mean() for v in fam])
    counts=np.array([len(v) for v in fam],dtype=float)
    nears=np.array([v.within_10m.sum() for v in fam],dtype=float)
    ix=RNG.integers(0,n,size=(8000,n))
    mf=seedfracs[ix].mean(axis=1)
    sf=nears[ix].sum(axis=1)/counts[ix].sum(axis=1)
    return {"mother_equal_CI95":[float(a) for a in np.quantile(mf,[.025,.975])],
            "offspring_equal_CI95":[float(a) for a in np.quantile(sf,[.025,.975])]}

def summarize(d):
    out=[]
    for (year,info),x in d.groupby(["Year","Info"],dropna=False):
        groups=x.groupby("Mother")
        mom=groups.within_10m.mean()
        v={"year":int(year),"info":str(info),
          "assigned_seed_with_distance":len(x),
          "distinct_mothers":int(x.Mother.nunique()),
          "distinct_fathers":int(x.Father.nunique()),
          "offspring_equal_frac_within_10m":float(x.within_10m.mean()),
          "mother_equal_frac_within_10m":float(mom.mean()),
          "mean_distance_m_offspring_equal":float(x.distance_m.mean()),
          "median_distance_m_offspring_equal":float(x.distance_m.median()),
          "maternal_family_sensitivity":mother_cluster_boot(x)}
        out.append(v)
    return sorted(out,key=lambda x:(x["year"],x["info"]))

def matched_mothers(d):
    a=d[d.Info.isin(["Day open","Night open"])].groupby(
       ["Year","Mother","Info"],as_index=False).agg(
        offspring=("within_10m","size"),near=("within_10m","sum"))
    a["near_fraction"]=a.near/a.offspring
    p=a.pivot(index=["Year","Mother"],columns="Info",values=["offspring","near_fraction"]).dropna()
    p.columns=["_".join(c) for c in p.columns]
    p=p.reset_index()
    p["diff_night_day"]=p["near_fraction_Night open"]-p["near_fraction_Day open"]
    out=[]
    for y,z in p.groupby("Year"):
        x=z.diff_night_day.to_numpy()
        n=len(x)
        estimates=np.mean(x[RNG.integers(0,n,size=(8000,n))],axis=1) if n>=3 else None
        out.append({"year":int(y),"matched_mothers":n,
            "mean_mother_equal_diff_night_minus_day":float(np.mean(x)),
            "n_night_higher":int((x>0).sum()),"n_day_higher":int((x<0).sum()),"n_equal":int((x==0).sum()),
            "mother_pair_boot_CI95":[float(v) for v in np.quantile(estimates,[.025,.975])] if estimates is not None else None,
            "note":"only mothers with >=1 high-confidence assigned offspring in BOTH windows; not random sample"})
    return out

def main():
    OUT.mkdir(exist_ok=True)
    f=load()
    allseed,assigned,joinqc=validate_and_join(f)
    parent,coordqc=distances(assigned,f)
    stat=summarize(parent)
    paired=matched_mothers(parent)
    out={"status":"original_zenodo_2023_paternity_secondary_finished",
      "source_article_doi":"10.1002/ajb2.16239",
      "source_data_doi":"10.5281/zenodo.8253921",
      "zip_sha256":ZIP_SHA,
      "original_file_sha256":SOURCE_SHA,
      "joint_source_integrity":joinqc,
      "coordinate_quality":coordqc,
      "per_group":stat,
      "paired_day_night_maternal_families":paired,
      "comparisons_require_assignment_and_known_coordinates":True,
      "not_all_fruit_or_pollen_received_produce_genetically_assigned_seeds":True,
      "not_a_direct_fitness_mediation_from_v18_fruit_stages":True,
      "no_anonymized_geolocations_genotypes_or_plant_ids_released":True}
    (OUT/"shibata_kudo2023_paternity_maternal_cluster_v19.json").write_text(json.dumps(out,indent=2)+"\n")
    # Final output contains aggregated numbers only.
    print("DAPHNE_PATERNITY_FAMILY_V19",json.dumps(out),flush=True)

if __name__=="__main__":main()
