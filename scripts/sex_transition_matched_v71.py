#!/usr/bin/env python3
"""v71: source-pinned locality join and genus-held-out species prediction.
Exploratory reanalysis; no individual coordinates, site identifiers or rows exported.
"""
import io, json, hashlib, re, math, concurrent.futures, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import sex_transition_empirical_v70 as old
OUT=Path("analysis_results_sex_transition_v71")
OUT.mkdir(exist_ok=True)
EXPECTED={
"Fig_Share_FemaleFrequency.xlsx":"02e38c74bd71c45553a8359bf4c682a1ccdc4f55dc1d2fb141d73fa209c8a0c6",
"Dataset_S1.xls":"faea1cab456c4165e0ea58ffacd3c81545d1dfcce8272dd98aeb004ad3270458",
"ECE3-16-e74093-s001.xlsx":"0f236915458b6f2fbba49032621e1874b38f27cd5fef62d47a95ec80960db39c"}
BOOKS={}
def capture(name,b,src):
    name=Path(name).name
    if name not in EXPECTED:return {"file":name,"ignored":True}
    sha=hashlib.sha256(b).hexdigest()
    if sha!=EXPECTED[name]:raise ValueError("source hash changed: "+name)
    book=pd.ExcelFile(io.BytesIO(b))
    BOOKS[name]={s:book.parse(s) for s in book.sheet_names}
    return {"file":name,"sha256":sha,"sheets":{s:len(x) for s,x in BOOKS[name].items()}}
old.inspect_file=capture
def fetch_retry(label,fn):
    a=old.job(label,fn)
    if not any(x.get("sha256") for x in a[1]["files"]):
        b=old.job(label,fn)
        if any(x.get("sha256") for x in b[1]["files"]):return b
        a[1]["errors"]+=b[1]["errors"]
    return a
def corr(x,y):
    x=np.asarray(x);y=np.asarray(y)
    if len(x)<3 or len(set(x))<2 or len(set(y))<2:return None
    return float(spearmanr(x,y).statistic)

def cv(d,models,group,alpha=1.0):
    """Equal observation weight; standardize only train data; group holdout."""
    if len(d)<10 or d[group].nunique()<8:return {"status":"underpowered_gate","n":len(d),"groups":int(d[group].nunique())}
    y=d["female_frequency"].to_numpy(float)
    pred={k:np.zeros(len(d)) for k in models}
    gg=d[group].to_numpy()
    for g in sorted(set(gg)):
        train=gg!=g;test=~train
        for k,cols in models.items():
            if not cols:pred[k][test]=y[train].mean();continue
            X=d[cols].to_numpy(float)
            mu=X[train].mean(0);sd=X[train].std(0);sd[sd==0]=1
            a=(X[train]-mu)/sd;z=(X[test]-mu)/sd
            b=np.linalg.solve(a.T@a+alpha*np.eye(len(cols)),a.T@(y[train]-y[train].mean()))
            pred[k][test]=np.clip(y[train].mean()+z@b,0,1)
    out={}
    for k,p in pred.items():
        err=np.abs(y-p)
        out[k]={"MAE":float(err.mean()),"RMSE":float(np.sqrt(np.mean((y-p)**2))),
                "group_equal_MAE":float(np.mean([err[gg==g].mean() for g in set(gg)]))}
    return {"status":"exploratory","n":len(d),"groups":int(d[group].nunique()),"alpha":alpha,"metrics":out}

def macro_tables(m,q):
    d=m["ALL DATA"].copy();d["system"]=d.system.astype(str).str.strip()
    # Varga & Soulsbury 2020, Table 1 and methods define dioecy, gynodioecy,
    # subdioecy; workbook G label corresponds to gynodioecy; not pooled.
    assert set(d.system.unique())=={"D","G","SbD"}
    d=d[d.system=="G"].copy()
    d["sp"]=d.species.map(old.taxon)
    d["freq"]=pd.to_numeric(d.freq,errors="coerce")
    assert d.freq.max()>1 and d.freq.max()<=100 and d.freq.min()>=0
    d["female_frequency"]=d.freq/100
    for c in ["bio9","bio11","lat","long"]:d[c]=pd.to_numeric(d[c],errors="coerce")
    g=d.dropna(subset=["sp","female_frequency","lat","long"])
    # Exact repeated coordinates within species are averaged, not independent sites.
    site=g.groupby(["sp","lat","long"],as_index=False)[["female_frequency","bio9","bio11"]].mean()
    sp=site.groupby("sp").agg(female_frequency=("female_frequency","mean"),bio9=("bio9","mean"),
        bio11=("bio11","mean"),sites=("female_frequency","size"),
        freq_min=("female_frequency","min"),freq_max=("female_frequency","max")).reset_index()
    sp["genus"]=sp.sp.str.split().str[0]
    review=json.loads(Path("data/2026-11-27_sex_function_transition_v68.json").read_text())["rows"][:28]
    prior=pd.DataFrame({"sp":[old.taxon(x["species"]) for x in review],"log_quantity":[math.log(x["A"]) for x in review],"sex_genetics":[x["genetics"] for x in review]})
    v=q["female advantage"].copy()
    v["sp"]=v.species.map(old.taxon);v["yi"]=pd.to_numeric(v.yi,errors="coerce")
    v["vi"]=pd.to_numeric(v.vi,errors="coerce")
    v["stage"]=v["life stage"].astype(str).str.strip().str.lower().str.replace(r"[\s_-]+","",regex=True)
    v=v.dropna(subset=["sp","reference id","yi","vi"]);v=v[v.vi>0]
    ref=v.groupby(["sp","reference id"],as_index=False).yi.mean()
    quality=ref.groupby("sp",as_index=False).agg(quality=("yi","mean"),quality_studies=("reference id","nunique"))
    vs=v.groupby(["sp","reference id","stage"],as_index=False).yi.mean()
    vs=vs.pivot_table(index=["sp","reference id"],columns="stage",values="yi").dropna(subset=["pregermination","postgermination"])
    vs["stage_delta"]=vs["postgermination"]-vs["pregermination"]
    delta=vs.groupby(level="sp")["stage_delta"].mean().reset_index()
    merged=sp.merge(prior,on="sp",how="left").merge(quality,on="sp",how="left").merge(delta,on="sp",how="left")
    report={"G_rows":len(d),"G_exact_species":int(d.sp.nunique()),"valid_rows":len(g),
        "exact_unique_sites":len(site),"species_with_valid_freq_sites":len(sp),
        "species_with_3plus_sites":int((sp.sites>=3).sum()),
        "species_equal_mean_frequency":float(sp.female_frequency.mean()),
        "range_of_species_means":[float(sp.female_frequency.min()),float(sp.female_frequency.max())],
        "quality_outcome_exact_species":int(v.sp.nunique()),"quality_ref_species":len(ref),
        "coverage":{"quantity":int(merged.log_quantity.notna().sum()),"quality":int(merged.quality.notna().sum()),
         "both":int((merged.log_quantity.notna()&merged.quality.notna()).sum()),"paired_stage_delta":int(merged.stage_delta.notna().sum())},
        "scope":"species-level ecological association only: sources may reuse studies, unmatched time/site, no phylogenetic covariance; equal sites then equal species; G subset only"}
    # Same-species panel for each model contrast; do not compare error across panels.
    panels={}
    for label,cols in [("quantity",["log_quantity"]),("quality",["quality"]),("both",["log_quantity","quality"]),("stage",["quality","stage_delta"])]:
        z=merged.dropna(subset=cols+["bio9","bio11"]).copy()
        models={"M0":[], "climate":["bio9","bio11"],"reproductive":cols,"climate_reproductive":["bio9","bio11"]+cols}
        if label=="both":models.update({"climate_quantity":["bio9","bio11","log_quantity"],"climate_quality":["bio9","bio11","quality"]})
        if label=="stage":models["climate_quality"]=["bio9","bio11","quality"]
        panels[label]={"n":len(z),"genera":int(z.genus.nunique()),"species":sorted(z.sp.tolist()),
            "spearman":{c:corr(z[c],z.female_frequency) for c in cols},
            "cv":cv(z,models,"genus")}
        if len(z)>=10 and z.genus.nunique()>=8:
            # Predeclared penalty and leave-one-genus-out; alpha sensitivity.
            panels[label]["alpha_sensitivity"]={str(a):cv(z,models,"genus",a) for a in [.1,10]}
    report["panels"]=panels
    return report

DMS=re.compile(r"(\d{1,3})\s*[°º]\s*(\d{1,2})\s*['′’]\s*(\d+(?:[.,]\d+)?)\s*(?:[\"″”]|'')?\s*([NSEW])",re.I)
def dms(x):
    z={}
    for a,b,c,d in DMS.findall(str(x)):
        a=int(a);b=int(b);c=float(c.replace(",","."))
        if b>=60 or c>=60:return None
        v=a+b/60+c/3600
        z[d.upper()]=v
    if len(z)!=2 or not ({"N","S"}&set(z) and {"E","W"}&set(z)):return None
    return (z["N"] if "N" in z else -z["S"],z["E"] if "E" in z else -z["W"])
def dist(a,b):
    r=math.pi/180
    dlat=(a[0]-b[0])*r;dlon=(a[1]-b[1])*r
    v=math.sin(dlat/2)**2+math.cos(a[0]*r)*math.cos(b[0]*r)*math.sin(dlon/2)**2
    return 6371000*2*math.asin(min(1,math.sqrt(v)))
def census_coord(r):
    try:
        a=float(r.latitude);b=float(r.longitude)
        if math.isfinite(a) and math.isfinite(b) and abs(a)<=90 and abs(b)<=180:return (a,b)
    except (ValueError,TypeError):pass
    return dms(r.coordinates)
def locality_stage(b):
    cols=["Florets per plant","Developed achenes per plant","Viable achenes per plant"]
    b=b.copy();b["_sex"]=b.Sex.map(old.sex);b["_local"]=b.Locality.map(old.key)
    for c in cols:b[c]=pd.to_numeric(b[c],errors="coerce")
    good=b._sex.notna()&b._local.notna()&b[cols].notna().all(1)&(b[cols]>=0).all(1)
    mono=(b[cols[0]]>=b[cols[1]])&(b[cols[1]]>=b[cols[2]])
    return b,cols,good,mono
def stage_rows(b,cols):
    out={};zeros=0;small=0;direct=[]
    for loc,g in b.groupby("_local"):
        f=g[g._sex=="F"];h=g[g._sex=="H"]
        if min(len(f),len(h))<5:small+=1;continue
        F=f[cols].sum().to_numpy(float);H=h[cols].sum().to_numpy(float)
        # Total achene ratio, including zero H denominator, kept separately.
        direct.append({"F_viable":float(F[2]/len(f)),"H_viable":float(H[2]/len(h))})
        if min(F.min(),H.min())<=0:zeros+=1;continue
        vals=[math.log((F[0]/len(f))/(H[0]/len(h))),math.log((F[1]/F[0])/(H[1]/H[0])),
              math.log((F[2]/F[1])/(H[2]/H[1])),math.log((F[2]/len(f))/(H[2]/len(h)))]
        assert abs(sum(vals[:3])-vals[3])<1e-10
        out[loc]=dict(zip(["log_display","log_developed_conversion","log_viability_conversion","log_total"],vals))
    return out,{"n_positive":len(out),"zero_component_sites":zeros,"small_sites":small,
                "zero_H_viable_sites":sum(x["H_viable"]==0 for x in direct),
                "zero_F_viable_sites":sum(x["F_viable"]==0 for x in direct)}
def cirsium_tables(t):
    result={}
    for name,cs,ts in [("C. heterophyllum","Set 1","Set 3"),("C. palustre","Set 2","Set 4")]:
        a=t[cs].copy();raw=t[ts]
        b,cols,good,mono=locality_stage(raw)
        # Explicit coordinate lookup appended below trait rows; only keys present in actual trait groups.
        validloc=set(b.loc[good,"_local"])
        lookup={};conflicts=set()
        for _,r in raw.iterrows():
            k=old.key(r.Locality);p=dms(r.Sex)
            if k not in validloc or p is None:continue
            if k in lookup and lookup[k]!=p:conflicts.add(k)
            lookup[k]=p
        for k in conflicts:lookup.pop(k,None)
        a["_coord"]=a.apply(census_coord,axis=1)
        a["_F"]=pd.to_numeric(a["Number of females"],errors="coerce")
        a["_H"]=pd.to_numeric(a["Number of cosexuals"],errors="coerce")
        a["_N"]=pd.to_numeric(a["population_sample_size"],errors="coerce")
        a=a[a[["_F","_H","_N"]].notna().all(1)&(a._F>=0)&(a._H>=0)&(a._N>0)&(a._F+a._H==a._N)&a._coord.notna()]
        records=[]
        for key,p in lookup.items():
            ds=sorted((dist(p,r._coord),i) for i,r in a.iterrows())
            if ds:records.append((key,ds))
        # One and only one candidate within threshold, no repeated census mapping.
        def match(threshold):
            pairs={k:ds[0][1] for k,ds in records if ds[0][0]<=threshold and (len(ds)==1 or ds[1][0]>threshold)}
            counts=pd.Series(list(pairs.values()),dtype=object).value_counts()
            pairs={k:i for k,i in pairs.items() if counts[i]==1}
            return pairs
        pairs=match(250)
        res={"trait_rows_valid":int(good.sum()),"trait_rows_nonmonotonic":int((good&~mono).sum()),
             "valid_trait_localities":len(validloc),"coordinate_lookup_localities":len(lookup),
             "lookup_conflicts":len(conflicts),"valid_census_rows":len(a),
             "unique_matches_by_max_m":{str(x):len(match(x)) for x in [25,100,250,1000]},
             "primary_match_max_m":250,
             "nearest_distance_quantiles_m":np.quantile([ds[0][0] for _,ds in records],[0,.5,1]).tolist() if records else None,
             "caution":"coordinates used internally only; censuses and traits may be different years; census sex ratio independent of trait sampling; source estimated per-plant quantities; no causal inference"}
        res["analyses"]={}
        for label,mask in [("all_complete",good),("monotonic_only",good&mono)]:
            stage,stats=stage_rows(b.loc[mask],cols)
            rows=[]
            for loc,s in stage.items():
                if loc not in pairs:continue
                r=a.loc[pairs[loc]]
                rows.append({**s,"female_frequency":float(r._F/r._N)})
            z=pd.DataFrame(rows)
            stats["matched_positive_localities"]=len(z)
            if stage:
                d=pd.DataFrame(stage.values())
                stats["geometric_ratios"]={c:float(np.exp(d[c].mean())) for c in d.columns}
                stats["display_lower_total_higher"]=int(((d.log_display<0)&(d.log_total>0)).sum())
                if len(d)>=3:
                    rng=np.random.default_rng(711)
                    ar=d.to_numpy();bs=ar[rng.integers(len(ar),size=(9999,len(ar)))].mean(1)
                    stats["locality_bootstrap_intervals"]={c:np.exp(np.quantile(bs[:,i],[.025,.975])).tolist() for i,c in enumerate(d.columns)}
            if len(z)>=8:
                z["locality_holdout"]=np.arange(len(z))
                stats["correlations"]={c:corr(z[c],z.female_frequency) for c in ["log_total","log_display","log_developed_conversion","log_viability_conversion"]}
                models={"M0":[],"M1_total":["log_total"],"M2_total_plus_stage":["log_total","log_display","log_developed_conversion"]}
                stats["cv"]=cv(z,models,"locality_holdout")
            res["analyses"][label]=stats
        result[name]=res
    return result

def main():
    assert dms("49°27'48.2\"N 16°01'36.0\"E") is not None
    assert abs(dist((50,15),(50,15)))<1e-8
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        sources=dict(ex.map(lambda x:fetch_retry(*x),[("figshare",old.figshare),("quality",old.dryad),("cirsium",old.cirsium)]))
    # Only source-level summaries; no original row, coords or site name.
    res={"sources":{k:{"files":v["files"],"error_count":len(v["errors"])} for k,v in sources.items()}}
    for k,v in res["sources"].items():print("SOURCE",k,json.dumps(v,ensure_ascii=False),flush=True)
    errors={}
    tasks=[
        ("macro",{"Fig_Share_FemaleFrequency.xlsx","Dataset_S1.xls"},lambda:macro_tables(BOOKS["Fig_Share_FemaleFrequency.xlsx"],BOOKS["Dataset_S1.xls"])),
        ("cirsium",{"ECE3-16-e74093-s001.xlsx"},lambda:cirsium_tables(BOOKS["ECE3-16-e74093-s001.xlsx"]))
    ]
    for k,required,fn in tasks:
        if required<=BOOKS.keys():
            try:
                value=fn()
                message=json.dumps(value,ensure_ascii=False,allow_nan=False)
                res[k]=value
                print("RESULT",k,message,flush=True)
            except Exception as e:
                import traceback
                errors[k]=type(e).__name__+": "+str(e)
                traceback.print_exc()
        else:errors[k]="missing hash-pinned inputs"
        res["analysis_errors"]=errors
        (OUT/"summary.json").write_text(json.dumps(res,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    if errors:raise RuntimeError("analysis gates incomplete: "+str(errors))
if __name__=="__main__":main()
