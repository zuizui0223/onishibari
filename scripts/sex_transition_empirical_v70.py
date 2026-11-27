#!/usr/bin/env python3
"""v69: bounded, public, read-only source acquisition and coverage audit.
Do not emit individual records/coordinates. No credential or access bypass.
Schema discovery is not evidence of sex system or matching cohorts.
"""
from pathlib import Path
import io, json, hashlib, re, zipfile, urllib.request, urllib.parse, concurrent.futures
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

OUT=Path("analysis_results_sex_transition_v70")
OUT.mkdir(exist_ok=True)
MAX=20_000_000
def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"Scientific-reproducibility-audit/0.1"})
    with urllib.request.urlopen(req,timeout=30) as r:
        b=r.read(MAX+1)
        if len(b)>MAX: raise ValueError("download cap")
        return b
def jget(url): return json.loads(get(url))
def norm(x): return re.sub(r"[^a-z0-9]","",str(x).lower())
def taxon(x):
    s=re.sub(r"[_\s]+"," ",str(x).strip()).lower()
    return s if re.fullmatch(r"[a-z][a-z-]+ [a-z][a-z-]+",s) else None
def inspect_file(name,b,src):
    suffix=Path(name).suffix.lower()
    rec={"file":name,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"sheets":[]}
    tables=[]
    if suffix in (".xls",".xlsx",".xlsm"):
        book=pd.ExcelFile(io.BytesIO(b))
        tables=[(s,book.parse(s)) for s in book.sheet_names]
    elif suffix in (".csv",".tsv",".txt"):
        try:
            df=pd.read_csv(io.BytesIO(b),sep=None,engine="python",on_bad_lines="error",encoding="utf-8-sig")
            if len(df.columns)>1: tables=[("text",df)]
        except Exception:
            pass
    if name.endswith("Dataset_S1.xls"):
        src["quality_paired_analysis"]=quality_analysis(dict(tables))
    if name.endswith("ECE3-16-e74093-s001.xlsx"):
        src["cirsium_population_analysis"]=cirsium_analysis(dict(tables))
    for sheet,df in tables:
        cols=[str(c) for c in df.columns]
        row={"sheet":sheet,"rows":len(df),"columns":cols}
        # Only explicitly recognized binomial column; no automatic genus-only or synonym match.
        speciescols=[c for c in df.columns if norm(c) in {"species","speciesname","taxon","scientificname","specieslatinname","speciesid","speciesnames","specie"}]
        if len(speciescols)==1:
            names=set(filter(None,(taxon(x) for x in df[speciescols[0]].dropna())))
            row["binomial_count"]=len(names);row["binomial_column"]=str(speciescols[0])
            src.setdefault("sets",{})[name+"::"+sheet]=sorted(names)
            # Only inspect distribution of sexual system if column explicitly named.
            syscols=[c for c in df.columns if norm(c) in {"sexualsystem","sexsystem","breedingsystem","system","sexualsystems"}]
            if len(syscols)==1:
                sc=syscols[0];row["system_counts"]={str(k):int(v) for k,v in df[sc].value_counts(dropna=False).items()}
                mask=df[sc].astype(str).str.strip().str.lower().isin(["gynodioecy","gynodioecious"])
                src.setdefault("gynodioecy_sets",{})[name+"::"+sheet]=sorted(set(filter(None,(taxon(x) for x in df.loc[mask,speciescols[0]].dropna()))))
                row["gynodioecy_rows_exact_text"]=int(mask.sum())
        rec["sheets"].append(row)
    return rec
def inspect(name,b,src):
    # Never extract archive contents onto disk. Public source integrity only.
    if zipfile.is_zipfile(io.BytesIO(b)) and Path(name).suffix.lower() not in (".xlsx",".xlsm"):
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            if sum(i.file_size for i in z.infolist())>MAX*3: raise ValueError("archive cap")
            for i in z.infolist():
                if i.is_dir():continue
                if Path(i.filename).suffix.lower() not in (".xlsx",".xls",".xlsm",".csv",".tsv",".txt"):continue
                try:src["files"].append(inspect_file(i.filename,z.read(i),src))
                except Exception as e:src["errors"].append(i.filename+": "+str(e)[:250])
    else:
        src["files"].append(inspect_file(name,b,src))

def quality_analysis(tables):
    d=tables["female advantage"].copy()
    d["stage"]=d["life stage"].astype(str).str.strip().str.lower().str.replace(r"[\s_-]+","",regex=True)
    counts=d["stage"].value_counts(dropna=False).to_dict()
    d["yi"]=pd.to_numeric(d["yi"],errors="coerce")
    d["vi"]=pd.to_numeric(d["vi"],errors="coerce")
    d=d.dropna(subset=["reference id","species","yi","vi"])
    d=d[(d.vi>0)&d.stage.isin(["pregermination","postgermination"])]
    result={"stage_counts":counts,"valid_records":len(d),"estimand":"within reference x species, post minus pre standardized female advantage; not paired individuals, not a log RR; equal study weight; cluster bootstrap by reference"}
    if not len(d):return result
    pairs=[]
    for (ref,sp),g in d.groupby(["reference id","species"]):
        a=g[g.stage=="pregermination"];b=g[g.stage=="postgermination"]
        if len(a) and len(b):
            pairs.append({"ref":str(ref),"species":str(sp),"pre_n":len(a),"post_n":len(b),"diff":float(b.yi.mean()-a.yi.mean()),"ivdiff":float(np.average(b.yi,weights=1/b.vi)-np.average(a.yi,weights=1/a.vi))})
    result["pairs"]=len(pairs)
    if not pairs:return result
    p=pd.DataFrame(pairs)
    study=p.groupby("ref")[["diff","ivdiff"]].mean()
    result.update({"studies":len(study),"species":p.species.nunique(),"species_study_pairs":pairs,"mean_equal_study":float(study["diff"].mean()),"inverse_variance_within_stage_sensitivity":float(study["ivdiff"].mean()),"positive_studies":int((study["diff"]>0).sum())})
    if len(study)>=3:
        rng=np.random.default_rng(20261127)
        a=study["diff"].to_numpy()
        bs=rng.choice(a,(10000,len(a)),replace=True).mean(axis=1)
        result["descriptive_study_bootstrap_95"]=np.quantile(bs,[.025,.975]).tolist()
        result["leave_one_study_out_range"]=[float((a.sum()-a).min()/(len(a)-1)),float((a.sum()-a).max()/(len(a)-1))]
    return result

def key(v):
    if pd.isna(v):return None
    if isinstance(v,(float,np.floating)) and v.is_integer():return str(int(v))
    return re.sub(r"\s+"," ",str(v).strip().casefold()) or None
def sex(v):
    z=key(v)
    if z in {"female","f","females"}:return "F"
    if z in {"hermaphrodite","hermaphrodites","hermaphroditic","h","cosexual","cosexuals","c","co-sexual","co-sexuals","bisexual","b"}:return "H"
    return None
def loocv(rows):
    # Population-equal-weight, held-out locality, ridge penalty fixed at 1.
    # The contrast predicts an independent outcome, not a quantity used in its definition.
    q=np.array([r["female_frequency"] for r in rows])
    X=np.array([[r["log_total"],r["log_display"],r["log_developed_conversion"],r["log_viability_conversion"]] for r in rows])
    pred={"intercept":[],"total":[],"total_plus_stage":[]}
    for i in range(len(rows)):
        train=np.arange(len(rows))!=i
        pred["intercept"].append(float(q[train].mean()))
        # Two independent log-stage coordinates plus log total; the omitted third coordinate is determined by the total.
        for label,cols in [("total",[0]),("total_plus_stage",[0,1,2])]:
            a=X[train][:,cols];z=X[i,cols];m=a.mean(axis=0);sd=a.std(axis=0);sd[sd==0]=1
            a=(a-m)/sd;z=(z-m)/sd
            beta=np.linalg.solve(a.T@a+np.eye(len(cols)),a.T@(q[train]-q[train].mean()))
            pred[label].append(float(np.clip(q[train].mean()+z@beta,0,1)))
    return {k:{"MAE":float(np.abs(np.asarray(v)-q).mean()),"RMSE":float(np.sqrt(((np.asarray(v)-q)**2).mean()))} for k,v in pred.items()}

def cirsium_analysis(tables):
    allres={}
    for name,census,traits in [("C. heterophyllum","Set 1","Set 3"),("C. palustre","Set 2","Set 4")]:
        a=tables[census].copy();b=tables[traits].copy()
        b["_sex"]=b["Sex"].map(sex);b["_local"]=b["Locality"].map(key)
        cols=["Florets per plant","Developed achenes per plant","Viable achenes per plant"]
        for c in cols:b[c]=pd.to_numeric(b[c],errors="coerce")
        complete=b["_sex"].notna()&b["_local"].notna()&b[cols].notna().all(axis=1)&(b[cols]>=0).all(axis=1)
        monotonic=(b[cols[0]]>=b[cols[1]])&(b[cols[1]]>=b[cols[2]])
        result={"census_rows":len(a),"trait_rows":len(b),"sex_labels":{str(k):int(v) for k,v in b.Sex.value_counts(dropna=False).items()},"complete_nonnegative_rows":int(complete.sum()),"nonmonotonic_complete_rows":int((complete&~monotonic).sum()),"localities_complete":b.loc[complete,"_local"].nunique(),"site_matching":"exact normalized code or locality description, only unique census row; no fuzzy geographic matching; survey years/cohort equality not established"}
        # The identities telescope even for non-monotonic quantities, but such rows may indicate
        # mismatched measurement or units; exclude and report instead of silently truncating.
        b=b.loc[complete&monotonic].copy()
        a["_F"]=pd.to_numeric(a["Number of females"],errors="coerce")
        a["_H"]=pd.to_numeric(a["Number of cosexuals"],errors="coerce")
        a["_N"]=pd.to_numeric(a["population_sample_size"],errors="coerce")
        a=a[a[["_F","_H","_N"]].notna().all(axis=1)&(a._F>=0)&(a._H>=0)&(a._N>0)&(a._F+a._H==a._N)]
        maps={}
        for _,r in a.iterrows():
            for col in ["code"," locality description"]:
                k=key(r[col])
                if k is not None:maps.setdefault(k,[]).append(r.name)
        rows=[]; matched=[];unmatched=0;ambiguous=0;nonpositive=0
        for loc,g in b.groupby("_local"):
            F=g[g._sex=="F"];H=g[g._sex=="H"]
            if min(len(F),len(H))<5:continue
            f=F[cols].sum().to_numpy(dtype=float);h=H[cols].sum().to_numpy(dtype=float)
            if min(f.min(),h.min())<=0:
                nonpositive+=1;continue
            l_display=np.log((f[0]/len(F))/(h[0]/len(H)))
            l_dev=np.log((f[1]/f[0])/(h[1]/h[0]))
            l_v=np.log((f[2]/f[1])/(h[2]/h[1]))
            l_tot=np.log((f[2]/len(F))/(h[2]/len(H)))
            assert abs(l_display+l_dev+l_v-l_tot)<1e-10
            z={"nF":len(F),"nH":len(H),"log_display":float(l_display),"log_developed_conversion":float(l_dev),"log_viability_conversion":float(l_v),"log_total":float(l_tot)}
            rows.append(z)
            ids=set(maps.get(loc,[]))
            if len(ids)==1:
                r=a.loc[next(iter(ids))]
                # Do not output the locality identifiers or coordinates.
                matched.append({**z,"female_frequency":float(r._F/r._N),"census_n":int(r._N)})
            elif not ids:unmatched+=1
            else:ambiguous+=1
        result.update({"valid_census_rows":len(a),"trait_localities_with_positive_components":len(rows),"zero_component_localities":nonpositive,"matched_localities":len(matched),"unmatched_localities":unmatched,"ambiguous_localities":ambiguous})
        if rows:
            d=pd.DataFrame(rows)
            result["locality_equal_weight_geometric_ratios"]={c:float(np.exp(d[c].mean())) for c in ["log_display","log_developed_conversion","log_viability_conversion","log_total"]}
            result["display_lower_but_total_higher_localities"]=int(((d.log_display<0)&(d.log_total>0)).sum())
        if len(matched)>=8:
            d=pd.DataFrame(matched)
            result["correlations"]={c:{"rho":float(spearmanr(d[c],d.female_frequency).statistic),"p_unadjusted":float(spearmanr(d[c],d.female_frequency).pvalue)} for c in ["log_total","log_display","log_developed_conversion","log_viability_conversion"] if d[c].nunique()>1 and d.female_frequency.nunique()>1}
            result["leave_one_locality_out"]=loocv(matched)
        result["matched_summary"]=matched
        allres[name]=result
    return allres


def figshare(src):
    j=jget("https://api.figshare.com/v2/articles/12280298")
    src["version"]=j.get("version");src["title"]=j.get("title");src["license"]=j.get("license")
    for f in j.get("files",[]):
        if f.get("size",MAX+1)>MAX:continue
        try:inspect(f["name"],get(f.get("downloadUrl") or f.get("download_url") or ("https://api.figshare.com/v2/file/download/"+str(f["id"]))),src)
        except Exception as e:src["errors"].append(f["name"]+": "+str(e)[:250])
def dryad(src):
    base="https://datadryad.org/api/v2"
    doi="10.5061/dryad.w9ghx3g62"
    url=base+"/datasets/"+urllib.parse.quote("doi:"+doi,safe="")
    j=jget(url+"/versions")
    versions=j.get("_embedded",{}).get("stash:versions",[])
    if not versions:raise ValueError("no versions")
    v=max(versions,key=lambda x:x.get("versionNumber",0))
    href=v.get("_links",{}).get("self",{}).get("href","")
    m=re.search(r"/versions/(\d+)",href)
    vid=m.group(1) if m else v.get("id")
    if not vid:raise ValueError("no version identifier")
    src["version_id"]=vid
    try:
        inspect("dryad.zip",get(base+"/versions/"+str(vid)+"/download"),src)
    except Exception as e:
        src["errors"].append("version archive: "+str(e)[:250])
    if not src["files"]:
        fs=jget(base+"/versions/"+str(vid)+"/files")
        for f in fs.get("_embedded",{}).get("stash:files",[]):
            name=f.get("path") or f.get("filename") or f.get("name") or ""
            if not name.lower().endswith((".xls",".xlsx",".csv",".tsv")):continue
            match=re.search(r"/files/(\d+)",f.get("_links",{}).get("self",{}).get("href",""))
            fid=f.get("id") or (match.group(1) if match else None)
            if not fid:continue
            try:inspect(name,get(base+"/files/"+str(fid)+"/download"),src)
            except Exception as e:src["errors"].append(name+": "+str(e)[:250])
def cirsium(src):
    b=get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13400990/supplementaryFiles")
    if not zipfile.is_zipfile(io.BytesIO(b)):raise ValueError("no supplementary zip")
    inspect("cirsium_supp.zip",b,src)
def job(label,fn):
    src={"files":[],"errors":[]}
    try:fn(src)
    except Exception as e:src["errors"].append(str(e)[:300])
    return label,src
def main():
    assert taxon("Cirsium palustre")=="cirsium palustre"
    assert taxon("Cirsium_heterophyllum")=="cirsium heterophyllum"
    assert taxon("Daphne jezoensis")=="daphne jezoensis"
    assert taxon("Cirsium sp. 1") is None
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        res=dict(ex.map(lambda x:job(*x),[("figshare",figshare),("dryad_quality",dryad),("cirsium",cirsium)]))
    q=set().union(*map(set,res["dryad_quality"].get("sets",{}).values()))
    macro=set().union(*map(set,res["figshare"].get("sets",{}).values()))
    gd=set().union(*map(set,res["figshare"].get("gynodioecy_sets",{}).values()))
    anchors=json.loads(Path("data/2026-11-27_sex_function_transition_v68.json").read_text())["rows"][:28]
    prior={taxon(x["species"]) for x in anchors}
    report={"sources":res,"coverage":{"quality_exact_binomials":len(q),"macro_all_exact_binomials":len(macro),"macro_gynodioecy_explicit_text_exact_binomials":len(gd),"quality_macro_all_overlap":len(q&macro) if q and macro else None,"quality_macro_gynodioecy_overlap":len(q&gd) if q and gd else None,"review_quality_overlap":len(prior&q) if q else None,"review_macro_gynodioecy_overlap":len(prior&gd) if gd else None,"triple_overlap":sorted(q&gd&prior) if q and gd else None},"caution":"Zero matched names without source access or recognized headers is NOT evidence of biological non-overlap. No coordinate/site matching, no meta-analysis. Inspect original codebook before interpreting system codes."}
    (OUT/"source_coverage.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    for k,v in res.items():
        for name in ["quality_paired_analysis","cirsium_population_analysis"]:
            if name in v:
                obj=v[name]
                if name=="quality_paired_analysis":obj={k:v for k,v in obj.items() if k!="species_study_pairs"}
                else:obj={k:{x:y for x,y in v.items() if x!="matched_summary"} for k,v in obj.items()}
                print("ANALYSIS",name,json.dumps(obj,ensure_ascii=False),flush=True)
        print("SOURCE",k,"FILES",len(v["files"]),"ERRORS",json.dumps(v["errors"],ensure_ascii=False),flush=True)
        for f in v["files"]:print("FILE",f["file"],f["sha256"],json.dumps(f["sheets"],ensure_ascii=False)[:8000],flush=True)
    print("COVERAGE",json.dumps(report["coverage"],ensure_ascii=False),flush=True)
if __name__=="__main__":main()
