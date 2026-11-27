#!/usr/bin/env python3
"""v69: bounded, public, read-only source acquisition and coverage audit.
Do not emit individual records/coordinates. No credential or access bypass.
Schema discovery is not evidence of sex system or matching cohorts.
"""
from pathlib import Path
import io, json, hashlib, re, zipfile, urllib.request, urllib.parse, concurrent.futures
import pandas as pd

OUT=Path("analysis_results_sex_transition_v69")
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
    for sheet,df in tables:
        cols=[str(c) for c in df.columns]
        row={"sheet":sheet,"rows":len(df),"columns":cols}
        # Only explicitly recognized binomial column; no automatic genus-only or synonym match.
        speciescols=[c for c in df.columns if norm(c) in {"species","speciesname","taxon","scientificname","specieslatinname","speciesid","speciesnames"}]
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
def figshare(src):
    j=jget("https://api.figshare.com/v2/articles/12280298")
    src["version"]=j.get("version");src["title"]=j.get("title");src["license"]=j.get("license")
    for f in j.get("files",[]):
        if f.get("size",MAX+1)>MAX:continue
        try:inspect(f["name"],get(f["downloadUrl"]),src)
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
        print("SOURCE",k,"FILES",len(v["files"]),"ERRORS",json.dumps(v["errors"],ensure_ascii=False),flush=True)
        for f in v["files"]:print("FILE",f["file"],json.dumps(f["sheets"],ensure_ascii=False)[:8000],flush=True)
    print("COVERAGE",json.dumps(report["coverage"],ensure_ascii=False),flush=True)
if __name__=="__main__":main()
