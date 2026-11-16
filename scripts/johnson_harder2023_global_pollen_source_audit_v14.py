#!/usr/bin/env python3
"""Audit Johnson & Harder 2023 public ORIGINAL species-level pollen economy files.

Two author sources are inspected, not assumed equivalent:
  Proc R Soc B: Johnson and Harder (2023) 228 spp., 10.1098/rspb.2023.1148,
  supplementary DOI 10.6084/m9.figshare.c.6837249;
  AJB: Harder and Johnson (2023) pollen:ovule data and PTE,
  Dryad DOI 10.5061/dryad.w6m905qtq.
Official Figshare /v2/collections and Dryad /api/v2 published archives ONLY.
Produces manifests, column names, small redacted category / data-type summaries;
not full original species records, locality, or unlicensed republishing.
"""
from __future__ import annotations
import csv,hashlib,io,json,re,urllib.request,urllib.error,urllib.parse,zipfile
from pathlib import Path
import pandas as pd

OUT=Path("analysis_results_global_pollen_economy_original_v14")
OUT.mkdir(exist_ok=True)
FIGSHARE_COLLECTION=6837249
DRYAD_DOI="10.5061/dryad.w6m905qtq"
CAP=10_000_000

def pull(url,cap=CAP):
    req=urllib.request.Request(url,headers={
      "Accept":"application/json, application/octet-stream, application/zip",
      "User-Agent":"Open-pollen-economy-source-2026-11-16/0.1"})
    try:
        with urllib.request.urlopen(req,timeout=70) as r:
            b=r.read(cap+1)
            if len(b)>cap:return {"status":int(r.status),"error":"above_byte_limit"}
            return {"status":int(r.status),"body":b,"final_url":r.geturl(),
                    "content_type":r.headers.get("Content-Type")}
    except urllib.error.HTTPError as e:return {"status":int(e.code),"error":"HTTP_error"}
    except Exception as e:return {"error":type(e).__name__+":"+str(e)[:200]}

def js(o):
    if o.get("status")!=200:return None
    try:return json.loads(o["body"])
    except (ValueError,TypeError,UnicodeError):return None

def report_table(name,b):
    out={"name":name,"bytes":len(b),"sha256_raw":hashlib.sha256(b).hexdigest()}
    name=name.lower()
    if name.endswith((".csv",".tsv",".txt",".xlsx",".xls")):
        if name.endswith((".xlsx",".xls")):
            try:
                xls=pd.ExcelFile(io.BytesIO(b))
                out["worksheets"]={}
                for sheet in xls.sheet_names[:15]:
                    data=pd.read_excel(xls,sheet_name=sheet)
                    out["worksheets"][str(sheet)]=profile(data)
            except Exception as e:out["parse_error"]=type(e).__name__+":"+str(e)[:200]
        elif name.endswith((".csv",".tsv",".txt")):
            for encoding in ("utf-8-sig","latin-1"):
                try:
                    txt=b.decode(encoding)
                    sep="\t" if name.endswith(".tsv") else None
                    data=pd.read_csv(io.StringIO(txt),sep=sep,engine="python",on_bad_lines="skip")
                    out["table"]=profile(data)
                    out["decode"]=encoding
                    break
                except Exception as e:out["parse_error"]=type(e).__name__+":"+str(e)[:200]
    if name.endswith((".rmd",".md",".r",".readme")):
        text=b.decode("utf-8-sig","replace")
        out["source_comment_head"]=text[:1200]
    return out

def profile(df):
    df.columns=[str(c).strip() for c in df.columns]
    o={"rows":len(df),"columns":list(df.columns),
       "missing_by_column":{str(k):int(v) for k,v in df.isna().sum().items() if v},
       "categorical_levels":{},"numeric_profiles":{}}
    for c in df.columns:
        x=df[c]
        if str(c).lower() in ("latitude","longitude","lat","lon","site","location","reference","doi","url"):
            continue
        uniq=x.astype("string").value_counts(dropna=False)
        if 1<len(uniq)<=10:
            o["categorical_levels"][str(c)]={str(k):int(v) for k,v in uniq.items()}
        nu=pd.to_numeric(x,errors="coerce")
        if len(df)>0 and nu.notna().sum()>max(4,.80*len(df)):
            o["numeric_profiles"][str(c)]={"n":int(nu.notna().sum()),
                "min":float(nu.min()),"median":float(nu.median()),
                "max":float(nu.max())}
    return o

def inspect_blob(name,b,depth=0):
    result=report_table(name,b)
    if zipfile.is_zipfile(io.BytesIO(b)) and depth<2:
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            members=[]
            for entry in z.infolist():
                if entry.is_dir():continue
                if entry.file_size>8_000_000:
                    members.append({"name":Path(entry.filename).name,"error":"file_size_cap"})
                    continue
                members.append(inspect_blob(Path(entry.filename).name,z.read(entry),depth+1))
            result["zip_entries"]=members
    return result

def figshare():
    base=f"https://api.figshare.com/v2/collections/{FIGSHARE_COLLECTION}"
    result={"doi":"10.6084/m9.figshare.c.6837249","public_collection_id":FIGSHARE_COLLECTION}
    metadata=pull(base,2_000_000)
    result["meta_http_status"]=metadata.get("status")
    result["meta_error"]=metadata.get("error")
    j=js(metadata)
    if isinstance(j,dict):
        result["title"]=j.get("title")
        result["resource_doi"]=j.get("doi")
        result["public_date"]=j.get("published_date")
        result["meta_keys"]=list(j)[:35]
    listres=pull(base+"/articles?page_size=100",2_000_000)
    result["articles_list_http_status"]=listres.get("status")
    articles=js(listres)
    if isinstance(articles,dict):
        articles=articles.get("items",articles.get("articles",[]))
    if not isinstance(articles,list):
        result["articles_error"]=listres.get("error") or "non_array"
        articles=[]
    result["n_article_records"]=len(articles)
    result["articles"]=[]
    for article in articles[:30]:
        if not isinstance(article,dict):continue
        item={"article_id":article.get("id"),"listed_title":article.get("title")}
        aid=article.get("id")
        if not isinstance(aid,int):result["articles"].append(item);continue
        detail=js(pull(f"https://api.figshare.com/v2/articles/{aid}",2_000_000))
        if not isinstance(detail,dict):result["articles"].append(item);continue
        item["doi"]=detail.get("doi")
        item["license"]=detail.get("license")
        item["n_files"]=len(detail.get("files",[]))
        item["files"]=[]
        for file in detail.get("files",[])[:45]:
            name=file.get("name") or ""
            size=file.get("size")
            f={"name":name,"id":file.get("id"),"size":size,
               "md5":file.get("supplied_md5") or file.get("computed_md5")}
            if isinstance(size,int) and size<=CAP and file.get("download_url"):
                b=pull(file["download_url"],CAP)
                f["download_http_status"]=b.get("status")
                f["download_error"]=b.get("error")
                if b.get("status")==200:
                    f["content"]=inspect_blob(name,b["body"])
            item["files"].append(f)
        result["articles"].append(item)
    return result

def dryad():
    base="https://datadryad.org/api/v2/datasets/"+urllib.parse.quote("doi:"+DRYAD_DOI,safe="")
    result={"doi":DRYAD_DOI}
    meta=js(pull(base,1_000_000))
    if meta:
        result["title"]=meta.get("title")
        result["license"]=meta.get("license")
        result["visibility"]=meta.get("visibility")
    versions=js(pull(base+"/versions",1_000_000))
    if not isinstance(versions,dict):
        result["error"]="cannot_read_official_versions"
        return result
    versions=versions.get("_embedded",{}).get("stash:versions",[])
    ids=[]
    for v in versions:
        link=v.get("_links",{}).get("self",{}).get("href","")
        m=re.search(r"/versions/(\d+)",link)
        if m:ids.append((int(v.get("versionNumber",0)),int(m.group(1))))
    if not ids:
        result["error"]="no_public_version"
        return result
    version=max(ids)[1]
    result["version_id"]=version
    source=pull(f"https://datadryad.org/api/v2/versions/{version}/download",CAP)
    result["archive_status"]=source.get("status")
    if source.get("status")!=200:
        result["error"]=source.get("error")
        return result
    if not zipfile.is_zipfile(io.BytesIO(source["body"])):
        result["error"]="not_zip"
        return result
    with zipfile.ZipFile(io.BytesIO(source["body"])) as z:
        result["files"]=[]
        for info in z.infolist():
            if info.is_dir():continue
            result["files"].append(inspect_blob(Path(info.filename).name,z.read(info)))
    return result

def main():
    o={"status":"sources_audited",
       "figshare_228_pollen_fates":figshare(),
       "dryad_2023_pollen_ovules_pte":dryad(),
       "not_claiming_228_all_have_both_measures":True,
       "original_species_location_rows_not_released":True}
    (OUT/"johnson_harder2023_figshare_dryad_source_inventory_v14.json").write_text(
        json.dumps(o,indent=2,ensure_ascii=False)+"\n")
    print("GLOBAL_POLLEN_SOURCE_INVENTORY",json.dumps(o,ensure_ascii=False)[:70000],flush=True)

if __name__=="__main__":main()
