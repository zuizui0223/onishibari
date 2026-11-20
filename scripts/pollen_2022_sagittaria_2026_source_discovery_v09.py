#!/usr/bin/env python3
"""Discover two independent, public source datasets without inventing analyses.

A. Leme da Cunha et al 2022 GEB, Dryad doi:10.5061/dryad.bg79cnpc7
   pollen/ovule phenotypes and latitude (predates authors 2023 paper).
B. Tang & Dai 2026 Ecology and Evolution, Figshare article 27610539
   Sagittaria trifolia sequential inflorescences/light+pollinator access.

Read-only official public APIs. Export source file names, checksums, sizes and
metadata only, never private IDs. Record blocked/unavailable paths explicitly.
"""
from __future__ import annotations
import json,re,urllib.parse,urllib.request,zipfile,hashlib,io
from pathlib import Path
OUT=Path("analysis_results_pollen_allocation_multi_source_v09");OUT.mkdir(exist_ok=True)
BASE="https://datadryad.org"
DRYAD="10.5061/dryad.bg79cnpc7"
FIGSHARE_ARTICLE=27610539
def grab(url,maxsize=12_000_000):
    q=urllib.request.Request(url,headers={"User-Agent":"Plant-pollen-supply-fitness-source-discovery/0.9","Accept":"application/json, application/octet-stream"})
    try:
        with urllib.request.urlopen(q,timeout=60) as r:
            b=r.read(maxsize+1)
            if len(b)>maxsize:return {"status":r.status,"error":"size_cap","bytes":len(b)}
            return {"status":r.status,"content":b,"url_final":r.geturl(),"content_type":r.headers.get("Content-Type")}
    except Exception as e:return {"error":type(e).__name__+":"+str(e)[:180]}
def bodyjson(r):
    if r.get("status")!=200:return None
    try:return json.loads(r["content"])
    except (ValueError,TypeError):return None
def main():
    report={"status":"source_discovery_complete","pollen_2022":{},"sagittaria_2026":{},"no_row_level_interpretation":True}
    api=BASE+"/api/v2/datasets/"+urllib.parse.quote("doi:"+DRYAD,safe="")
    md=bodyjson(grab(api,700000))
    v=bodyjson(grab(api+"/versions",700000))
    pollen=report["pollen_2022"]
    pollen["doi"]=DRYAD
    if md:
        pollen["metadata_fields"]=list(md)[:25]
        pollen["title"]=md.get("title")
        pollen["versionNumber"]=md.get("versionNumber")
        pollen["visibility"]=md.get("visibility")
        pollen["license"]=md.get("license")
    versions=v.get("_embedded",{}).get("stash:versions",[]) if isinstance(v,dict) else []
    pollen["n_published_versions_listed"]=len(versions)
    listed=[]
    for x in versions:
        m=re.search(r"/versions/(\d+)",x.get("_links",{}).get("self",{}).get("href",""))
        if m:listed.append((x.get("versionNumber",0),int(m.group(1))))
    if listed:
        vid=max(listed)[1];pollen["version_id"]=vid
        z=grab(BASE+f"/api/v2/versions/{vid}/download",8_000_000)
        if z.get("status")==200 and zipfile.is_zipfile(io.BytesIO(z["content"])):
            with zipfile.ZipFile(io.BytesIO(z["content"])) as f:
                pollen["files"]=[]
                for i in f.infolist():
                    if i.is_dir():continue
                    b=f.read(i)
                    pollen["files"].append({"name":Path(i.filename).name,"bytes":len(b),
                                             "sha256_uncompressed":hashlib.sha256(b).hexdigest()})
        else:pollen["download_error"]=z.get("error") or z.get("status")
    fs=report["sagittaria_2026"]
    fs["article_id"]=FIGSHARE_ARTICLE
    fs["doi"]="10.6084/m9.figshare.27610539"
    info=bodyjson(grab(f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}",2_000_000))
    if info:
        fs["title"]=info.get("title")
        fs["published_date"]=info.get("published_date")
        fs["doi_api"]=info.get("doi")
        fs["license"]=info.get("license")
        fs["files"]=[{"name":f.get("name"),"id":f.get("id"),"size":f.get("size"),
            "md5":f.get("supplied_md5") or f.get("computed_md5"),
            "download_url":f.get("download_url")} for f in info.get("files",[])]
        for file in fs["files"][:10]:
            if file.get("size") is not None and int(file["size"])<=900000 and file.get("download_url"):
                r=grab(file["download_url"],1_000_000)
                file["probe"]={"status":r.get("status"),"error":r.get("error"),
                               "bytes":len(r.get("content",b""))}
                if r.get("status")==200:
                    file["sha256"]=hashlib.sha256(r["content"]).hexdigest()
                    if file["name"].lower().endswith((".csv",".tsv",".txt",".md")):
                        file["head"]=r["content"][:700].decode("utf-8","replace")
    else:
        fs["api_error"]="figshare_public_article_metadata_not_available"
    (OUT/"public_original_pollen_2022_sagittaria_2026_sources_v09.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("SOURCE_DISCOVERY",json.dumps(report,ensure_ascii=False)[:20000],flush=True)
if __name__=="__main__":main()
