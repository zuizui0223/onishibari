#!/usr/bin/env python3
"""Check anonymous officially documented Dryad REST metadata/file access.

No credential, token, user-agent spoofing, proxy, mirror, bypass, or automated
login. Published data license may allow reuse; metadata alone is not a dataset.

Original docs: https://github.com/datadryad/dryad-app/tree/main/documentation/apis
Run only on known small explicitly public datasets, cap downloads at 700 KB.
When no anonymous file access exists, report 'inaccessible' and stop.
"""
from __future__ import annotations
import json, re, urllib.request, urllib.error, urllib.parse
from pathlib import Path

OUT=Path("analysis_results_dryad_access_20261115")
OUT.mkdir(exist_ok=True)
BASE="https://datadryad.org/api/v2"
DOIS={
    "Silene_latifolia_2011":"10.5061/dryad.20k62r83",
    "Brassica_rapa_2026":"10.5061/dryad.hx3ffbgtq",
    "Silene_dioica_2022":"10.5061/dryad.0rxwdbs4d",
    "Veratrum_grandiflorum_2025":"10.5061/dryad.9s4mw6mw2",
    "Clarkia_reproduction_2024":"10.5061/dryad.k3j9kd56r"
}
MAX_BYTES=700_000

def fetch_json(url, max_bytes=500_000):
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"Scientific-reproducibility-audit/0.1"})
    try:
        with urllib.request.urlopen(req,timeout=18) as r:
            data=r.read(max_bytes+1)
            if len(data)>max_bytes:return {"status":r.status,"error":"response too large"}
            try:obj=json.loads(data)
            except (ValueError,UnicodeDecodeError):return {"status":r.status,"error":"unexpected non-json content","content_type":r.headers.get("Content-Type")}
            return {"status":r.status,"json":obj}
    except urllib.error.HTTPError as e:
        return {"status":e.code,"error": "http_error"}
    except Exception as e:
        return {"error":type(e).__name__}

def fetch_small_file(url):
    req=urllib.request.Request(url,headers={"Accept":"text/csv, text/plain, application/octet-stream","User-Agent":"Scientific-reproducibility-audit/0.1"})
    try:
        with urllib.request.urlopen(req,timeout=20) as r:
            b=r.read(MAX_BYTES+1)
            ctype=r.headers.get("Content-Type","")
            if len(b)>MAX_BYTES:return {"status":r.status,"error":"more_than_file_cap","content_type":ctype}
            looks_binary=(b.count(b"\x00")>0)
            return {"status":r.status,"bytes":len(b),"content_type":ctype,
                "looks_csv_or_text":not looks_binary,
                "preview_header":b[:180].decode("utf-8","replace").replace("\r"," ").replace("\n"," ") if not looks_binary else None}
    except urllib.error.HTTPError as e:return {"status":e.code,"error":"http_error"}
    except Exception as e:return {"error":type(e).__name__}

def find_version_id(obj):
    if isinstance(obj,dict):
        for key in ("id","versionId","versionID"):
            v=obj.get(key)
            if isinstance(v,int) or (isinstance(v,str) and v.isdigit()):
                return int(v)
        # Dryad API versions commonly use canonical links rather than an ID field.
        links=obj.get("_links",{})
        for link_key in ("self","stash:version","stash:files"):
            v=links.get(link_key,{}) if isinstance(links,dict) else {}
            href=v.get("href","") if isinstance(v,dict) else ""
            match=re.search(r"/versions/(\d+)",str(href))
            if match:return int(match.group(1))
    return None

def main():
    report={"status":"probe_completed","method":"anonymous_read_only_official_Dryad_v2_REST","dataset_probes":{}}
    for key,doi in DOIS.items():
        doc={"doi":doi}
        dataset_url=BASE+"/datasets/"+urllib.parse.quote("doi:"+doi,safe="")
        metadata=fetch_json(dataset_url)
        doc["dataset_status"]=metadata.get("status")
        doc["dataset_error"]=metadata.get("error")
        if metadata.get("status")==200:
            d=metadata["json"]
            doc["dataset_meta_keys"]=list(d)[:35]
            doc["canonical_dataset_id"]=d.get("id")
            doc["versionNumber"]=d.get("versionNumber")
        versions=fetch_json(dataset_url+"/versions")
        doc["versions_status"]=versions.get("status")
        doc["versions_error"]=versions.get("error")
        vid=None
        if versions.get("status")==200:
            j=versions["json"]
            doc["version_meta_keys"]=list(j)[:35] if isinstance(j,dict) else ["array"]
            cand=j.get("_embedded",{}).get("stash:versions",[]) if isinstance(j,dict) else j
            if isinstance(cand,list) and cand:
                doc["n_versions_returned"]=len(cand)
                doc["version_candidate_public_keys"]=[list(v)[:20] if isinstance(v,dict) else str(type(v)) for v in cand[:3]]
                doc["version_candidate_links"]=[v.get("_links",{}) if isinstance(v,dict) else None for v in cand[:3]]
                vid=find_version_id(cand[-1]) or find_version_id(cand[0])
            if not vid and isinstance(j,dict):vid=find_version_id(j)
        if not vid and metadata.get("status")==200:
            v=metadata["json"].get("versionNumber")
            # versionNumber is not the numeric version ID; do not misuse it.
        doc["chosen_version_id"]=vid
        if vid:
            files=fetch_json(BASE+"/versions/"+str(vid)+"/files")
            doc["files_status"]=files.get("status")
            doc["files_error"]=files.get("error")
            if files.get("status")==200:
                obj=files["json"]
                doc["files_keys"]=list(obj)[:20] if isinstance(obj,dict) else ["array"]
                candidates=obj.get("_embedded",{}).get("stash:files",[]) if isinstance(obj,dict) else obj
                doc["file_count"]=len(candidates) if isinstance(candidates,list) else 0
                doc["file_overview"]=[]
                if isinstance(candidates,list):
                    for f in candidates[:14]:
                        if not isinstance(f,dict):continue
                        fid=f.get("id")
                        name=f.get("path") or f.get("filename") or f.get("name") or ""
                        meta={"file_id":fid,"name":name,"size":f.get("size")}
                        if isinstance(fid,int) and isinstance(f.get("size"),int) and f["size"]<=MAX_BYTES and str(name).lower().endswith((".csv",".txt",".md")):
                            # Public endpoint documented by official Dryad API.
                            meta["public_download_probe"]=fetch_small_file(BASE+"/files/"+str(fid)+"/download")
                        doc["file_overview"].append(meta)
        report["dataset_probes"][key]=doc
        print("DATASET",key,json.dumps(doc,ensure_ascii=False)[:3500],flush=True)
    (OUT/"official_dryad_anonymous_access_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("PROBE_DONE",len(DOIS),"No dataset rows interpreted from metadata",flush=True)

if __name__=="__main__":
    main()
