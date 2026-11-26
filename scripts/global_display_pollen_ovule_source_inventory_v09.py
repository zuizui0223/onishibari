#!/usr/bin/env python3
"""Discover authors' published global flower display, pollen, ovule data.

Calviño & Aizen 2023 AJB doi:10.1002/ajb2.16180, original published data
doi:10.5061/dryad.dz08kps20. Official anonymous Dryad v2 dataset and
published-version ZIP endpoints only. No scraping, mirror, credentials or
material re-distribution. Report SHA256 of uncompressed originals, CSV schema,
source-derived species grain and missingness; do not export full source rows.
"""
from __future__ import annotations
import csv, hashlib, io, json, re, urllib.request, urllib.parse, zipfile
from collections import Counter
from pathlib import Path
OUT=Path("analysis_results_global_flower_pollen_inventory_v09");OUT.mkdir(exist_ok=True)
BASE="https://datadryad.org"
DOI="10.5061/dryad.dz08kps20"
MAX_ZIP=12_000_000
MAX_UNPACK=55_000_000
SPECIES_KEYS=("species","species_name","taxon","taxa","binomial","sp")

def get(url,limit):
    request=urllib.request.Request(url,headers={"User-Agent":"Public-original-plant-trait-analysis/0.9","Accept":"application/json, application/zip, text/csv"})
    with urllib.request.urlopen(request,timeout=85) as r:
        data=r.read(limit+1)
        if len(data)>limit:raise ValueError("file larger than allowed cap")
        return data,int(r.status)

def parse_csv(content):
    txt=content.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
    sample=txt[:6000]
    try:dialect=csv.Sniffer().sniff(sample,delimiters=",;\t")
    except csv.Error:dialect=csv.excel
    reader=csv.DictReader(io.StringIO(txt),dialect=dialect)
    return reader.fieldnames or [],list(reader)

def provenance_info():
    dataset_url=BASE+"/api/v2/datasets/"+urllib.parse.quote("doi:"+DOI,safe="")
    d,_=get(dataset_url,1_000_000)
    meta=json.loads(d)
    versions,_=get(dataset_url+"/versions",1_000_000)
    vdoc=json.loads(versions)
    v=vdoc.get("_embedded",{}).get("stash:versions",[])
    if not v:
        limited={"status":"no_published_version_listed_by_official_api",
          "doi":DOI,"dataset_meta_keys":list(meta)[:50],
          "dataset_visibility":meta.get("visibility"),
          "dataset_version_number":meta.get("versionNumber"),
          "dataset_link_rels":list(meta.get("_links",{}).keys()),
          "dataset_stash_version_link":meta.get("_links",{}).get("stash:version"),
          "versions_total":vdoc.get("total"),"versions_count":vdoc.get("count"),
          "versions_embedded_keys":list(vdoc.get("_embedded",{}).keys()),
          "not_new_empirical_source_rows":True}
        (OUT/"global_flower_pollen_ovule_source_inventory_v09.json").write_text(json.dumps(limited,indent=2)+"\n")
        print("SOURCE_UNAVAILABLE_PUBLIC_VERSION",json.dumps(limited),flush=True)
        return
    eligible=[]
    for item in v:
        href=item.get("_links",{}).get("self",{}).get("href","")
        m=re.search(r"/versions/(\d+)",href)
        if m:eligible.append((int(item.get("versionNumber",0)),int(m.group(1))))
    if not eligible:raise ValueError("no published version ID")
    version=max(eligible)[1]
    archive_url=f"{BASE}/api/v2/versions/{version}/download"
    zipdata,_=get(archive_url,MAX_ZIP)
    if not zipfile.is_zipfile(io.BytesIO(zipdata)):raise ValueError("not published ZIP")
    report={
      "status":"failed",
      "doi":DOI,"primary_publication":"doi:10.1002/ajb2.16180",
      "source":archive_url,"version_id":version,
      "dataset_visibility":meta.get("visibility"),"license":meta.get("license"),
      "versionNumber":meta.get("versionNumber"),"archive_bytes":len(zipdata),
      "zip_sha256_informational_not_stable":hashlib.sha256(zipdata).hexdigest(),
      "original_source_tables":[],"full_rows_not_exported":True}
    with zipfile.ZipFile(io.BytesIO(zipdata)) as z:
        if sum(i.file_size for i in z.infolist())>MAX_UNPACK:
            raise ValueError("uncompressed cap exceeded")
        for x in z.infolist():
            fn=Path(x.filename).name
            if not fn or x.is_dir():continue
            b=z.read(x)
            f={"name":fn,"bytes":len(b),"sha256_uncompressed":hashlib.sha256(b).hexdigest()}
            if fn.lower().endswith((".csv",".txt",".tsv")):
                cols,rows=parse_csv(b)
                f.update({"columns":cols,"rows":len(rows)})
                species_cols=[c for c in cols if any(k in str(c).lower() for k in SPECIES_KEYS)]
                f["candidate_species_columns"]=species_cols
                f["candidate_species_unique_n"]={}
                for c in species_cols[:8]:
                    ids={str(r.get(c,"")).strip() for r in rows if str(r.get(c,"")).strip()}
                    f["candidate_species_unique_n"][c]=len(ids)
                f["numeric_like_columns"]=[]
                f["numeric_missing_audit"]={}
                for col in cols:
                    nonblank=[str(r.get(col,"")).strip() for r in rows if str(r.get(col,"")).strip().lower() not in ("","na","nan","null","none")]
                    if not nonblank:continue
                    try:
                        vals=[float(v) for v in nonblank]
                        ratio=len(vals)/len(nonblank)
                    except ValueError:
                        parsed=0
                        for vv in nonblank:
                            try:float(vv);parsed+=1
                            except ValueError:pass
                        ratio=parsed/len(nonblank)
                    if ratio>=.9:
                        f["numeric_like_columns"].append(col)
                        f["numeric_missing_audit"][col]={"n_nonmissing_candidate":len(nonblank),"n_missing":len(rows)-len(nonblank)}
            elif fn.lower().endswith((".r",".rmd",".md",".readme")):
                # Short source descriptions only; the original content remains on Dryad.
                f["source_text_head"]=b[:1700].decode("utf-8","replace")
            report["original_source_tables"].append(f)
            print("SOURCE_FILE",json.dumps(f,ensure_ascii=False)[:12000],flush=True)
    report["status"]="success_original_public_archive_structural_audit"
    (OUT/"global_flower_pollen_ovule_source_inventory_v09.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("GLOBAL_SOURCE_AUDIT_SUCCESS",len(report["original_source_tables"]),flush=True)
if __name__=="__main__":provenance_info()
