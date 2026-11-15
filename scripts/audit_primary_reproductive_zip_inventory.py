#!/usr/bin/env python3
"""Read published Dryad version ZIPs via the OFFICIAL anonymous REST endpoint.

This is the officially documented /api/v2/versions/<id>/download endpoint.
It is not a workaround for file-level 401: the public published-version ZIP
explicitly responds 200 without credentials. No proxies, tokens, mirrors or
cookies are used. We parse only named public research datasets and export
aggregate schema/provenance, not plant or site IDs.
"""
from __future__ import annotations
import csv
from collections import Counter
import hashlib
import io
import json
import urllib.request
import zipfile
from pathlib import Path

OUT=Path("analysis_results_primary_repro_csv_inventory")
OUT.mkdir(exist_ok=True)
HOST="https://datadryad.org"
VERSIONS={
  "Silene_latifolia_2011":{"version_id":11856,"doi":"10.5061/dryad.20k62r83"},
  "Brassica_rapa_2026":{"version_id":434166,"doi":"10.5061/dryad.hx3ffbgtq"},
  "Silene_dioica_2022":{"version_id":251371,"doi":"10.5061/dryad.0rxwdbs4d"},
  "Veratrum_grandiflorum_2025":{"version_id":405186,"doi":"10.5061/dryad.9s4mw6mw2"},
  "Clarkia_2024":{"version_id":332601,"doi":"10.5061/dryad.k3j9kd56r"}
}
MAX_COMPRESSED=800_000
MAX_TOTAL_UNCOMPRESSED=7_000_000
MAX_FILE=1_500_000
GENET_KEYS={"Individual","Individual.hand","Plant","plant","plant_id","id","Id","ID","plantID","ind"}
CATEGORY_KEYS={"Treatment","treatment","Pollination","Level","Sex","sex","line","Line","focal_sp","Raceme","Frt","fertilized","develop","abort","herbivory"}
NO_ROWS_EXPORT=True

def retrieve(version_id):
    url=f"{HOST}/api/v2/versions/{version_id}/download"
    request=urllib.request.Request(url,headers={"User-Agent":"Published-reproductive-data-verification/1.0","Accept":"application/zip"})
    with urllib.request.urlopen(request,timeout=70) as r:
        blob=r.read(MAX_COMPRESSED+1)
        if len(blob)>MAX_COMPRESSED:raise RuntimeError("public zip exceeds size cap")
        if r.status!=200:raise RuntimeError("unexpected public ZIP HTTP "+str(r.status))
        return blob, r.geturl()

def open_csv(data:bytes):
    txt=data.decode("utf-8-sig","replace")
    try:dialect=csv.Sniffer().sniff(txt[:8192],delimiters=",;\t")
    except csv.Error:dialect=csv.excel
    dr=csv.DictReader(io.StringIO(txt),dialect=dialect)
    return list(dr.fieldnames or []),list(dr)

def provenance():
    out={"status":"failed","access":"official published-version ZIP returned 200 in earlier source audit","restricted_to_public_data":True,"published_versions":{},"no_individual_id_or_geography_exported":True}
    for name,spec in VERSIONS.items():
        blob,redirect=retrieve(spec["version_id"])
        if not zipfile.is_zipfile(io.BytesIO(blob)):raise ValueError(f"{name} content not ZIP")
        info={"doi":spec["doi"],"version_id":spec["version_id"],"source_url":f"{HOST}/api/v2/versions/{spec['version_id']}/download",
              "archive_bytes":len(blob),"archive_sha256":hashlib.sha256(blob).hexdigest(),"files":[]}
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            names=z.namelist()
            total=sum(x.file_size for x in z.infolist())
            if total>MAX_TOTAL_UNCOMPRESSED:raise ValueError("uncompressed file size cap")
            if any(x.file_size>MAX_FILE for x in z.infolist()):raise ValueError("per-file size cap")
            for member in z.infolist():
                path=member.filename
                if member.is_dir():continue
                if path.startswith("/") or ".." in path.split("/"):raise ValueError("unexpected ZIP entry path")
                b=z.read(member)
                f={"name":path.split("/")[-1],"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}
                if path.lower().endswith((".csv",".txt")):
                    headers,rows=open_csv(b)
                    f["rows"]=len(rows)
                    f["headers"]=headers
                    f["categorical_audit"]={}
                    f["id_column_unique_counts"]={}
                    for col in headers:
                        if col in GENET_KEYS:
                            f["id_column_unique_counts"][col]=len({r.get(col) for r in rows if r.get(col)})
                        if col in CATEGORY_KEYS:
                            freq=Counter(r.get(col,"").strip() for r in rows)
                            if len(freq)<=14:
                                f["categorical_audit"][col]={k:int(v) for k,v in freq.items()}
                    f["missing_by_column"]={c:sum(r.get(c) in ("",None,"NA","NaN","nan") for r in rows) for c in headers if c in CATEGORY_KEYS or c in GENET_KEYS}
                info["files"].append(f)
        out["published_versions"][name]=info
        print("RAW_INVENTORY",name,json.dumps(info,ensure_ascii=False)[:13000],flush=True)
    out["status"]="success_original_published_raw_zip_parsed"
    (OUT/"published_zip_inventory_20261115.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("SOURCE_ARCHIVES_PARSED",len(out["published_versions"]),flush=True)

if __name__=="__main__":
    provenance()
