#!/usr/bin/env python3
"""Audit Bonnaudet et al. 2026 openly deposited Angiosperms353 phylogenies.

Discovery & license constraints:
- Zenodo https://doi.org/10.5281/zenodo.17517526
- Inspect actual taxon-level names; do not equate family phylogeny with full coverage
  of the P/C/K East Asian focal trio.
- Download only small tree/readme files, never bulk supercontig archives.
"""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import re
import time
from urllib.request import Request,urlopen
from Bio import Phylo

RECORD="17517526"
API=f"https://zenodo.org/api/records/{RECORD}"
OUT=Path("analysis_results_zenodo_thymelaeaceae")
OUT.mkdir(exist_ok=True)
HEADERS={"User-Agent":"Daphne-onishibari-source-inventory/0.3 academic research"}
MAX_DOWNLOAD=25*1024*1024

def get(url,limit=MAX_DOWNLOAD):
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers=HEADERS),timeout=90) as r:
                b=r.read(limit+1)
                if len(b)>limit: raise ValueError(f"oversized response: {url}")
                return b
        except Exception:
            if attempt==3: raise
            time.sleep(2**attempt)

def main():
    report={"source_doi":"10.5281/zenodo.17517526",
       "citation":"Bonnaudet et al. 2026 Annals of Botany doi:10.1093/aob/mcag206",
       "api":API,"has_raw_353_loci":False,"focal_taxa_present":{},
       "records":[],"trees":[]}
    try:
        raw=get(API,limit=5_000_000)
        meta=json.loads(raw)
        (OUT/"zenodo_record.json").write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n")
        report["record_id"]=meta.get("id")
        report["record_updated"]=meta.get("updated")
        report["record_status"]=meta.get("status")
        report["license"]=((meta.get("metadata") or {}).get("license"))
        files=meta.get("files",[])
        if not files:raise ValueError("Zenodo record has no file listing")
        for f in files:
            name=f.get("key") or f.get("filename") or ""
            size=f.get("size")
            links=f.get("links") or {}
            url=links.get("self") or links.get("download") or links.get("content")
            file_row={"key":name,"size":size,"url":url,
               "checksum":f.get("checksum"),"downloaded":False}
            report["records"].append(file_row)
            if re.search(r"(loci|supercontig|aligned)",name,re.I):
                report["has_raw_353_loci"]=True
            if not url or not size or size>MAX_DOWNLOAD:continue
            if not (name.lower().endswith((".tre",".tree",".nwk",".txt",".nex"))
                    or "readme" in name.lower()):continue
            rawfile=get(url,limit=MAX_DOWNLOAD)
            if size is not None and len(rawfile)!=size:
                raise ValueError(f"{name}: remote bytes size mismatch {len(rawfile)} vs {size}")
            path=OUT/Path(name).name
            path.write_bytes(rawfile)
            file_row.update(downloaded=True,sha256=hashlib.sha256(rawfile).hexdigest())
            if not name.lower().endswith((".tre",".tree",".nwk")):continue
            try:
                tree=Phylo.read(io.StringIO(rawfile.decode("utf8")),"newick")
                tips=[t.name or "" for t in tree.get_terminals()]
                if len(tips)<3: raise ValueError("too few tree tips")
            except Exception as exc:
                file_row["tree_parse_error"]=str(exc)
                continue
            (OUT/(path.stem+"_tips.txt")).write_text("\n".join(tips)+"\n")
            daphne=[t for t in tips if re.search(r"(^|[_\s-])daphne([_\s-]|$)",t,re.I)]
            # also search for derived filenames or genus prefix within composite labels
            if not daphne:
                daphne=[t for t in tips if "daphne" in t.lower()]
            focal_terms=[
                "pseudomezereum","kamtschatica","koreana","jezoensis",
                "miyabeana","kiusiana","laureola","reichsteinii","rebunensis"
            ]
            hits={term:[t for t in tips if term in t.lower()] for term in focal_terms}
            x={"tree_file":name,"n_taxa":len(tips),
               "n_daphne_named_tips":len(daphne),
               "daphne_named_tips":daphne,
               "focal_name_matches":hits,
               "first_five_tips":tips[:5]}
            report["trees"].append(x)
        report["completed"]=True
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        report["completed"]=False
    (OUT/"inventory.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"completed":report["completed"],"error":report.get("error"),
      "record_count":len(report["records"]), "records":report["records"],
      "trees":report["trees"]},ensure_ascii=False,indent=2))
    if not report["completed"]:raise RuntimeError(report.get("error"))
if __name__=="__main__":main()
