#!/usr/bin/env python3
"""NCBI SRA index audit for focal Daphne taxa, synonyms, and phylogenomic comparators.

Queries indexed SRA database at NCBI, not just GenBank nucleotide flatfile DBLINK.
No results in this index query is not proof that unpublished reads do not exist.
"""
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERIES = {
    "P_original": '"Daphne pseudomezereum"[All Fields]',
    "P_var_nominate": '"Daphne pseudomezereum var. pseudomezereum"[All Fields]',
    "K_original": '"Daphne kamtschatica"[All Fields]',
    "C_historical": '"Daphne pseudomezereum var. koreana"[All Fields]',
    "C_former_species": '"Daphne koreana"[All Fields]',
    "jezoensis": '"Daphne jezoensis"[All Fields]',
    "rebunensis": '"Daphne rebunensis"[All Fields]',
    "kiusiana_reference": '"Daphne kiusiana"[All Fields]',
    "reichsteinii_reference": '"Daphne reichsteinii"[All Fields]',
    "new_hybseq_project": '"PRJNA1493253"[All Fields]',
}
HEAD = {"User-Agent": "Daphne-comparative-SRA-audit/0.3 (academic)"}

def get(op, params):
    u=ROOT+op+"?"+urllib.parse.urlencode(params)
    for i in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(u,headers=HEAD),timeout=90) as h:
                return h.read()
        except (TimeoutError,urllib.error.URLError,urllib.error.HTTPError):
            if i==3: raise
            time.sleep(2**i)
    raise RuntimeError("request retries exhausted")

def main():
    out=Path("analysis_results_sra")
    out.mkdir(exist_ok=True)
    report={
      "question":"Do indexed NCBI SRA study/run records exist for focal P/K/C or useful comparators?",
      "date_utc":"2026-11-15",
      "scope":"NCBI SRA esearch exact quoted phrases, All Fields; retmax 35 metadata each",
      "queries":{},
      "caveat":"A null result does NOT imply reads unsubmitted anywhere; paper says 2025 Illumina P/K reads generated."
    }
    for name,term in QUERIES.items():
        xml=get("esearch.fcgi",dict(db="sra",term=term,retmode="xml",retmax="35"))
        root=ET.fromstring(xml)
        count=int(root.findtext("Count","0"))
        ids=[e.text for e in root.findall("./IdList/Id") if e.text]
        result={"term":term,"count":count,"ids":ids,"example_records":[]}
        if ids:
            records=get("esummary.fcgi",dict(db="sra",id=",".join(ids[:15]),retmode="xml"))
            (out/(name+"_esummary.xml")).write_bytes(records)
            rs=ET.fromstring(records)
            for row in list(rs):
                uid=row.attrib.get("uid",row.attrib.get("Id"))
                chunks=[(x.attrib.get("Name",""), (x.text or "")[:1500]) for x in list(row)]
                result["example_records"].append({"uid":uid, "fields":chunks[:8]})
        report["queries"][name]=result
        print(name,"SRA_count",count,"IDs",ids[:5],flush=True)
        time.sleep(0.45)
    (out/"index_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")

if __name__=="__main__":
    main()
