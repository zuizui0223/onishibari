#!/usr/bin/env python3
"""Audit direct GenBank BioProject/BioSample/SRA links of the 2025 focal plastomes.

Lack of an SRA link in a GenBank record does not prove raw reads are unpublished.
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
ACC = {"P":"PV752465", "K":"PV752464", "C":"ON244034"}
out=Path("analysis_results_metadata")
out.mkdir(exist_ok=True)
report={"source":"NCBI GenBank flat-file links; not a comprehensive SRA index search",
        "date":"2026-11-15","records":{}}
for label,acc in ACC.items():
    params=urllib.parse.urlencode(dict(db="nuccore",id=acc,rettype="gb",retmode="text"))
    request=urllib.request.Request(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+params,
        headers={"User-Agent":"Daphne-systematics-audit/0.2 (scientific)"})
    with urllib.request.urlopen(request,timeout=90) as r:
        s=r.read().decode("utf8")
    if not s.startswith("LOCUS") or acc not in s[:3000]:
        raise ValueError("Invalid GenBank flatfile response: "+acc)
    # Leave full record available as source artifact, but display only identifiers.
    (out/(acc+".gb")).write_text(s,encoding="utf8")
    links=re.findall(r"(BioProject|BioSample|SRA):\s*([A-Za-z0-9_]+)",s,re.IGNORECASE)
    report["records"][label]={"accession":acc,"direct_links":links,
      "db_link_section":re.findall(r"(?m)^DBLINK\s+.*(?:\n {12}.*)*",s)[:3],
      "record_length_chars":len(s),
      "no_direct_SRA_link_in_flatfile":not any(k.lower()=="sra" for k,v in links)}
    time.sleep(0.4)
(out/"audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
print(json.dumps(report,indent=2))
