#!/usr/bin/env python3
"""Source-checked ENA/NCBI archive coverage audit for Daphne P-C-K.

Does NOT claim absence of unpublished WGS or nuclear data. Distinguishes
NCBI project raw-read quantities from originally sequenced raw volumes.
"""
import csv
import io
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from xml.etree import ElementTree

OUT=Path("analysis_results_crossarchive")
OUT.mkdir(exist_ok=True)
UA={"User-Agent":"Daphne-cytonuclear-archive-audit/0.5 (research)"}
PROJECT="PRJNA835716"
RUN="SRR19117239"
TERMS={
 "P_voucher": '"Lee-Dpp230303"[All Fields]',
 "K_voucher": '"K005971094"[All Fields]',
 "C_2025_voucher": '"Byeon-Dpk250701"[All Fields]',
 "P_plastome": '"PV752465"[All Fields]',
 "K_plastome": '"PV752464"[All Fields]',
 "C_plastome": '"ON244034"[All Fields]',
 "C_2023_project": '"PRJNA835716"[All Fields]',
}
FIELDS="run_accession,study_accession,sample_accession,experiment_accession,scientific_name,instrument_platform,library_strategy,library_source,library_selection,read_count,base_count,fastq_bytes,fastq_ftp,submitted_ftp"
def get(url,maxbytes=5_000_000):
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers=UA),timeout=75) as r:
                buf=r.read(maxbytes+1)
                if len(buf)>maxbytes:raise ValueError("Response exceeds safety ceiling")
                return buf.decode("utf8")
        except (OSError,TimeoutError):
            if attempt==3:raise
            time.sleep(2**attempt)
def ena_report(accession):
    u="https://www.ebi.ac.uk/ena/portal/api/filereport?"+urlencode({
        "accession":accession,"result":"read_run","fields":FIELDS,"format":"tsv"})
    txt=get(u)
    if txt.lstrip().startswith(("{","<html","ERROR")):raise ValueError("ENA non-TSV: "+txt[:250])
    rows=list(csv.DictReader(io.StringIO(txt),delimiter="\t"))
    if rows and "run_accession" not in rows[0]:raise ValueError("ENA header mismatch "+txt[:300])
    return {"url":u,"n_runs":len(rows),"runs":rows}
def ncbi_ids(term):
    u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"+urlencode({
        "db":"sra","term":term,"retmode":"xml","retmax":25})
    xml=ElementTree.fromstring(get(u).encode())
    count=int(xml.findtext("Count","0"))
    ids=[n.text for n in xml.findall("./IdList/Id")]
    return {"term":term,"sra_record_count":count,"ids":ids}
def main():
    summary={"source":"ENA portal API filereport for exact BioProject/BioSample/run; NCBI ESearch SRA for vouchers",
      "date":"2026-11-15","record":{},"ncbi_sra_queries":{},
      "expected_source_2023_raw_bases":7340000000,
      "notes":["Raw paper 7.34Gbp is reported pre-filtering, compare to public SRA bases cautiously.",
       "An empty archive index or voucher query does not establish no archived private data.",
       "SRR19117239 public data was previously observed as >99.9% plastid-mapped and 0 ITS-mapped."]}
    try:
        for key in [PROJECT,"SAMN28108557",RUN]:
            summary["record"][key]=ena_report(key)
            time.sleep(0.4)
        runs=summary["record"][PROJECT]["runs"]
        known=[r for r in runs if r.get("run_accession")==RUN]
        if len(known)!=1:raise ValueError("Expected public Korean run absent from ENA PRJNA835716")
        if summary["record"][PROJECT]["n_runs"]!=summary["record"][RUN]["n_runs"]:
            raise ValueError("ENA project vs run report is not trivially comparable")
        reported=int(known[0].get("base_count") or 0)
        summary["visible_bases"]=reported
        summary["visible_fraction_of_paper_7_34Gbp"]=round(reported/7340000000,6)
        for label,term in TERMS.items():
            summary["ncbi_sra_queries"][label]=ncbi_ids(term)
            time.sleep(0.4)
        summary["complete"]=True
    except Exception as exc:
        summary["complete"]=False
        summary["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"archive_crosswalk.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    if not summary["complete"]:raise RuntimeError(summary["error"])
if __name__=="__main__":main()
