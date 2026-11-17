#!/usr/bin/env python3
"""v0.42: attempt archival recovery of original 2010 Japanese Daphne barcode loci.

Original study: Kawahara et al., 121st Japanese Forest Society (2010),
DOI 10.11519/jfsc.121.0.332.0: examined rbcL, matK, trnH-psbA,
psaB-rpl4, rpl14-rps16, trnL-trnF (plastid) + ITS2 (multi-copy nrDNA).

Important: marker targets in an abstract are NOT GenBank depositions or an
independent multilocus species tree. NCBI's literature/year indexes and
author metadata are incomplete, so zero query hits DO NOT prove that
original sequences, other accession names, or privately held alignments
do not exist. Does not download or redistribute sequence data.
"""
from __future__ import annotations
import json,time,urllib.request,urllib.parse,xml.etree.ElementTree as ET
from pathlib import Path

OUT=Path("analysis_results_daphne_2010_barcode_archive_v42")
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SEARCHES={
 "Kawahara_author_year_2008_2013":
    'Daphne[Organism] AND Kawahara[Author] AND 2008:2013[PDAT]',
 "Kawahara_plus_team_author_year_2008_2015":
    'Daphne[Organism] AND (Kawahara[Author] OR Yoshimura[Author] OR Yamashita[Author]) AND 2008:2015[PDAT]',
 "focal_named_taxa_2010":
    '(Daphne pseudomezereum[Organism] OR Daphne jezoensis[Organism] OR Daphne kamtschatica[Organism]) AND 2008:2013[PDAT]',
 "Daphne_plastid_spacer_2010":
    'Daphne[Organism] AND (trnH-psbA[All Fields] OR psaB-rpl4[All Fields] OR rpl14-rps16[All Fields] OR trnL-trnF[All Fields]) AND 2008:2013[PDAT]',
 "Daphne_ITS2_2010":
    'Daphne[Organism] AND ITS2[All Fields] AND 2008:2013[PDAT]',
}
EXPECTED_PUBLISHED_MARKERS={
 "plastid":[ "rbcL", "matK", "trnH-psbA", "psaB-rpl4", "rpl14-rps16", "trnL-trnF" ],
 "multicopy_nuclear_ribosomal":["ITS2"]
}

def get(operation,params):
    url=BASE+operation+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={"User-Agent":"Daphne-2010-archival-sequence-audit/0.42 scientific"})
    last=None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req,timeout=55) as res:
                b=res.read(1_200_001)
            if len(b)>1_200_000:raise ValueError("NCBI response too large")
            root=ET.fromstring(b)
            return root,url
        except (OSError,ET.ParseError) as e:
            last=e
            if attempt==3:raise
            time.sleep(2**attempt)
    raise last

def search(term):
    node,url=get("esearch.fcgi",{"db":"nuccore","term":term,
                                 "retmode":"xml","retmax":30})
    assert node.tag=="eSearchResult"
    count=int(node.findtext("Count","0"))
    ids=[q.text for q in node.findall("./IdList/Id")]
    status=node.findtext("ERROR")
    if status:raise RuntimeError("NCBI search error: "+status[:200])
    return {"search_query":term,"count_indexed_records":count,
      "ids":ids,"query_url":url,
      "capped_at_30":count>30}

def describe(ids):
    if not ids:return []
    node,url=get("esummary.fcgi",{"db":"nuccore","id":",".join(ids),
                                   "retmode":"xml"})
    assert node.tag in ("eSummaryResult","eSummaryResultV2")
    records=[]
    for ds in node.findall(".//DocSum"):
        d={"nuccore_id":ds.findtext("Id","")}
        for item in ds.findall("Item"):
            name=item.attrib.get("Name","")
            if name in ("Caption","Title","TaxId","Length","CreateDate","UpdateDate"):
                d[name]=item.text
        records.append(d)
    if len(records)==0:raise RuntimeError("No parsable NCBI esummary DocSum metadata")
    return records

def main():
    OUT.mkdir(exist_ok=True)
    summary={"status":"NCBI_EUTILS_2010_DAPHNE_ARCHIVE_INDEX_AUDIT",
      "original_published_conference_doi":"10.11519/jfsc.121.0.332.0",
      "original_conference_url":"https://www.jstage.jst.go.jp/article/jfsc/121/0/121_0_332/_pdf/-char/ja",
      "historical_KAKEN":"20380093",
      "published_barcode_target_markers":EXPECTED_PUBLISHED_MARKERS,
      "n_independent_low_copy_nuclear_marker_targets_in_2010_abstract":0,
      "n_multicopy_nuclear_ribosomal_marker_targets":1,
      "n_plastid_marker_targets_not_independent_geneaologies":6,
      "queries":{},
      "bibliographic_limits":{
        "author_name_index_may_be_absent_on_sequence_submissions":True,
        "2010_presentation_date_not_sequence_deposit_date":True,
        "negative_indexed_query_does_not_prove_original_data_absent":True,
        "sequences_in_markers_not_yet_linked_to_original_voucher":True,
        "plastid_regions_share_one_organelle_genealogy":True,
        "ITS2_is_multicopy_not_independent_low_copy_gene":True,
        "not_a_multilocus_P_J_C_K_species_tree":True,
        "not_same_individual_P_C_K_sequence_evidence":True}}
    for key,term in SEARCHES.items():
        d=search(term)
        d["sample_summary"]=describe(d.pop("ids")[:20])
        summary["queries"][key]=d
        time.sleep(.55)
    # Index counts are complete for each query but displayed accession
    # summaries are deliberately capped to avoid a fake exhaustive catalog.
    if len(summary["queries"])!=len(SEARCHES):raise ValueError("Incomplete query set")
    outfile=OUT/"daphne_2010_barcode_genbank_recovery_index_v42.json"
    outfile.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    concise={
      "status":summary["status"],
      "search_counts":{k:v["count_indexed_records"] for k,v in summary["queries"].items()},
      "sample_accessions":{k:[x.get("Caption") for x in v["sample_summary"][:10]]
       for k,v in summary["queries"].items()},
      "seven_markers_still_only_one_plastid_and_one_multicopy_nuclear_genealogy":True,
      "no_conclusion_from_zero_hits":True}
    print("DAPHNE_2010_NCBI_ARCHIVE_V42",json.dumps(concise,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
