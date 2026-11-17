#!/usr/bin/env python3
"""v0.43: discriminate original Kawahara et al. (2010) Daphne barcode data
from unrelated 2008–2013 GenBank ITS2 and plastid hits by full provenance.

The earlier source-author search returned 0 hits, broad ITS2 index returned
60 records, plastid spacer query 2. No hits are automatically assigned to
the Japanese original project without authors, deposited taxon and voucher.
Outputs metadata only, not raw nucleotide sequences or restricted vouchers.
"""
from __future__ import annotations
import json,re,time,urllib.request,urllib.parse
from pathlib import Path
from daphne_2010_barcode_genbank_recovery_v42 import search,SEARCHES

OUT=Path("analysis_results_daphne_2010_barcode_archive_provenance_v43")
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_flat(ids):
    assert 0<len(ids)<=15
    args={"db":"nuccore","id":",".join(ids),"rettype":"gb","retmode":"text"}
    req=urllib.request.Request(BASE+"?"+urllib.parse.urlencode(args),
      headers={"User-Agent":"Japanese-Daphne-barcode-provenance-v43 scientific"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req,timeout=80) as h:
                b=h.read(350_001)
            if len(b)>350_000:raise ValueError("GenBank record too large")
            s=b.decode("utf8")
            if "LOCUS" not in s or "ACCESSION" not in s:raise ValueError("Not a valid GenBank response")
            return s
        except OSError:
            if attempt>=3:raise
            time.sleep(2**attempt)

def extract(text):
    parts=re.split(r"(?m)^//\s*$",text)
    out=[]
    for source in parts:
        if "LOCUS" not in source:continue
        acc=re.search(r"(?m)^ACCESSION\s+([A-Z]{1,4}\d+)",source)
        definition=re.search(r"(?ms)^DEFINITION\s+(.+?)(?=^ACCESSION\s)",source)
        organism=re.search(r"(?m)^\s+ORGANISM\s+(.+)$",source)
        refblock=re.search(r"(?ms)^REFERENCE\s+\d+.+?(?=^REFERENCE\s+\d+|^COMMENT\s|^FEATURES\s|^ORIGIN\s)",source)
        authors=re.search(r"(?ms)^\s+AUTHORS\s+(.+?)(?=^\s+TITLE\s|^\s+JOURNAL\s)",refblock.group(0)) if refblock else None
        title=re.search(r"(?ms)^\s+TITLE\s+(.+?)(?=^\s+JOURNAL\s|^\s+PUBMED\s)",refblock.group(0)) if refblock else None
        journal=re.search(r"(?m)^\s+JOURNAL\s+(.+)$",refblock.group(0)) if refblock else None
        country=re.search(r'/country="([^"]+)"',source)
        coll=re.search(r'/collected_by="([^"]+)"',source)
        gene=re.search(r'/(?:gene|product)="([^"]+)"',source)
        date=re.search(r"(?m)^LOCUS\s+.+?\s+(\d{2}-[A-Z]{3}-\d{4})\s*$",source)
        out.append({
          "accession":acc.group(1) if acc else None,
          "organism":organism.group(1).strip() if organism else None,
          "short_definition":re.sub(r"\s+"," ",definition.group(1)).strip()[:180] if definition else None,
          "first_reference_authors":re.sub(r"\s+"," ",authors.group(1)).strip()[:250] if authors else None,
          "first_reference_title":re.sub(r"\s+"," ",title.group(1)).strip()[:230] if title else None,
          "first_reference_journal":journal.group(1).strip()[:150] if journal else None,
          "locus_last_update_date":date.group(1) if date else None,
          "reported_country":country.group(1) if country else None,
          "public_collected_by":coll.group(1) if coll else None,
          "marker_annotation_example":gene.group(1) if gene else None,
        })
    return out

def team_match(d):
    s=" ".join(str(d.get(k) or "") for k in (
       "first_reference_authors","first_reference_title","public_collected_by")).lower()
    return bool(re.search(r"\b(kawahara|yamashita|yoshimura|m(i)?yaura|lei|kondoh|azuma)\b",s))

def main():
    OUT.mkdir(exist_ok=True)
    candidates=[]
    for label in ("Daphne_plastid_spacer_2010","Daphne_ITS2_2010"):
        q=search(SEARCHES[label])
        # Keep broad count and only first records of the search.
        ids=q["ids"]
        assert len(ids)==q["count_indexed_records"],"Some historical index hits were not retrieved; abort before completeness claims"
        if not ids:continue
        for i in range(0,len(ids),12):
            batch=ids[i:i+12]
            sub=extract(fetch_flat(batch))
            if len(sub)!=len(batch):
                raise ValueError("GenBank flat-file records missing: "+label+" "+str((len(sub),len(batch))))
            for row in sub:
                row["source_search_label"]=label
                row["team_author_match_in_first_reference"]=team_match(row)
                candidates.append(row)
            time.sleep(.55)
        time.sleep(.5)
    assert len(candidates)>=60
    summary={
       "status":"GenBank_first_record_provenance_original_Daphne_2010_v43",
       "barcoding_presentation_doi":"10.11519/jfsc.121.0.332.0",
       "sequence_record_count_inspected":len(candidates),
       "records_with_any_team_first_reference_name_match":
           sum(v["team_author_match_in_first_reference"] for v in candidates),
       "record_metadata":candidates,
       "record_first_reference_title_tally":{name:sum(x["first_reference_title"]==name for x in candidates)
            for name in sorted({x["first_reference_title"] for x in candidates}, key=str)},
       "team_matching_accessions":[x["accession"] for x in candidates if x["team_author_match_in_first_reference"]],
       "strict_guard":{
         "original_seven_barcode_markers_six_plastid_one_multicopy_nuclear":True,
         "first_reference_metadata_not_exhaustive_source_provenance":True,
         "no_original_2010_sequence_accession_claim_without_voucher_crosswalk":True,
         "search_subset_not_all_historical_NCBI_records":True,
         "submission_year_is_not_conference_year":True,
         "not_2010_low_copy_independent_nuclear_loci":True,
         "no_focal_Onishibari_new_genetic_measurement":True}}
    (OUT/"historical_daphne_2010_barcode_record_provenance_v43.json").write_text(
       json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_2010_NCBI_GENBANK_PROVENANCE_V43",json.dumps({
        "status":summary["status"],
        "n_records_inspected":len(candidates),
        "n_2010_team_author_name_first_reference":summary["records_with_any_team_first_reference_name_match"],
        "sample_candidates":[{"acc":a["accession"],"org":a["organism"],
          "authors":a["first_reference_authors"],"title":a["first_reference_title"]}
          for a in candidates[:18]],
        "team_name_matching_accessions":summary["team_matching_accessions"],
        "reference_title_tally":summary["record_first_reference_title_tally"],
        "not_confirmed_2010_source_without_voucher":True
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
