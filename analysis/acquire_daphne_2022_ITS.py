#!/usr/bin/env python3
"""Accession-verified targeted retrieval of published 2022 Daphne nrITS records.

Source: Lee et al. (2022) Plant Diversity Table 1 (19 accessions; not 19 species).
No inferred sexuality, phylogenetic species-tree support or sample linkage.
"""
import argparse,csv,hashlib,io,json,time,urllib.parse,urllib.request
from datetime import datetime,timezone
from pathlib import Path
from Bio import SeqIO

GROUPS={
 "Daphne kiusiana var. kiusiana (Kumamoto)":(list(range(684,688)),"Daphne kiusiana"),
 "Daphne kiusiana var. kiusiana (Jeju)":(list(range(688,692)),"Daphne kiusiana"),
 "Daphne kiusiana var. atrocaulis (Hunan)":(list(range(692,695)),"Daphne kiusiana"),
 "Daphne odora (Guangdong)":(list(range(695,699)),"Daphne odora"),
 "Daphne papyracea (Guangdong)":(list(range(699,703)),"Daphne papyracea"),
}
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DOI="https://doi.org/10.1016/j.pld.2021.11.002"

def accession_to_label():
    result={}
    for label,(nums,expected) in GROUPS.items():
        for n in nums:
            acc="MT623"+str(n)
            if acc in result: raise ValueError("Duplicate accession")
            result[acc]=(label,expected)
    assert len(result)==19
    return result

def validate_record(rec,accession,expected):
    if rec.id.split(".")[0]!=accession:
        raise ValueError("Wrong accession record: "+rec.id)
    organism=rec.annotations.get("organism","")
    if not organism.startswith(expected):
        raise ValueError(f"Taxon mismatch {rec.id}: {organism}")
    desc=rec.description.lower()
    if not any(term in desc for term in ("internal transcribed","ribosomal","its1","its2")):
        raise ValueError(f"Not annotated as nuclear ribosomal ITS: {rec.id}")
    seq=str(rec.seq).upper()
    if not 350<=len(seq)<=1200 or set(seq)-set("ACGTRYSWKMBDHVN"):
        raise ValueError(f"Length/characters outside ITS assay scope: {rec.id}")
    return {"accession":rec.id,"organism":organism,"length":len(seq),
            "description":rec.description,"seq":seq}

def fetch(accessions):
    url=BASE+"?"+urllib.parse.urlencode({
        "db":"nuccore","id":",".join(accessions),
        "rettype":"gb","retmode":"text"})
    last=None
    for tries in range(3):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"onishibari-phylogeny-literature-audit/2026"})
            with urllib.request.urlopen(req,timeout=60) as s:
                return s.read().decode("utf-8"),s.geturl()
        except Exception as exc:
            last=exc
            time.sleep(2**tries)
    raise RuntimeError(str(last))

def test():
    mapping=accession_to_label()
    assert mapping["MT623684"][0].endswith("Kumamoto)")
    assert mapping["MT623690"][0].endswith("Jeju)")
    assert mapping["MT623693"][0].endswith("Hunan)")
    assert mapping["MT623701"][1]=="Daphne papyracea"
    assert "MT623683" not in mapping
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec=SeqRecord(Seq("A"*550),id="MT623684.1",
                  description="internal transcribed spacer of Daphne kiusiana")
    rec.annotations["organism"]="Daphne kiusiana var. kiusiana"
    assert validate_record(rec,"MT623684","Daphne kiusiana")["length"]==550
    try:validate_record(rec,"MT623684","Daphne pseudomezereum")
    except ValueError:pass
    else:raise AssertionError("Taxon guard failed")
    print("LEE2022_ITS_SELFTEST=SUCCESS mapping 19 accessions and accession/taxon guard")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,
                   default=Path("data/derived/lee2022_ITS_evidence"))
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    test()
    if args.self_test:return
    output=args.out
    if output.exists() and any(output.iterdir()):
        raise SystemExit("Refuse to overwrite existing dated evidence snapshot")
    output.mkdir(parents=True,exist_ok=True)
    targets=accession_to_label()
    recs,error=[],[]
    # Split into three limited batches to guard individual NCBI size/rate limits.
    items=list(targets)
    for offset in range(0,len(items),7):
        batch=items[offset:offset+7]
        raw,url=fetch(batch)
        found={}
        for rec in SeqIO.parse(io.StringIO(raw),"genbank"):
            acc=rec.id.split(".")[0]
            if acc not in batch:
                raise ValueError("Unexpected accession returned: "+rec.id)
            if acc in found: raise ValueError("Duplicate downloaded accession "+acc)
            found[acc]=rec
        for acc in batch:
            try:
                if acc not in found: raise ValueError("Accession not returned")
                group,expected=targets[acc]
                parsed=validate_record(found[acc],acc,expected)
                parsed["study_group"]=group
                parsed["source_url"]=url
                parsed["source_paper"]=DOI
                parsed["sha256"]=hashlib.sha256(parsed["seq"].encode()).hexdigest()
                recs.append(parsed)
            except Exception as e:
                error.append({"accession":acc,"error":str(e)})
        time.sleep(.5)
    if len(recs)+len(error)!=19: raise AssertionError("Loss of accession audit row")
    with (output/"published_2022_ITS.fasta").open("w") as f:
        for r in recs:
            f.write(">"+r["accession"]+"|"+r["organism"].replace(" ","_")+"\n")
            seq=r["seq"]
            f.write("\n".join(seq[x:x+80] for x in range(0,len(seq),80))+"\n")
    header=["accession","organism","length","description","study_group",
            "source_url","source_paper","sha256"]
    with (output/"inventory.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=header)
        writer.writeheader()
        for r in recs:
            writer.writerow({k:r[k] for k in header})
    summary={"retrieved_utc":datetime.now(timezone.utc).isoformat(),
             "source_paper":DOI,"expected_accessions":19,
             "verified_accessions":len(recs),"failed":error,
             "biology_note":"Registered accessions are individuals, not distinct species; 2022 Kumamoto and Jeju observations must not be equated to 2023 eastern/western selfing estimates",
             "no_tree_or_trait_assignment":True}
    (output/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print("LEE2022_ITS_RESULTS="+json.dumps(summary,ensure_ascii=False),flush=True)
    if not recs: raise SystemExit("No verified accessions returned; do not claim acquisition")
if __name__=="__main__": main()
