#!/usr/bin/env python3
"""Exploratory 4-taxon whole-plastome alignment audit for Han et al. 2025.

Chloroplast DNA is one linked genealogy; site counts do NOT identify
gene flow, chloroplast capture, or incomplete lineage sorting.
Plastomes include duplicated inverted repeat regions; report them as-is,
consistent with initial full-plastome alignment sensitivity checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from daphne_its_audit import fetch,write_fasta,read_fasta,audit_sites

ACCESSIONS = {
    "P": ("PV752465", "Daphne pseudomezereum, Mie, Japan"),
    "C": ("ON244034", "D. kamtschatica var. koreana, Gangneung, Korea (different individual from nrITS)"),
    "K": ("PV752464", "Daphne kamtschatica, Kamchatka, Russia"),
    "O1": ("MN201546", "Daphne laureola, plastome"),
    "O2": ("KY991380", "Daphne kiusiana, Jeju, Korea, plastome"),
}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="analysis_results_cp")
    args=ap.parse_args()
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/"raw_plastomes"
    raw.mkdir(exist_ok=True)
    report=dict(source="Han et al. 2025 doi:10.11110/kjpt.2025.55.3.157 and Yoo et al. 2023",
      reanalyzed=False,accessions={},quartet_audits={},
      caveats=["Single nonrecombining plastid inheritance unit; sites are not independent gene genealogies",
        "Large duplicated IR regions in Daphne can make raw site counts pseudo-replicated",
        "Taxon C nrITS and plastome are from different individuals",
        "Local 4-taxon test is not original full-taxon study or species-tree inference",
        "No identification of introgression or ILS from plastome alone"])
    try:
        seqs={}
        for label,(accession,description) in ACCESSIONS.items():
            header,seq=fetch(accession)
            report["accessions"][label]=dict(accession=accession,origin=description,
                genbank_header=header,n_bases=len(seq),
                sha256=hashlib.sha256(seq.encode()).hexdigest())
            if not (140000 <= len(seq) <= 190000):
                raise ValueError(f"Unexpected plastome length: {accession}, {len(seq)}")
            seqs[label]=seq
            write_fasta(raw/(label+"_"+accession+".fa"),[(label,seq)])
            time.sleep(0.4)
        for og in ("O1","O2"):
            inf=out/("plastome_"+og+"_raw.fa")
            aln=out/("plastome_"+og+"_aligned.fa")
            write_fasta(inf,[(label,seqs[label]) for label in ("P","C","K",og)])
            with aln.open("w",encoding="utf8") as handle:
                subprocess.run(["mafft","--auto","--quiet","--thread","2",str(inf)],
                    check=True,stdout=handle,stderr=subprocess.PIPE,timeout=900)
            report["quartet_audits"][og]=audit_sites(read_fasta(aln),og)
        report["reanalyzed"]=True
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        (out/"audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
        raise
    (out/"audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
    print(json.dumps(report["quartet_audits"],indent=2))

if __name__=="__main__":
    main()
