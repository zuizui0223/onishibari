#!/usr/bin/env python3
"""Evaluate focal ITS quartet under three competing unrooted topologies.

This is a MARKER GENE TREE TEST, not introgression-vs-ILS inference.
P/C/K = Japanese D. pseudomezereum / Korean var koreana / D. kamtschatica.
Rooting with the fourth sequence can imply an ancestral split, not a full
species tree with genome-wide support.
"""
from pathlib import Path
import json
import subprocess
from daphne_its_audit import fetch, write_fasta, read_fasta

OUT=Path("analysis_results_topology")
OUT.mkdir(exist_ok=True)
FOCAL_ITS={"P":"PX056852","C":"PX056848","K":"PX056858","O":"MT623684"}
FOCAL_CP={"P":"PV752465","C":"ON244034","K":"PV752464","O":"MN201546"}
TREES="""((P,C),(K,O));
((C,K),(P,O));
((P,K),(C,O));
"""
def run(command,log):
    with log.open("w") as fh:
        subprocess.run(command,check=True,stdout=fh,stderr=subprocess.STDOUT,timeout=1200)
def analyze(marker,accessions,model):
    d=OUT/marker; d.mkdir(exist_ok=True)
    seqs=[]
    for label,accession in accessions.items():
        header,seq=fetch(accession)
        seqs.append((label,seq))
    raw=d/"input.fa"; aln=d/"alignment.fa"; alt=d/"three_topologies.tree"
    write_fasta(raw,seqs)
    with aln.open("w") as handle:
        subprocess.run(["mafft","--auto","--quiet","--thread","2",str(raw)],
            check=True,stdout=handle,stderr=subprocess.PIPE,timeout=1200)
    alt.write_text(TREES)
    # Per IQ-TREE topology testing documentation:
    # -z list of user trees; -n 0 evaluate user trees only; -zb RELL replicates; -au AU.
    run(["iqtree2","-s",str(aln),"-m",model,"-z",str(alt),
         "-n","0","-zb","10000","-au","-T","2","-pre",str(d/"topotest")],
        d/"topotest_stdout.txt")
    ir_control = None
    if marker == "plastome":
        aligned = read_fasta(aln)
        # Yoo et al. 2023: ON244034 is ordered LSC(84,963), IRb(41,725),
        # SSC(2,739), IRa(41,725). Remove final IRa's C-reference columns.
        ncp = len(aligned["C"].replace("-", ""))
        if ncp != 171152:
            raise ValueError("Unexpected C plastome length for IR coordinates: "+str(ncp))
        pos = 0
        keep = []
        for i,base in enumerate(aligned["C"]):
            if base != "-":
                pos += 1
            if pos <= 129427:
                keep.append(i)
        masked = d/"alignment_IRa_excluded.fa"
        write_fasta(masked,[(taxon,"".join(seq[i] for i in keep))
                            for taxon,seq in aligned.items()])
        run(["iqtree2","-s",str(masked),"-m",model,"-z",str(alt),
             "-n","0","-zb","10000","-au","-T","2",
             "-pre",str(d/"topotest_IRa_excluded")],
            d/"topotest_IRa_excluded_stdout.txt")
        ir_control={"C_reference":"ON244034","reference_length":ncp,
            "excluded_C_reference_positions_1based":"129428-171152",
            "excluded_one_inverted_repeat_copy_bp":41725,
            "unaligned_columns_retained":len(keep),
            "caution":"IR coordinates based on published plastome order LSC+IRb+SSC+IRa; cp gene tree still a single linked locus"}
    return {"accessions":accessions,"alignment":str(aln),"model":model,
      "plastome_ir_exclusion_control":ir_control,
      "tree_order":["PC|KO","CK|PO","PK|CO"],
      "note":"A reject/not-reject result concerns one ITS marker or one linked plastome, not ILS/introgression."}
def main():
    summary={"question":"Compare alternative unrooted gene trees within each marker",
             "nuclear_is_multicopy_rDNA":True,"plastome_is_linked_one_genealogy":True,
             "markers":{}}
    for marker,acc,model in [
        ("nrITS",FOCAL_ITS,"GTR+F+I+G4"),
        ("plastome",FOCAL_CP,"K3Pu+F+I+G4")]:
        try:
            summary["markers"][marker]=analyze(marker,acc,model)
        except Exception as e:
            summary["markers"][marker]={"error":type(e).__name__+": "+str(e)}
            (OUT/"index.json").write_text(json.dumps(summary,indent=2)+"\n")
            raise
    (OUT/"index.json").write_text(json.dumps(summary,indent=2)+"\n")

if __name__=="__main__":main()
