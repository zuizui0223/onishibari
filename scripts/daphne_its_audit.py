#!/usr/bin/env python3
"""Exploratory Daphne nrITS quartet audit; does NOT distinguish ILS from introgression."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import urllib.request
import urllib.parse

# Focal accessions: Han et al. 2025; D. kiusiana outgroups: Lee et al. 2022.
ACCESSIONS = {
    "P": ("PX056852", 550, 650, "Daphne pseudomezereum, Japan"),
    "C": ("PX056848", 550, 650, "D. kamtschatica var. koreana, Korea"),
    "K": ("PX056858", 550, 650, "Daphne kamtschatica, Russia"),
    "O1": ("MT623684", 400, 900, "Daphne kiusiana, Kumamoto, Japan"),
    "O2": ("MT623688", 400, 900, "Daphne kiusiana, Jeju, Korea"),
    "O3": ("MT623695", 400, 900, "Daphne odora, Guangdong, China"),
    "O4": ("MT623699", 400, 900, "Daphne papyracea, Guangdong, China"),
}
IUPAC = set("ACGTRYSWKMBDHVN")
ACGT = set("ACGT")

def fetch(accession):
    q = urllib.parse.urlencode(dict(db="nuccore", id=accession, rettype="fasta", retmode="text"))
    req = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + q,
        headers={"User-Agent": "Daphne-systematics-audit/0.2 (scientific)"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                s = r.read().decode("utf8")
            lines = [x.strip() for x in s.splitlines() if x.strip()]
            if not lines or not lines[0].startswith(">") or accession not in lines[0]:
                raise ValueError("Unexpected FASTA header for " + accession + ": " + s[:100])
            seq = "".join(lines[1:]).upper().replace("U","T")
            if len(seq) == 0 or not set(seq) <= IUPAC:
                raise ValueError("Invalid nucleotide string for " + accession)
            return lines[0], seq
        except (OSError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)

def write_fasta(p, pairs):
    p.write_text("".join(">" + label + "\n" +
        "\n".join(seq[i:i+80] for i in range(0,len(seq),80)) + "\n"
        for label,seq in pairs), encoding="utf8")

def read_fasta(p):
    result = {}
    label = None
    for line in p.read_text(encoding="utf8").splitlines():
        if line.startswith(">"):
            label = line[1:].split()[0]
            if label in result: raise ValueError("Duplicate FASTA header")
            result[label] = ""
        elif label is not None:
            result[label] += line.strip().upper()
    return result

def audit_sites(aln, out):
    ids = ("P","C","K",out)
    if set(aln) != set(ids) or len(set(map(len,aln.values()))) != 1:
        raise ValueError("Invalid 4-taxon alignment")
    informative = {"PC":0,"CK":0,"PK":0}
    complete = variable = 0
    for p,c,k,o in zip(*(aln[x] for x in ids)):
        if not {p,c,k,o} <= ACGT: continue
        complete += 1
        if len({p,c,k,o}) > 1: variable += 1
        if p==c and k==o and p!=k: informative["PC"] += 1
        if c==k and p==o and c!=p: informative["CK"] += 1
        if p==k and c==o and p!=c: informative["PK"] += 1
    pairs = {}
    for a,b in [("P","C"),("P","K"),("C","K")]:
        sites = [(x,y) for x,y in zip(aln[a],aln[b]) if x in ACGT and y in ACGT]
        pairs[a+"-"+b] = dict(compared_sites=len(sites),
            identity=(round(sum(x==y for x,y in sites)/len(sites),6) if sites else None))
    return dict(outgroup=out, alignment_length=len(aln["P"]), complete_acgt=complete,
                variable_complete_acgt=variable, quartet_2plus2_sites=informative,
                pairwise_ungapped_identity=pairs)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out", default="analysis_results")
    args=ap.parse_args()
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    raw=out/"raw_its"
    raw.mkdir(exist_ok=True)
    report={
      "source":"Han et al. 2025 doi:10.11110/kjpt.2025.55.3.157",
      "reanalyzed":False,"accessions":{},"quartet_audits":{},
      "limitations":["One focal individual per taxon",
         "nrITS is multicopy and not independent low-copy nuclear loci",
         "The published C plastome is from a different individual than the C ITS",
         "Four-taxon sensitivity test is not a reproduction of the original 12 Daphne + outgroup ML analyses",
         "No introgression-versus-ILS inference is possible from this dataset alone"]}
    try:
        seqs={}
        for key,(acc,mn,mx,description) in ACCESSIONS.items():
            header,seq=fetch(acc)
            report["accessions"][key]=dict(accession=acc,origin=description,
                genbank_header=header,length=len(seq),expected_length_pass=mn<=len(seq)<=mx,
                sha256=hashlib.sha256(seq.encode()).hexdigest())
            if not mn<=len(seq)<=mx:
                raise ValueError("Unexpected length for "+acc+": "+str(len(seq)))
            seqs[key]=seq
            write_fasta(raw/(key+"_"+acc+".fa"),[(key,seq)])
            time.sleep(0.4)
        for og in ("O1","O2","O3","O4"):
            input_file=out/("quartet_"+og+"_raw.fa")
            aligned_file=out/("quartet_"+og+"_aligned.fa")
            write_fasta(input_file,[(x,seqs[x]) for x in ("P","C","K",og)])
            with aligned_file.open("w",encoding="utf8") as fh:
                subprocess.run(["mafft","--auto","--quiet",str(input_file)],
                    check=True,stdout=fh,stderr=subprocess.PIPE,timeout=180)
            report["quartet_audits"][og]=audit_sites(read_fasta(aligned_file),og)
        report["reanalyzed"]=True
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
        (out/"audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")
        raise
    (out/"audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf8")

if __name__=="__main__":
    main()
