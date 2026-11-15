#!/usr/bin/env python3
"""Daphne ITS taxon-sampling stress test using GenBank series MT623670–MT623710.

This is an audit of ITS gene tree stability, not a phylogenomic species tree.
Series membership is validated by GenBank FASTA headers. Historic names and
duplicate species are tracked; only one sequence per nonfocal species.
"""
import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path
from Bio import Phylo
from Bio import SeqIO
from daphne_its_audit import fetch,write_fasta,read_fasta

OUT=Path("analysis_results_multitaxon_its")
OUT.mkdir(exist_ok=True)
IUPAC=set("ACGTRYSWKMBDHVN")
KNOWN={
    "P":("PX056852","Daphne pseudomezereum"),
    "C":("PX056848","Daphne kamtschatica var. koreana"),
    "K":("PX056858","Daphne kamtschatica"),
}
PREFER={
    "kiusiana":"MT623684", "odora":"MT623695","papyracea":"MT623699"
}
def fetch_batch():
    ids=",".join("MT"+str(i) for i in range(623670,623711))
    params=urllib.parse.urlencode(dict(db="nuccore",id=ids,rettype="fasta",retmode="text"))
    req=urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+params,
       headers={"User-Agent":"Daphne-ITS-taxon-audit/0.3 (academic research)"})
    for k in range(4):
        try:
            with urllib.request.urlopen(req,timeout=120) as r: content=r.read().decode("utf8")
            break
        except OSError:
            if k==3: raise
            time.sleep(2**k)
    out={}
    for record in SeqIO.parse(__import__("io").StringIO(content),"fasta"):
        m=re.match(r"^(MT\d+)\.\d+\s+Daphne\s+([a-z-]+)\b",record.description)
        if not m: continue
        accession,species=m.groups()
        seq=str(record.seq).upper().replace("U","T")
        if not 600<=len(seq)<=1050 or not set(seq)<=IUPAC:
            continue
        if "internal transcribed spacer" not in record.description.lower():continue
        out.setdefault(species,[]).append((accession,seq,record.description))
    return out
def filter_alignment(aligned):
    seqs=read_fasta(aligned)
    n=len(seqs); length={len(v) for v in seqs.values()}
    if len(length)!=1:raise ValueError("not aligned")
    keep=[i for i in range(next(iter(length)))
       if sum(s[i] in "ACGT" for s in seqs.values())>=.8*n]
    p=aligned.with_name(aligned.stem+"_occupancy80.fa")
    write_fasta(p,[(k,"".join(s[i] for i in keep)) for k,s in seqs.items()])
    return p,len(keep),next(iter(length))
def run_tree(aligned,prefix):
    subprocess.run(["iqtree2","-s",str(aligned),"-m","MFP","-B","1000","-alrt","1000",
        "-T","2","-pre",str(prefix)],check=True,timeout=480)
    tree=Phylo.read(str(prefix)+".treefile","newick")
    taxa=[t.name for t in tree.get_terminals()]
    anchor="D_odora"
    if anchor not in taxa: anchor=next(t for t in taxa if t not in ("P","C","K"))
    tree.root_with_outgroup(anchor)
    result={"rooting_anchor_for_display_only":anchor,"n_taxa":len(taxa),
       "P_C_are_cherry":False,"C_K_are_cherry":False,"P_K_are_cherry":False,
       "ml_treefile":str(prefix)+".treefile",
       "iqtree_report":str(prefix)+".iqtree"}
    for x,y in (("P","C"),("C","K"),("P","K")):
        clade=tree.common_ancestor(x,y)
        tips=sorted(t.name for t in clade.get_terminals())
        result[x+"_"+y+"_are_cherry"]=tips==sorted([x,y])
        result[x+"_"+y+"_MRCA_tips"]=tips
        result[x+"_"+y+"_MRCA_support"]=str(clade.confidence or clade.name)
    return result

def main():
    report={"research_question":"Does addition of independent named Daphne ITS records affect the P,C,K gene-tree positions?",
      "locus":"nrITS multicopy rDNA only; NOT low-copy nuclear loci",
      "candidate_search":"GenBank MT623670–MT623710 queried by nucleotide record, validated exact headers",
      "genbank_series":{},"retained_accessions":{},"results":{}}
    try:
        series=fetch_batch()
        report["genbank_series"]={s:[a for a,seq,h in values] for s,values in sorted(series.items())}
        if len(series)<8:raise ValueError("Fewer than 8 distinct Daphne species in accession series")
        selected=[]
        for label,(acc,desc) in KNOWN.items():
            h,s=fetch(acc)
            selected.append((label,s))
            report["retained_accessions"][label]={"accession":acc,"description":h}
        for species,choices in sorted(series.items()):
            if species in ("kamtschatica","pseudomezereum"):continue
            acc,seq,head=sorted(choices,key=lambda x:(x[0]!=PREFER.get(species,""),x[0]))[0]
            label="D_"+species
            if label in ("P","C","K"):continue
            selected.append((label,seq))
            report["retained_accessions"][label]={"accession":acc,"description":head}
        if not any(lab=="D_odora" for lab,seq in selected):
            h,s=fetch("MT623695")
            selected.append(("D_odora",s))
            report["retained_accessions"]["D_odora"]={"accession":"MT623695","description":h}
        unaligned=OUT/"daphne_multitaxon_raw.fa"
        aligned=OUT/"daphne_multitaxon_aligned.fa"
        write_fasta(unaligned,selected)
        with aligned.open("w") as handle:
            subprocess.run(["mafft","--auto","--quiet","--thread","2",str(unaligned)],
                check=True,stdout=handle,stderr=subprocess.PIPE,timeout=500)
        filtered,kept,total=filter_alignment(aligned)
        report["n_named_taxa"]=len(selected)
        report["alignment_cols"]=total
        report["occupancy_80pct_alignment_cols"]=kept
        for name,aln in (("all_sites",aligned),("occupancy_80pct",filtered)):
            report["results"][name]=run_tree(aln,OUT/name)
        report["complete"]=True
    except Exception as e:
        report["complete"]=False;report["error"]=type(e).__name__+": "+str(e)
        (OUT/"taxon_sampling_audit.json").write_text(json.dumps(report,indent=2)+"\n")
        raise
    (OUT/"taxon_sampling_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report["results"],indent=2))
if __name__=="__main__":main()
