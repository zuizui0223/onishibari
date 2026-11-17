#!/usr/bin/env python3
"""Independent plastid barcode alignment for *D. jezoensis* and closest Daphne.

Uses voucher nakata074 matK/rbcL partial accessions and three full plastomes
with identity checks. The two linked chloroplast genes are ONE genealogy, not
two independent nuclear loci. Not a formal species tree or sex evolution test.
"""
import argparse,csv,hashlib,io,json,time,urllib.parse,urllib.request
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from Bio import SeqIO, Align

BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PLASTOMES={"PV752465":"Daphne pseudomezereum",
           "PV752464":"Daphne kamtschatica",
           "ON244034":"Daphne pseudomezereum var. koreana"}
GENES=("matK","rbcL")
LABELS={"PV752465":"D_pseudomezereum","PV752464":"D_kamtschatica",
        "ON244034":"koreana","jezoensis":"D_jezoensis"}
EXPECTED_JEZ={"LC680075":"matK","LC693863":"rbcL"}
EXPECTED_DISTANCE={
  ("D_pseudomezereum","D_kamtschatica"):3,
  ("D_pseudomezereum","koreana"):2,
  ("D_pseudomezereum","D_jezoensis"):1,
  ("D_kamtschatica","koreana"):1,
  ("D_kamtschatica","D_jezoensis"):4,
  ("koreana","D_jezoensis"):3
}

def fetch():
    url=BASE+"?"+urllib.parse.urlencode(dict(db="nuccore",id=",".join(PLASTOMES),
                                  rettype="gbwithparts",retmode="text"))
    error=None
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Daphne-plastid-comparison/2026"})
            with urllib.request.urlopen(req,timeout=65) as r:
                text=r.read(9000000).decode("utf-8")
            break
        except Exception as e:
            error=e;time.sleep(2**retry)
    else:raise RuntimeError("Unable to fetch core plastomes: "+str(error))
    fetched={}
    for r in SeqIO.parse(io.StringIO(text),"genbank"):
        acc=r.id.split(".")[0]
        if acc not in PLASTOMES or acc in fetched:
            raise ValueError("Unexpected/duplicate accession "+r.id)
        org=r.annotations.get("organism","")
        if not org.startswith(PLASTOMES[acc]) or (
                acc=="PV752465" and "koreana" in org.lower()):
            raise ValueError("Taxon accession mismatch "+r.id+" "+org)
        genes={}
        for g in GENES:
            seqs={str(f.extract(r.seq)).upper() for f in r.features
                  if f.type=="CDS" and g in f.qualifiers.get("gene",[])}
            if len(seqs)!=1:raise ValueError("Gene missing/discordant "+acc+" "+g)
            genes[g]=seqs.pop()
        fetched[acc]=genes
    if set(fetched)!=set(PLASTOMES):
        raise ValueError("Required chloroplast taxon absent")
    return fetched

def read_jezoensis(path):
    genes={}
    for r in SeqIO.parse(str(path),"fasta"):
        a=r.id.split("|")
        if len(a)<4 or a[0].split(".")[0] not in EXPECTED_JEZ:
            continue
        g=EXPECTED_JEZ[a[0].split(".")[0]]
        if g!=a[1] or a[2]!="Daphne_jezoensis" or a[3]!="nakata074":
            raise ValueError("D. jezoensis voucher/marker mismatch")
        if g in genes:raise ValueError("Repeated Dje locus")
        genes[g]=str(r.seq).upper()
    if set(genes)!=set(GENES):
        raise ValueError("Both same-voucher Dje matK and rbcL required")
    return genes

def map_fragment(reference,partial):
    aligner=Align.PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,
                                 open_gap_score=-8,extend_gap_score=-.5)
    aln=aligner.align(reference,partial)[0]
    positions={}
    for (x0,x1),(y0,y1) in zip(aln.aligned[0],aln.aligned[1]):
        for x,y in zip(range(x0,x1),range(y0,y1)):
            if reference[x] in "ACGT" and partial[y] in "ACGT":
                positions[x]=partial[y]
    if len(positions)!=len(partial):
        raise ValueError("Partial marker has unaligned or ambiguous bases")
    return positions

def audit(fragments,core):
    aln={}
    per_gene={}
    for g in GENES:
        reference=core["PV752465"][g]
        query=fragments[g]
        matched=map_fragment(reference,query)
        sites=sorted(matched)
        seqs={LABELS[a]:"".join(core[a][g][i] for i in sites)
              for a in PLASTOMES}
        seqs["D_jezoensis"]="".join(matched[i] for i in sites)
        if any(len(s)!=len(sites) for s in seqs.values()):
            raise AssertionError("DNA matrix length inconsistent")
        per_gene[g]={
          "sites":len(sites),"reference_start_1based":sites[0]+1,
          "reference_end_1based":sites[-1]+1,
          "variable_sites":[dict(reference_position_1based=sites[i]+1,
                    bases={n:s[i] for n,s in seqs.items()})
                    for i in range(len(sites)) if len({s[i] for s in seqs.values()})>1]
        }
        for name,s in seqs.items():aln.setdefault(name,{})[g]=s
    distances=[]
    for a,b in combinations(LABELS.values(),2):
        genes={g:sum(x!=y for x,y in zip(aln[a][g],aln[b][g]))
               for g in GENES}
        distances.append(dict(taxon_a=a,taxon_b=b,matK_diff=genes["matK"],
             rbcL_diff=genes["rbcL"],combined_diff=sum(genes.values()),
             common_bp=sum(x["sites"] for x in per_gene.values())))
    return per_gene,distances,aln

def selftest():
    pos=map_fragment("TTACGATCCG","ACGATC")
    assert sorted(pos)==[2,3,4,5,6,7]
    assert "".join(pos[x] for x in sorted(pos))=="ACGATC"
    print("JEZOENSIS_PLASTID_COMPARISON_SELFTEST=SUCCESS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--jezoensis",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args();selftest()
    if args.self_test:return
    if not args.jezoensis or not args.out:p.error("Require --jezoensis and --out")
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError("Use empty output directory")
    args.out.mkdir(parents=True,exist_ok=True)
    genes=read_jezoensis(args.jezoensis)
    core=fetch()
    per,rows,seqs=audit(genes,core)
    expected_bp=1354
    if sum(v["sites"] for v in per.values())!=expected_bp:
        raise ValueError("Shared marker length drift; verify accession version")
    for r in rows:
        k=(r["taxon_a"],r["taxon_b"])
        if k not in EXPECTED_DISTANCE or r["combined_diff"]!=EXPECTED_DISTANCE[k]:
            raise ValueError("Unexpected barcode distance; source drift "+str(r))
    with (args.out/"Daphne_core4_matK_rbcL_distances.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    for g in GENES:
        with (args.out/(g+"_four_taxa_common.fasta")).open("w") as f:
            for name in sorted(seqs):
                f.write(">"+name+"\n"+seqs[name][g]+"\n")
    report=dict(n_taxa=4,common_total_bp=expected_bp,genes=per,
                comparison=rows,
                accession_names=dict(core=PLASTOMES,
                    jezoensis={"LC680075":"matK","LC693863":"rbcL"}),
                note="Only 2 linked chloroplast gene fragments from a single Dje voucher. Low polymorphism cannot resolve species tree, selfing evolution or introgression.")
    (args.out/"Daphne_jezoensis_plastid_barcode_report.json").write_text(
          json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("JEZOENSIS_CHLOROPLAST_MATRIX="+json.dumps(
        dict(common_bp=expected_bp,gene_lengths={g:x["sites"] for g,x in per.items()},
             distances=[(x["taxon_a"],x["taxon_b"],x["combined_diff"]) for x in rows]),
        ensure_ascii=False),flush=True)
if __name__=="__main__":main()
