#!/usr/bin/env python3
"""Same-individual nrITS recovery from Yoo et al. 2023 chloroplast-library raw reads.

CRITICAL:
- ON244034 plastome source SRA: SRR19117239, voucher Lee 7891, Gangneung.
- 2025 C nrITS accession PX056848 is Pyeongchang, a DIFFERENT specimen.
- SRA was prepared from TOTAL leaf DNA, but nrITS read presence is unknown.
- ITS is multicopy: within-specimen mismatches may be repeat variation, not heterozygosity.
- Bowtie2 reference recruitment may miss divergent repeat copies. Do not infer
  ILS or chloroplast capture from one rDNA consensus.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys
import pysam
from daphne_its_audit import fetch, write_fasta, read_fasta, ACGT

OUT = Path("analysis_results_same_individual")
RUN = "SRR19117239"
TARGET = {"P": "PX056852", "C2025": "PX056848", "K": "PX056858", "O": "MT623684"}

def main():
    OUT.mkdir(exist_ok=True)
    refs = {}
    sources = {}
    for alias, acc in TARGET.items():
        header, seq = fetch(acc)
        refs[alias] = seq
        sources[alias] = {
            "accession": acc, "header": header, "bp": len(seq),
            "sha256": hashlib.sha256(seq.encode()).hexdigest()
        }
    target = OUT/"C2025_reference.fa"
    write_fasta(target,[("C2025",refs["C2025"])])
    paired = OUT/"C2025_mapped.bam"
    if not paired.exists():
        raise FileNotFoundError(f"Missing mapped BAM {paired}, run bowtie2 mapping first")
    af = pysam.AlignmentFile(str(paired),"rb")
    counts=[Counter() for _ in refs["C2025"]]
    for col in af.pileup("C2025", truncate=True, stepper="samtools",
         max_depth=100000, min_base_quality=25, min_mapping_quality=15):
        if col.reference_pos >= len(counts): continue
        for read in col.pileups:
            if read.is_del or read.is_refskip: continue
            r=read.alignment
            if r.is_secondary or r.is_supplementary or r.is_duplicate: continue
            p=read.query_position
            if p is None: continue
            base=r.query_sequence[p].upper()
            if base in ACGT: counts[col.reference_pos][base]+=1
    af.close()
    bases=[]
    rows=[]
    for i,(ref,cn) in enumerate(zip(refs["C2025"],counts),1):
        depth=sum(cn.values())
        major,major_n=cn.most_common(1)[0] if depth else ("N",0)
        frac=major_n/depth if depth else 0
        base=major if depth>=8 and frac>=0.90 else "N"
        bases.append(base)
        rows.append(dict(position=i,reference=ref,depth=depth,major=major,
             major_fraction=round(frac,4),A=cn.get("A",0),C=cn.get("C",0),
             G=cn.get("G",0),T=cn.get("T",0),confident_base=base))
    write_fasta(OUT/"C2023_recovered_consensus_masked.fa",[("C2023", "".join(bases))])
    tsv=OUT/"C2023_ITS_base_counts.tsv"
    fields=["position","reference","depth","major","major_fraction","A","C","G","T","confident_base"]
    tsv.write_text("\t".join(fields)+"\n"+"\n".join(
      "\t".join(str(row[f]) for f in fields) for row in rows)+"\n")
    # Validate the six P-C vs K-outgroup diagnostic nrITS sites with same-individual 2023 reads.
    quartet=OUT/"P_C2025_K_O_raw.fa"
    aln=OUT/"P_C2025_K_O_aligned.fa"
    write_fasta(quartet,[(x,refs[x]) for x in ("P","C2025","K","O")])
    with aln.open("w") as o:
        subprocess.run(["mafft","--auto","--quiet",str(quartet)],check=True,stdout=o,timeout=180)
    aligned=read_fasta(aln)
    cposition=0
    diagnostic=[]
    for col, (p,c,k,o) in enumerate(zip(*(aligned[x] for x in ("P","C2025","K","O"))),1):
        if c!='-': cposition+=1
        if p==c and k==o and p!=k and all(v in ACGT for v in (p,c,k,o)):
            d=rows[cposition-1]
            diagnostic.append({"alignment_column":col,"its_position":cposition,
              "P_C2025":p,"K_outgroup":k,
              "C2023_depth":d["depth"],"C2023_counts":{b:d[b] for b in "ACGT"},
              "C2023_major":d["major"],"C2023_major_fraction":d["major_fraction"],
              "C2023_confident_call":d["confident_base"],
              "matches_P_C2025":d["confident_base"]==p,
              "matches_K":d["confident_base"]==k})
    all_cov=[r["depth"] for r in rows]
    confident=[r for r in rows if r["confident_base"] in ACGT]
    mismatch=[r for r in confident if r["confident_base"]!=r["reference"]]
    diag_known=[d for d in diagnostic if d["C2023_confident_call"] in ACGT]
    samstat=subprocess.check_output(["samtools","flagstat",str(paired)],text=True)
    stats={
      "source_2023": {"run":RUN,"plastome":"ON244034","voucher":"Lee 7891",
        "site":"Gangneung, South Korea","library":"total DNA, Illumina NovaSeq 6000",
        "sra_total_bases":359842060},
      "independent_C2025_nrITS_source":{"accession":"PX056848","voucher":"Byeon-Dpk250701",
        "site":"Pyeongchang, South Korea"},
      "references":sources,
      "nucleotide_positions":len(rows),
      "positions_depth_at_least_1":sum(d>=1 for d in all_cov),
      "positions_depth_at_least_8":sum(d>=8 for d in all_cov),
      "positions_depth_at_least_20":sum(d>=20 for d in all_cov),
      "median_depth":statistics.median(all_cov),
      "mean_depth":round(statistics.mean(all_cov),3),
      "high_confidence_calls_depth8_major90pct":len(confident),
      "high_confidence_consensus_mismatch_vs_C2025":len(mismatch),
      "high_confidence_mismatch_positions":[r["position"] for r in mismatch],
      "PC_vs_K_outgroup_diagnostic_sites":diagnostic,
      "confident_diagnostic_sites":len(diag_known),
      "diagnostic_same_PC2025":sum(x["matches_P_C2025"] for x in diag_known),
      "diagnostic_same_K":sum(x["matches_K"] for x in diag_known),
      "other_known_diagnostic":sum(not x["matches_P_C2025"] and not x["matches_K"] for x in diag_known),
      "samtools_flagstat":samstat,
      "warnings":[
        "Nuclear rDNA is a multicopy locus: consensus does not indicate diploid genotype",
        "Reference recruitment can bias detection of divergent paralogous rDNA copies",
        "The 2023 material is from the ON244034 plastome individual, not the Han 2025 nrITS individual",
        "No nuclear species tree or introgression model was inferred"],
    }
    (OUT/"same_individual_nrits_audit.json").write_text(
       json.dumps(stats,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:v for k,v in stats.items() if k not in ["references","samtools_flagstat","PC_vs_K_outgroup_diagnostic_sites"]},indent=2))
    print("DIAGNOSTIC SITES", json.dumps(diagnostic,indent=2))

if __name__=="__main__":
    main()
