#!/usr/bin/env python3
"""Build an accession-verified, reference-anchored 24-record ITS audit.

One multicopy nuclear rDNA ITS region, NOT an independent multi-locus species
tree; specimen-level DNA is NOT linked to floral morphology or selfing rates.
"""
import argparse,csv,hashlib,io,json,time,urllib.parse,urllib.request
from collections import defaultdict,Counter
from pathlib import Path
from Bio import Align,SeqIO

CORE={
 "PX056852":"Daphne pseudomezereum",
 "PX056858":"Daphne kamtschatica",
 "PX056848":"Daphne pseudomezereum",
 "OQ756479":"Daphne arisanensis",
 "OQ756480":"Daphne arisanensis",
}
SRC="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def read_lee(path):
    result={}
    for r in SeqIO.parse(str(path),"fasta"):
        acc=r.id.split("|")[0]
        if acc in result:raise ValueError("Duplicate accession in Lee 2022 file")
        result[acc]=str(r.seq).upper()
    if len(result)!=19:raise ValueError("Expected 19 published Lee 2022 accession records")
    if not all(x.startswith("MT623") for x in result):
        raise ValueError("Non-Lee accession found")
    return result

def retrieve_core():
    url=SRC+"?"+urllib.parse.urlencode(dict(db="nuccore",id=",".join(CORE),
                            rettype="gb",retmode="text"))
    err=None
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Daphne-ITS-2026-11-audit"})
            with urllib.request.urlopen(req,timeout=60) as stream:
                txt=stream.read(4000000).decode("utf-8")
            break
        except Exception as exc:
            err=exc
            time.sleep(2**retry)
    else:raise RuntimeError("NCBI efetch failed: "+str(err))
    records=list(SeqIO.parse(io.StringIO(txt),"genbank"))
    found={}
    for r in records:
        acc=r.id.split(".")[0]
        if acc not in CORE or acc in found:
            raise ValueError("Unexpected/duplicate GenBank accession")
        org=r.annotations.get("organism","")
        if not org.startswith(CORE[acc]):
            raise ValueError("Wrong organism "+r.id+": "+org)
        if acc=="PX056852" and "koreana" in org.lower():
            raise ValueError("Target cannot be substituted with Korean variety")
        found[r.id]=str(r.seq).upper()
    if {x.split(".")[0] for x in found}!=set(CORE):
        raise ValueError("Missing any required core accession")
    return found

def assign_group(acc):
    core=acc.split(".")[0]
    if core=="PX056852":return "Dpm"
    if core=="PX056858":return "Dkamtschatica"
    if core=="PX056848":return "koreana"
    if core.startswith("OQ756"):return "D_arisanensis"
    num=int(core.removeprefix("MT623"))
    if 684<=num<=687:return "D_kiusiana_Kumamoto"
    if 688<=num<=691:return "D_kiusiana_Jeju"
    if 692<=num<=694:return "D_kiusiana_atrocaulis_Hunan"
    if 695<=num<=698:return "D_odora"
    if 699<=num<=702:return "D_papyracea"
    raise ValueError("Unexpected accession: "+acc)

def align_common(records,reference):
    aligner=Align.PairwiseAligner(mode="global",match_score=2,mismatch_score=-1,
                                 open_gap_score=-10,extend_gap_score=-.5)
    mapped={reference:list(records[reference])}
    for key,seq in records.items():
        if key==reference:continue
        a=aligner.align(records[reference],seq)[0]
        mapped[key]=[alt for ref,alt in zip(str(a[0]),str(a[1])) if ref!="-"]
        if len(mapped[key])!=len(records[reference]):
            raise AssertionError("Bad coordinate alignment")
    keep=[i for i in range(len(records[reference]))
          if all(mapped[k][i] in "ACGT" for k in mapped)]
    return {k:"".join(mapped[k][i] for i in keep) for k in records},keep

def cluster(records):
    d=defaultdict(list)
    for acc,seq in records.items():
        d[seq].append(acc)
    ordered=sorted(d.items(),key=lambda p:(-len(p[1]),sorted(p[1])))
    mapping={}
    patterns=[]
    for i,(seq,accs) in enumerate(ordered,1):
        hap=f"H{i}"
        groups=Counter(assign_group(x) for x in accs)
        patterns.append(dict(haplotype=hap,count=len(accs),
                   groups=dict(groups),accessions=sorted(accs),
                   sha256=hashlib.sha256(seq.encode()).hexdigest()))
        mapping.update({x:hap for x in accs})
    return patterns,mapping

def selftest():
    records={"PX056852.1":"AAAA","PX056858.1":"AAAT",
             "MT623684.1":"TTTT","MT623688.1":"TTTT"}
    patterns,mapping=cluster(records)
    assert len(patterns)==3
    assert mapping["MT623684.1"]==mapping["MT623688.1"]
    assert patterns[0]["count"]==2
    assert assign_group("MT623688.1")=="D_kiusiana_Jeju"
    assert assign_group("MT623692.1")=="D_kiusiana_atrocaulis_Hunan"
    assert assign_group("PX056852.1")=="Dpm"
    print("DAPHNE_ITS24_COLLAPSE_SELFTEST=SUCCESS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--lee",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args();selftest()
    if args.self_test:return
    if args.lee is None or args.out is None:
        p.error("Need --lee published_2022_ITS.fasta and --out")
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Use new empty output folder for evidence snapshot")
    args.out.mkdir(parents=True,exist_ok=True)
    records=read_lee(args.lee)
    records.update(retrieve_core())
    if len(records)!=24:raise ValueError("Unexpected accessions after merging")
    ref=next(r for r in records if r.split(".")[0]=="PX056852")
    aligned,sites=align_common(records,ref)
    patterns,mapping=cluster(aligned)
    out=args.out
    with (out/"ITS24_reference_common.fasta").open("w") as fh:
        for acc,s in sorted(aligned.items()):
            fh.write(">"+acc+"|"+assign_group(acc)+"|"+mapping[acc]+"\n"+s+"\n")
    with (out/"ITS24_accession_haplotypes.csv").open("w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=("accession","study_group","haplotype"))
        writer.writeheader()
        for acc in sorted(records):
            writer.writerow(dict(accession=acc,study_group=assign_group(acc),
                                 haplotype=mapping[acc]))
    groups=defaultdict(set)
    for acc,hap in mapping.items():
        groups[assign_group(acc)].add(hap)
    kiugroups=("D_kiusiana_Kumamoto","D_kiusiana_Jeju",
               "D_kiusiana_atrocaulis_Hunan")
    shared=set.intersection(*(groups[g] for g in kiugroups))
    info=dict(input_lee2022=19,verified_core=5,accessions=len(records),
        original_reference_bp=len(records[ref]),common_bp=len(sites),
        unique_common_sequence_patterns=len(patterns),
        pattern_multiplicities=sorted([p["count"] for p in patterns],reverse=True),
        n_d_kiusiana_accessions=sum(assign_group(a) in kiugroups for a in records),
        n_d_kiusiana_unique_patterns=len(set.union(*(groups[g] for g in kiugroups))),
        shared_pattern_across_three_kiusiana_groups=sorted(shared),
        patterns=patterns,
        caution="Many accessions have identical 500-bp ITS; cannot treat 24 as 24 different alleles or independent species; regional samples have no measured mating phenotype",
        method="Biopython 1.86 PairwiseAligner global reference PX056852, ACGT positions shared by all records, gap-open -10, gap extension -0.5")
    (out/"ITS24_haplotype_collapse_summary.json").write_text(json.dumps(info,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_ITS24_HAPLOTYPES="+json.dumps({
        k:info[k] for k in ("accessions","common_bp",
          "unique_common_sequence_patterns","pattern_multiplicities",
          "n_d_kiusiana_accessions","n_d_kiusiana_unique_patterns",
          "shared_pattern_across_three_kiusiana_groups")}),flush=True)
    # These are independently reproduced from the downloaded accession set
    # and will fail rather than silently update earlier quantitative claims.
    if info["common_bp"]!=500 or info["unique_common_sequence_patterns"]!=9:
        raise SystemExit("Alignment changed; review before asserting 500bp/9-pattern result")
if __name__=="__main__":main()
