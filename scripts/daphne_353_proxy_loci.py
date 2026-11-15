#!/usr/bin/env python3
"""Identify reusable Daphne Angiosperms353 alignment loci from Bonnaudet (2026).

NOT the focal P/C/K sequences. Outputs a candidate *proxy* reference panel
for testing low-copy nuclear read recovery when original P/C/K FASTQ exists.
An orthology/paralogy screen remains essential before formal phylogenomics.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
from pathlib import Path
import statistics
import tarfile
import time
from urllib.request import Request,urlopen
from Bio import SeqIO

OUT=Path("analysis_results_proxy_loci")
OUT.mkdir(exist_ok=True)
RECORD="https://zenodo.org/api/records/17517526"
HEAD={"User-Agent":"Daphne-353-proxy-locus-inspection/0.4 (academic)"}
DAPHNE={"sureil":"27665_L001","tangutica":"27673_L001",
        "mucronata":"27691_L001","laureola":"07178_L001"}
NAMES=set(DAPHNE.values())
ACGT=set("ACGT")

def fetch(url, limit=15_000_000):
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers=HEAD),timeout=120) as h:
                b=h.read(limit+1)
                if len(b)>limit:raise ValueError("Exceeds download limit "+url)
                return b
        except OSError:
            if attempt==3:raise
            time.sleep(2**attempt)

def fasta_records(data):
    try:
        txt=data.decode("utf-8")
        if not txt.lstrip().startswith(">"):return {}
        records=list(SeqIO.parse(io.StringIO(txt),"fasta"))
        return {r.id:str(r.seq).upper().replace("U","T") for r in records if r.id}
    except (UnicodeError,ValueError):
        return {}

def extract_four(recs):
    result={}
    for name,ID in DAPHNE.items():
        matches=[(hdr,seq) for hdr,seq in recs.items() if hdr==ID
                 or hdr.startswith(ID+"|") or hdr.startswith(ID+"_")]
        if len(matches)==1:result[name]=matches[0][1]
        if len(matches)>1:raise ValueError(f"ambiguous multiple taxa for {ID}")
    return result

def main():
    report={"source":"Bonnaudet et al. (2026) Angiosperms353 dataset 10.5281/zenodo.17517526",
      "focal_P_C_K_in_reference":False,"proxy_taxa":DAPHNE,
      "files_examined":0,"valid_fasta_loci":0,"loci_all_four":0,
      "unexpected_formats":[],"provenance":{},"complete":False}
    rows=[]; best=[]
    try:
        meta=json.loads(fetch(RECORD,5_000_000))
        source=next(f for f in meta["files"] if f.get("key")=="3_AlignedLoci.tar.gz")
        raw=fetch(source["links"]["self"],15_000_000)
        if len(raw)!=source.get("size"):raise ValueError("Zenodo source file size mismatch")
        report["provenance"]={"filename":source["key"],"size":len(raw),
           "sha256":hashlib.sha256(raw).hexdigest(),"archive_items":0}
        with tarfile.open(fileobj=io.BytesIO(raw),mode="r:gz") as tar:
            files=[m for m in tar.getmembers() if m.isfile()]
            report["provenance"]["archive_items"]=len(files)
            for m in files:
                report["files_examined"]+=1
                if m.size>400_000:continue
                f=tar.extractfile(m)
                if f is None:continue
                dat=f.read()
                recs=fasta_records(dat)
                if not recs:
                    if len(report["unexpected_formats"])<12:
                        report["unexpected_formats"].append({"name":m.name,
                           "first_chars":dat[:100].decode("utf8","replace")})
                    continue
                report["valid_fasta_loci"]+=1
                d=extract_four(recs)
                if len(d)<4:continue
                report["loci_all_four"]+=1
                lens={len(x) for x in d.values()}
                if len(lens)!=1:continue
                alnlen=lens.pop()
                complete=[i for i in range(alnlen) if all(s[i] in ACGT for s in d.values())]
                complete_n=len(complete)
                var=sum(len({s[i] for s in d.values()})>1 for i in complete)
                ungapped=[sum(ch in ACGT for ch in seq) for seq in d.values()]
                item=dict(source_member=m.name,locus_id=Path(m.name).stem,
                  aligned_sites=alnlen,all_four_ACGT_sites=complete_n,
                  complete_fraction=round(complete_n/alnlen,5) if alnlen else 0,
                  variable_complete_sites=var,
                  min_Daphne_ungapped_bases=min(ungapped),
                  max_Daphne_ungapped_bases=max(ungapped),
                  proxy_taxa_count=4)
                rows.append(item)
                if complete_n>=200 and min(ungapped)>=250 and
                    0.001<=var/max(complete_n,1)<=0.15:
                    best.append((item,d))
        rows.sort(key=lambda x:x["locus_id"])
        report["qualifying_proxy_loci"]=len(best)
        if report["loci_all_four"]==0:
            raise ValueError("No four-taxon alignment rows; inspect archive formats or sample identifiers")
        if len(best)<10:
            raise ValueError(f"Only {len(best)} qualifying loci; no 50-locus proxy panel built")
        fields=["source_member","locus_id","aligned_sites","all_four_ACGT_sites",
                "complete_fraction","variable_complete_sites",
                "min_Daphne_ungapped_bases","max_Daphne_ungapped_bases","proxy_taxa_count"]
        with (OUT/"all_four_daphne_alignment_coverage.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        # Maximize aligned coverage first, vary loci by variable-site rate.
        best.sort(key=lambda x:(x[0]["all_four_ACGT_sites"],
            x[0]["variable_complete_sites"],x[0]["min_Daphne_ungapped_bases"]),reverse=True)
        top=best[:50]
        report["selected_loci"]=len(top)
        report["selected_four_daphne_complete_sites_min"]=min(x[0]["all_four_ACGT_sites"] for x in top)
        report["selected_four_daphne_complete_sites_median"]=statistics.median(
            x[0]["all_four_ACGT_sites"] for x in top)
        with (OUT/"candidate_proxy_loci_top50.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=fields)
            writer.writeheader()
            writer.writerows([x[0] for x in top])
        with (OUT/"candidate_proxy_fourDaphne_sequences.fa").open("w") as fh:
            for item,d in top:
                for species,aligned in d.items():
                    # The ungapped DNA fragments are proxy reference sequences.
                    seq="".join(ch for ch in aligned if ch in ACGT)
                    fh.write(f">{item['locus_id']}|Daphne_{species}|{DAPHNE[species]}\n")
                    for i in range(0,len(seq),80):
                        fh.write(seq[i:i+80]+"\n")
        report["complete"]=True
        report["limitations"]=[
          "All four Daphne are distant from the Japanese/Korean/Kamchatka P/C/K trio",
          "Reference panel does not infer P/C/K genealogy and does not certify orthology",
          "Original published alignment may include multiple high-copy sequences or artifacts",
          "Loci selected by completeness and variability in four other Daphne taxa",
          "No reference genome is required for an Angiosperms353 capture pilot"]
    except Exception as exc:
        report["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"proxy_panel_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False))
    if not report["complete"]:raise RuntimeError(report["error"])
if __name__=="__main__":main()
