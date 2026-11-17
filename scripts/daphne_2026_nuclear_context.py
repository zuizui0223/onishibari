#!/usr/bin/env python3
"""Gene-tree quartet audit in the Bonnaudet et al. (2026) 89-tip Angiosperms353 data.

Scope: only Daphne sureil, tangutica, mucronata, laureola; **not** focal East
Asian Daphne pseudomezereum/kamtschatica/koreana. Gene-tree quartet frequency
here is a contextual genus-wide benchmark and cannot adjudicate P/C/K ILS vs
introgression. No invented localities or sampling.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import io
import json
import re
import tarfile
import time
from pathlib import Path
from urllib.request import urlopen,Request
from Bio import Phylo

OUT=Path("analysis_results_daphne_2026_nuclear")
OUT.mkdir(exist_ok=True)
API="https://zenodo.org/api/records/17517526"
HEAD={"User-Agent":"Daphne-Angiosperms353-context/0.4 academic research"}
FOCAL={"D_sureil":"27665_L001","D_tangutica":"27673_L001",
       "D_mucronata":"27691_L001","D_laureola":"07178_L001"}
MAX_BYTES=20*1024*1024

def fetch(url,limit=MAX_BYTES):
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers=HEAD),timeout=120) as r:
                b=r.read(limit+1)
                if len(b)>limit:raise ValueError("File exceeds expected ceiling: "+url)
                return b
        except Exception:
            if attempt==3:raise
            time.sleep(2**attempt)

def splits(t,interest):
    """Extract unrooted quartet bipartition from any (possibly rooted) tree."""
    interest=set(interest)
    present=interest & {x.name for x in t.get_terminals() if x.name}
    if len(present)<4:return None
    parts=set()
    for clade in t.find_clades(order="postorder"):
        leafset={x.name for x in clade.get_terminals() if x.name} & present
        if len(leafset)==2:
            complement=present-leafset
            if len(complement)==2:
                a=tuple(sorted(leafset))
                b=tuple(sorted(complement))
                parts.add(tuple(sorted((a,b))))
    if len(parts)==1:return next(iter(parts))
    return "UNRESOLVED" if not parts else "CONFLICTING"

def parse_newick(x):
    return Phylo.read(io.StringIO(x.decode("utf8","replace").strip()),"newick")

def main():
    report={"dataset_doi":"10.5281/zenodo.17517526",
      "study_doi":"10.1093/aob/mcag206",
      "method":"Unrooted induced 4-tip quartet frequency from Angiosperms353 per-gene trees",
      "scope":"Daphne sureil, tangutica, mucronata, laureola only; P/C/K entirely absent",
      "four_daphne_sample_ids":FOCAL, "files":[], "quartet_counts":{},
      "gene_trees_total":0,"gene_trees_4taxa_covered":0,
      "gene_trees_fewer_4taxa":0,"gene_trees_parse_error":0,
      "tree_examples":[]}
    try:
        record=json.loads(fetch(API,5_000_000))
        files=record.get("files",[])
        by_key={x.get("key",""):x for x in files}
        req=["4_IQTree.tar.gz","5_iqtree-BS10_sp_V2.tre","6_ThymelaeaceaeDated.tre"]
        for name in req:
            if name not in by_key:raise ValueError("Expected Zenodo file missing "+name)
        def obtain(key):
            f=by_key[key]
            url=(f.get("links") or {}).get("self")
            if not url:raise ValueError("Missing direct link for "+key)
            b=fetch(url,MAX_BYTES)
            if len(b)!=f.get("size"):raise ValueError("Size mismatch "+key)
            report["files"].append({"key":key,"size":len(b),"sha256":hashlib.sha256(b).hexdigest()})
            return b
        astral=parse_newick(obtain("5_iqtree-BS10_sp_V2.tre"))
        dated=Phylo.read(io.StringIO(obtain("6_ThymelaeaceaeDated.tre").decode()),"nexus")
        daphne_taxa=[t.name for t in dated.get_terminals() if t.name and
            "_Daphne_" in t.name]
        report["dated_tree_tip_count"]=len(dated.get_terminals())
        report["dated_tree_Daphne_taxa"]=daphne_taxa
        astral_ids={x.name for x in astral.get_terminals()}
        if not all(x in astral_ids for x in FOCAL.values()):
            raise ValueError("Four published Daphne sample IDs absent from ASTRAL")
        report["astral_four_daphne_are_exclusive_monophyletic"] = (
            len(astral.common_ancestor(*FOCAL.values()).get_terminals())==4)
        report["astral_quartet_split"]=splits(astral,FOCAL.values())
        binary=obtain("4_IQTree.tar.gz")
        outcounts=Counter()
        per_gene=[]
        with tarfile.open(fileobj=io.BytesIO(binary),mode="r:gz") as arc:
            members=[m for m in arc.getmembers() if m.isfile()]
            report["archive_members_total"]=len(members)
            report["archive_sample_names"]=[m.name for m in members[:12]]
            # Treefile paths are nested; use extension and parse confirmation.
            for member in members:
                if member.size>2_000_000:continue
                if not re.search(r"(treefile|\\.tre$|\\.tree$|\\.nwk$)",member.name,re.I):
                    continue
                f=arc.extractfile(member)
                if not f:continue
                raw=f.read()
                try:
                    t=parse_newick(raw)
                    ntips=len(t.get_terminals())
                    if ntips<4:continue
                    report["gene_trees_total"]+=1
                    part=splits(t,FOCAL.values())
                    if part is None:
                        report["gene_trees_fewer_4taxa"]+=1
                        per_gene.append([member.name,ntips,0,"missing_1plus_Daphne"])
                    else:
                        report["gene_trees_4taxa_covered"]+=1
                        key=part if isinstance(part,str) else "/".join("-".join(x) for x in part)
                        outcounts[key]+=1
                        per_gene.append([member.name,ntips,1,key])
                        if len(report["tree_examples"])<8:
                            report["tree_examples"].append({"path":member.name,"quartet":key})
                except Exception as e:
                    report["gene_trees_parse_error"]+=1
                    if len(report["tree_examples"])<8:
                        report["tree_examples"].append({"path":member.name,"error":str(e)})
        report["quartet_counts"]=dict(outcounts)
        if len(per_gene)!=report["gene_trees_total"]:
            raise ValueError("Per-gene table and number of parsed gene trees disagree")
        (OUT/"per_gene_coverage_and_quartet.tsv").write_text(
            "gene_tree_path\\tn_gene_tree_tips\\tall_four_daphne_present\\tunrooted_quartet\\n".replace("\\t","\t").replace("\\n","\n")
            + "".join("\t".join(map(str,row))+"\n" for row in per_gene),
            encoding="utf8")
        denominator=report["gene_trees_4taxa_covered"]
        report["quartet_fractions"]={
            key:round(n/denominator,5) for key,n in outcounts.items()} if denominator else {}
        report["completed"]=True
        if report["gene_trees_total"]==0:
            report["warning"]="No gene-tree files parsed; inspect archive member paths and format"
    except Exception as e:
        report["completed"]=False
        report["error"]=type(e).__name__+": "+str(e)
    (OUT/"nuclear_context_report.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,default=str)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2,default=str))
    if not report["completed"]:raise RuntimeError(report["error"])
if __name__=="__main__":main()
