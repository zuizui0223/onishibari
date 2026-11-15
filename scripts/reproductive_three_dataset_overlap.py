#!/usr/bin/env python3
"""Cross-source coverage for angiosperm reproductive decoupling (pilot only).

Wang2021 72,459 records are the deposited species list, not 68,453 species
used in the publication's analyses. The latter includes analysis exclusions.
GloPL_with_id.csv is authors' curated working copy, 2,969 experiment cases;
each species may have multiple site/year/treatment/morph observations.
PROTEUS authors' cleaned data contains 364 species with varied missingness.
This script measures linkage feasibility, not evolutionary mechanisms.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import csv, hashlib, io, json, math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from reproductive_data_coverage import get,find_glopl,csvrows,clean_name,SOURCES

OUT=Path("analysis_results_overlap")
OUT.mkdir(exist_ok=True)
W_URL="https://zenodo.org/api/records/4895147"
P_URL=SOURCES["proteus"]

def fetch_sources():
    meta=json.loads(get(W_URL,limit=4_000_000)[0])
    file=next(x for x in meta["files"] if x.get("key")=="SplistWithSexSyst.csv")
    wang_blob,wang_url=get(file["links"]["self"],limit=6_000_000)
    prot_blob,prot_url=get(P_URL,limit=1_000_000)
    glopl_blob,glopl_meta=find_glopl()
    report={"wang":dict(url=wang_url,size=len(wang_blob),sha256=hashlib.sha256(wang_blob).hexdigest()),
            "proteus":dict(url=prot_url,size=len(prot_blob),sha256=hashlib.sha256(prot_blob).hexdigest()),
            "glopl":dict(url=glopl_meta["download"],release=glopl_meta["release"],
                         size=len(glopl_blob),sha256=hashlib.sha256(glopl_blob).hexdigest())}
    return csvrows(wang_blob)[1],csvrows(prot_blob)[1],csvrows(glopl_blob)[1],report

def canon(row,col):
    value=clean_name(row.get(col,""))
    if not value or value in ("na","none","nan","unknown"):return ""
    return value

def stat_sort(d):
    return dict(sorted(d.items(),key=lambda z:(-z[1],z[0])))

def main():
    report={"status":"not_started", "limitations":[
        "Exact casefolded scientific names only; synonym resolution and taxonomic ranks not applied",
        "Wang 72,459 deposited list includes species not in published 68,453 analysis; exclusions unknown",
        "Wang labels are sex system classes (hermaphroditism, dioecy, monoecy), while PROTEUS SexualSystem uses morphological monomorphism/dimorphism",
        "GloPL rows are unique experiment cases, NOT independent plant species or life-history trajectories",
        "No phylogenetic model / causal or speciation inference performed"] }
    try:
        wang,prot,glopl,source=fetch_sources()
        report["sources"]=source
        report["raw_rows"]={"Wang":len(wang),"PROTEUS":len(prot),"GloPL_cases":len(glopl)}
        W=defaultdict(set)
        for r in wang:
            n=canon(r,"Genus_species")
            if n:W[n].add(r.get("SexSyst","NA"))
        P={}
        for r in prot:
            n=canon(r,"")
            if n:
                if n in P: raise ValueError("Duplicate PROTEUS species "+n)
                P[n]=r
        G=defaultdict(list)
        for r in glopl:
            n=canon(r,"Species_accepted_names")
            if n:G[n].append(r)
        keyW,keyP,keyG=map(set,(W,P,G))
        triple=keyW&keyP&keyG
        wp=keyW&keyP
        wg=keyW&keyG
        pg=keyP&keyG
        report["unique_species"]={"Wang":len(keyW),"PROTEUS":len(keyP),"GloPL":len(keyG)}
        report["exact_species_intersections"]={
            "Wang_PROTEUS":len(wp),
            "Wang_GloPL":len(wg),
            "PROTEUS_GloPL":len(pg),
            "all_three":len(triple)}
        report["experiments_with_exact_wang"]=sum(len(G[n]) for n in wg)
        report["experiments_with_exact_proteus"]=sum(len(G[n]) for n in pg)
        report["experiments_with_exact_all_three"]=sum(len(G[n]) for n in triple)
        report["Wang_sex_system_categories"]=stat_sort(Counter(
            list(v)[0] if len(v)==1 else "MULTIPLE_CONFLICT" for v in W.values()))
        report["Wang_unique_with_multiple_sex_labels"]=sum(len(v)>1 for v in W.values())
        report["GloPL_repeat_case_distribution"]={
            "unique_species_single_case":sum(len(v)==1 for v in G.values()),
            "unique_species_multiple_cases":sum(len(v)>1 for v in G.values()),
            "max_cases_one_species":max(map(len,G.values())) if G else 0,
            "five_most_replicated_species":[[k,len(v)] for k,v in sorted(G.items(),key=lambda z:-len(z[1]))[:5]]}
        report["experiment_meta_categories"]={
            "included_meta_analysis":stat_sort(Counter(r.get("Included_Bennett_et_al_meta_analysis") for r in glopl)),
            "supplementation_level":stat_sort(Counter(r.get("Level_of_Supplementation") for r in glopl)),
            "outcome_FS_means_both_positive":sum(
                val(r,"Natural_X_FS")>0 and val(r,"Sup_X_FS")>0 for r in glopl)}
        report["PROTEUS_trait_counts"]={
            field:{"nonmissing":sum(str(r.get(field,"NA")) not in ("NA","","nan") for r in prot),
                   "categories":stat_sort(Counter(r.get(field,"NA") for r in prot))}
            for field in ("SexualSystem","Mating","Pollination","Lifespan","FlowerSex","Woodiness","FloralReward")}
        report["Wang_PROTEUS_sex_annotations"]={}
        W_match={}
        for n in wp:
            w=set(W[n])
            if len(w)!=1:continue
            k=(next(iter(w)),P[n].get("SexualSystem","NA"),P[n].get("FlowerSex","NA"))
            W_match[k]=W_match.get(k,0)+1
        report["Wang_PROTEUS_sex_annotations"]={
            "|".join(k):v for k,v in sorted(W_match.items(),key=lambda z:-z[1])}
        csvcols=["species","Wang_SexSyst","PROTEUS_SexualSystem","PROTEUS_Mating",
                 "PROTEUS_FlowerSex","PROTEUS_Lifespan","GloPL_case_count",
                 "GloPL_natural_sup_FS_positive_cases"]
        with (OUT/"triple_overlap_raw_axes.csv").open("w",newline="",encoding="utf8") as fh:
            cw=csv.DictWriter(fh,fieldnames=csvcols)
            cw.writeheader()
            for n in sorted(triple):
                r=P[n]; w=sorted(W[n])
                cw.writerow({"species":n, "Wang_SexSyst":";".join(w),
                             "PROTEUS_SexualSystem":r.get("SexualSystem"),
                             "PROTEUS_Mating":r.get("Mating"),
                             "PROTEUS_FlowerSex":r.get("FlowerSex"),
                             "PROTEUS_Lifespan":r.get("Lifespan"),
                             "GloPL_case_count":len(G[n]),
                             "GloPL_natural_sup_FS_positive_cases":sum(
                                 val(x,"Natural_X_FS")>0 and val(x,"Sup_X_FS")>0 for x in G[n])})
        # This is a compact source-index only; no new pollen limitation effect sizes inferred.
        with (OUT/"species_link_coverage.csv").open("w",newline="",encoding="utf8") as fh:
            cw=csv.writer(fh);cw.writerow(["species","Wang","PROTEUS","GloPL","GloPL_case_count"])
            for n in sorted(keyW|keyP|keyG):
                cw.writerow([n,int(n in W),int(n in P),int(n in G),len(G.get(n,[]))])
        report["status"]="success"
    except Exception as ex:
        report["status"]="failure"
        report["error"]=type(ex).__name__+": "+str(ex)
    (OUT/"overlap_summary.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,default=str)+"\n",encoding="utf8")
    print(json.dumps(report,indent=2,ensure_ascii=False,default=str)[:40000])
    if report["status"]!="success":raise RuntimeError(report["error"])

def val(row,key):
    try:
        x=float(row.get(key,"NA"))
        return x if math.isfinite(x) else float("nan")
    except Exception:
        return float("nan")
if __name__=="__main__":main()
