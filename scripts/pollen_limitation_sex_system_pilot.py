#!/usr/bin/env python3
"""Exploratory, taxonomically exact descriptive link of GloPL with Wang sexual systems.

Pollen limitation (PL) as published standardized log response ratio, not
visitor frequency. One row = experimental case; first condense to species
medians to avoid treating repeated site-years as independent species.
No phylogeny, no evolutionary inference, no p-value; trend directions are
hypothesis-generating only. Treat whole-plant vs floral supplement separately.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import csv, hashlib, json, math, statistics, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from reproductive_data_coverage import get,find_glopl,csvrows,clean_name

OUT=Path("analysis_results_pl_pilot")
OUT.mkdir(exist_ok=True)
WANG_META="https://zenodo.org/api/records/4895147"
SELECT="SplistWithSexSyst.csv"

def convert(x):
    try:
        v=float(x)
        if math.isfinite(v):return v
    except Exception:pass
    return None

def summarize_group(values):
    pairs=defaultdict(list)
    for n,v in values:pairs[n].append(v)
    med=[statistics.median(v) for v in pairs.values()]
    return {
        "species":len(med),"cases":sum(map(len,pairs.values())),
        "median_of_species_medians":round(statistics.median(med),5) if med else None,
        "mean_of_species_medians":round(statistics.mean(med),5) if med else None,
        "fraction_species_median_positive":round(sum(v>0 for v in med)/len(med),4) if med else None,
        "min_species":round(min(med),5) if med else None,
        "max_species":round(max(med),5) if med else None,
    }

def main():
    d={"method":"Exact-name Wang2021 sexual system joined to GloPL standardized PL_Effect_Size; descriptive species medians",
       "pre_analysis_warnings":[
           "Wang deposited 72,459 names, paper analysis used 68,453; filters not replicated",
           "GloPL ES may be defined via supplemental or bagged cross-pollen and different response types",
           "GloPL has publication/location/year/morph repeated cases, no phylogenetic correction",
           "This comparison does not measure visitation rarity, mating success of both sexes, or speciation",
           "Whole plant supplement is necessary for direct lifetime fitness extrapolation; even whole-plant cases are not necessarily lifetime observations",
           "Cross-sectional comparison only: no causal estimates or selection inference"]}
    try:
        metadata=json.loads(get(WANG_META,3_000_000)[0])
        f=next(x for x in metadata["files"] if x["key"]==SELECT)
        wb,wu=get(f["links"]["self"],6_000_000)
        gb,gm=find_glopl()
        W_rows=csvrows(wb)[1]
        G_rows=csvrows(gb)[1]
        cols=list(G_rows[0])
        required=["PL_Effect_Size","Species_accepted_names","PL_Effect_Size_Type1",
                  "PL_Effect_Size_Type2","Level_of_Supplementation",
                  "Included_Bennett_et_al_meta_analysis"]
        if not set(required)<=set(cols):raise ValueError("Missing final GloPL columns "+str(set(required)-set(cols)))
        W={}
        for r in W_rows:
            n=clean_name(r.get("Genus_species",""))
            if not n:continue
            s=r.get("SexSyst","")
            if n in W and W[n]!=s:raise ValueError("Conflicting Wang class for "+n)
            W[n]=s
        accepted={"hermaphroditism","dioecy","monoecy"}
        baskets=defaultdict(list)
        types=Counter();response=Counter();levels=Counter()
        data_rows=[]
        for row in G_rows:
            if row.get("Included_Bennett_et_al_meta_analysis","").lower()!="yes":continue
            n=clean_name(row.get("Species_accepted_names",""))
            value=convert(row.get("PL_Effect_Size"))
            if not n or value is None:continue
            cat=W.get(n)
            lev=row.get("Level_of_Supplementation","NA")
            response[row.get("PL_Effect_Size_Type1","NA")]+=1
            types[row.get("PL_Effect_Size_Type2","NA")]+=1
            levels[lev]+=1
            if cat not in accepted:continue
            typ=row.get("PL_Effect_Size_Type2","NA")
            data_rows.append({
               "species":n,"sex_system":cat,"effect_size_log_ratio":value,
               "level":lev,"outcome_type":row.get("PL_Effect_Size_Type1","NA"),
               "treatment_type":typ,
               "source_author":row.get("Author",""),
               "source_year":row.get("Year","")})
            for label in ("all", lev, "treatment:"+typ):
                baskets[(label,cat)].append((n,value))
        d["source_files"]={"Wang":dict(url=wu,sha256=hashlib.sha256(wb).hexdigest()),
                            "GloPL":dict(url=gm["download"],sha256=hashlib.sha256(gb).hexdigest(),release=gm["release"])}
        d["source_GloPL_records"]=len(G_rows)
        d["included_numeric_GloPL_records"]=sum(levels.values())
        d["included_numeric_case_response_types"]=dict(response)
        d["included_numeric_case_treatment_types"]=dict(types)
        d["included_numeric_case_supplementation_levels"]=dict(levels)
        d["matched_numeric_cases"]=len(data_rows)
        d["matched_numeric_species"]=len({r["species"] for r in data_rows})
        groups={}
        for k,v in sorted(baskets.items()):
            label,sex=k
            groups.setdefault(label,{})[sex]=summarize_group(v)
        d["descriptive_groups"]=groups
        cols=["species","sex_system","effect_size_log_ratio","level","outcome_type","treatment_type","source_author","source_year"]
        with (OUT/"linked_PLeffect_cases_source_graded.csv").open("w",newline="",encoding="utf8") as f:
            writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader();writer.writerows(data_rows)
        d["status"]="success"
    except Exception as e:
        d["status"]="failure";d["error"]=type(e).__name__+": "+str(e)
    (OUT/"descriptive_summary.json").write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(d,ensure_ascii=False,indent=2)[:50000])
    if d["status"]!="success":raise RuntimeError(d["error"])
if __name__=="__main__":main()
