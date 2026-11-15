#!/usr/bin/env python3
"""Test regional nrDNA marker resolution in the original Lee et al. 2022 panel.

The original 780-bp segment spans ITS1/5.8S/ITS2 with flanking rRNA.
Cross-taxon all-common columns can be more restrictive. These records
did not supply the 2023 population microsatellite selfing-rate estimates.
"""
import argparse
from collections import Counter
from itertools import product
from pathlib import Path
import json

GROUP={
 "Kumamoto":[f"MT623{i}.1" for i in range(684,688)],
 "Jeju":[f"MT623{i}.1" for i in range(688,692)],
 "Hunan":[f"MT623{i}.1" for i in range(692,695)]
}

def fasta(path):
    records={}
    key=None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            key=line[1:].split("|")[0].split()[0]
            if key in records:raise ValueError("Duplicate accession "+key)
            records[key]=""
        else:
            if key is None:raise ValueError("Sequence without accession")
            records[key]+=line.strip().upper()
    if not records:raise ValueError("Empty FASTA")
    return records

def measure(data):
    k=[data[acc] for acc in GROUP["Kumamoto"]]
    j=[data[acc] for acc in GROUP["Jeju"]]
    lens={len(s) for s in k+j}
    if len(lens)!=1:raise ValueError("Need same length/aligned reference positions")
    L=lens.pop()
    if not L or any(set(s)-set("ACGT") for s in k+j):
        raise ValueError("Only ACGT bases valid in shared columns")
    variants=[]
    fixed=[]
    for pos in range(L):
        a=Counter(seq[pos] for seq in k)
        b=Counter(seq[pos] for seq in j)
        if len(a|b)>1:
            variants.append({"position_1based":pos+1,"Kumamoto":dict(a),"Jeju":dict(b)})
        if len(a)==len(b)==1 and set(a)!=set(b):
            fixed.append(pos+1)
    between=[sum(a!=b for a,b in zip(x,y)) for x,y in product(k,j)]
    haplotypes=sorted(set(k+j))
    freqs=[{"haplotype":"H"+str(i+1),"Kumamoto":k.count(h),"Jeju":j.count(h)}
           for i,h in enumerate(haplotypes)]
    return {
      "n_kumamoto":len(k),"n_jeju":len(j),"aligned_bases":L,
      "n_variable_sites":len(variants),"variable_sites":variants,
      "n_fixed_diagnostic_sites":len(fixed),"fixed_diagnostic_sites":fixed,
      "n_distinct_haplotypes":len(haplotypes),
      "n_shared_haplotypes":len(set(k)&set(j)),
      "haplotype_counts":freqs,
      "between_min_snp":min(between),"between_mean_snp":sum(between)/len(between),
      "between_max_snp":max(between),
    }

def selftest():
    ids=GROUP["Kumamoto"]+GROUP["Jeju"]
    toy={id_:seq for id_,seq in zip(ids,["AAAA","AAAA","AAAT","AAAT",
                                            "AAAA","AAAT","AAAT","AAAT"])}
    a=measure(toy)
    assert a["aligned_bases"]==4
    assert a["n_variable_sites"]==1
    assert a["n_fixed_diagnostic_sites"]==0
    assert a["n_shared_haplotypes"]==2
    for acc in GROUP["Jeju"]:toy[acc]="TTTT"
    assert measure(toy)["n_fixed_diagnostic_sites"]>=3
    toy[GROUP["Kumamoto"][0]]="AAA"
    try:measure(toy)
    except ValueError:pass
    else:raise AssertionError("Length mismatch not detected")
    print("KIUSIANA_MARKER_RESOLUTION_SELFTEST=SUCCESS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--raw",type=Path)
    p.add_argument("--common24",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    selftest()
    if args.self_test:return
    if not args.raw or not args.output:
        p.error("Requires --raw and --output, optional --common24")
    original=fasta(args.raw)
    if len(original)!=19:
        raise ValueError("Expected 19 Lee et al. 2022 accessions")
    full=measure(original)
    if full["aligned_bases"]!=780:
        raise ValueError("Unexpected 2022 Kumamoto/Jeju original amplicon length")
    data={
      "source":"Lee et al. 2022: https://doi.org/10.1016/j.pld.2021.11.001",
      "original_2022_780bp_nrDNA_region":full,
      "scope_warning":"780bp includes flanking ribosomal RNA as well as ITS1/5.8S/ITS2",
      "2023_warning":"2023 eastern/western selfing values come from different samples and microsatellite markers; cannot assign to any 2022 accession",
      "branch_warning":"Zero fixed variants in this small ITS sample does not prove no geographic genetic differentiation or common mating system",
    }
    if args.common24:
        aligned=fasta(args.common24)
        if len(aligned)!=24:raise ValueError("Expected 24 aligned accession tips")
        subset=measure(aligned)
        if subset["aligned_bases"]>=full["aligned_bases"]:
            raise ValueError("Cross-taxon overlap should have fewer columns")
        data["twentyfour_taxa_common_reference_positions"]=subset
        data["trim_warning"]="Reference-anchored 24-tip site filtering removes columns; do not treat reduced polymorphism as an independent biological result"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("KIUSIANA_MARKER_RESOLUTION="+json.dumps({
      "full_bp":full["aligned_bases"],
      "full_variants":full["n_variable_sites"],
      "shared_haplotypes":full["n_shared_haplotypes"],
      "fixed_differences":full["n_fixed_diagnostic_sites"],
      "common24_bp":data.get("twentyfour_taxa_common_reference_positions",{}).get("aligned_bases"),
      "common24_variants":data.get("twentyfour_taxa_common_reference_positions",{}).get("n_variable_sites"),
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
