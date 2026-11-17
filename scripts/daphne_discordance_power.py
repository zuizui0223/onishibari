#!/usr/bin/env python3
"""Quantitative exploratory locus planning for focal Daphne P/C/K phylogenomics.

Under the *simple*, panmictic, neutral three-species multispecies coalescent
WITHOUT introgression, each discordant quartet topology has probability
exp(-t)/3 (t in coalescent units). The number of discordant loci is stochastic
and estimated gene trees are NOT necessarily true genealogies.

A conditional test compares two discordant topology counts under H0:p=.5
versus illustrative H1:p=2/3 (a 2:1 imbalance). This does not identify gene
flow; ancestor structure and gene-tree estimation error can produce asymmetry.
The output is a sensitivity/feasibility guide, NOT a sample-size guarantee.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

OUT=Path("analysis_results_nuclear_power")
OUT.mkdir(exist_ok=True)

def binom_cdf_05(k,n):
    if k<0:return 0.0
    if k>=n:return 1.0
    return sum(math.comb(n,i) for i in range(k+1))/2**n

def two_sided_exact_pvalue(k,n):
    # Under p=1/2 a probability-ordering two-sided test is just twice
    # the lower-tail probability for k<=n/2, capped at 1.
    lo=min(k,n-k)
    return min(1.0, 2*binom_cdf_05(lo,n))

def alt_power(n,q=2/3,alpha=.05):
    if not 0<q<1 or n<1:raise ValueError("invalid")
    # Compute pmfs recurrently to avoid very large factorials; n <=250 here.
    pmf=(1-q)**n
    probability=0.0
    for k in range(n+1):
        if two_sided_exact_pvalue(k,n)<=alpha:
            probability+=pmf
        if k<n:
            pmf=pmf*(n-k)/(k+1)*q/(1-q)
    return probability

def main():
    tuples=[]
    for nd in [20,40,60,80,100,150,200]:
        tuples.append({"discordant_gene_trees":nd,
                       "power_detect_2to1_imbalance":round(alt_power(nd),5),
                       "null_exact_two_sided_alpha":.05})
    t_rows=[]
    for t in [.25,.5,1.0,1.5]:
        each_discordant=math.exp(-t)/3
        total_discordant=2*each_discordant
        target=80
        t_rows.append({"internal_branch_coalescent_units":t,
           "P_concordant_species_tree":round(1-2*math.exp(-t)/3,5),
           "P_each_discordant":round(each_discordant,5),
           "P_any_discordant":round(total_discordant,5),
           "expected_loci_to_obtain_80_discordant_gene_trees":
               math.ceil(target/total_discordant),
           "note":"expected count only; never substitute for sequencing power"})
    with (OUT/"asymmetry_power_conditional.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(tuples[0]));w.writeheader();w.writerows(tuples)
    with (OUT/"ils_expected_discordance_by_branch.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(t_rows[0]));w.writeheader();w.writerows(t_rows)
    report={"question":"Is an arbitrary 100-locus pilot sufficient to distinguish ILS from a 2:1 discordant topology skew?",
      "null":"simple 3-taxon neutral MSC without gene flow, single unstructured ancestors, true gene trees",
      "alternative":"illustrative 2:1 imbalance among the TWO discordant gene-tree types (not a biological prediction)",
      "statistical_test":"conditional two-sided exact binomial test on discordant independent locus trees",
      "conditional_power":tuples,
      "branch_sensitivity":t_rows,
      "interpretation":[
       "Power ~0.82 requires approximately 80 OBSERVED discordant independent genealogies for a 2:1 imbalance; 100 total loci may be inadequate.",
       "At t=0.5, roughly 198 total independent loci would yield 80 discordants on average. At t=1.0, roughly 327 loci would be needed on average.",
       "These are illustrative expectations and have no allowance for gene-tree error, linkage, missingness, unequal lineage sampling or drift heterogeneity.",
       "A significant skew does NOT uniquely implicate gene flow; ancestral structure, ghost taxa and ascertainment bias remain.",
       "A symmetric gene-tree distribution cannot rule out organelle capture, which may leave no nuclear trace.",
       "First collect a 50-locus pilot to measure recoverability and estimate coalescent internal branch and gene-tree reliability; adapt total target after that."]
    }
    (OUT/"power_design_sensitivity.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
