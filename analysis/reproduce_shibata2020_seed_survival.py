#!/usr/bin/env python3
"""Reproduce reported Daphne jezoensis seed-stage fitness from Shibata (2020).

Data: original dissertation chapter 3, pooled seeds from 3 populations.
Never treat this analysis as observations of D. pseudomezereum.
No plant/genet-level SEs can be inferred from these pooled counts.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/derived/shibata2020_seed_pipeline.csv"
REPORT = ROOT / "docs/shibata2020_seed_pipeline_results.md"
URL = "https://doi.org/10.14943/doctoral.k14152"

def exact_fisher(a, b, c, d):
    """Conditional two-sided Fisher exact p-value, ties included."""
    n1, n2 = a+b, c+d
    successes = a+c
    total = n1+n2
    lo = max(0, successes-n2)
    hi = min(n1, successes)
    denominator = comb(total, n1)
    def hypergeom(x):
        return Fraction(comb(successes, x)*comb(total-successes, n1-x),
                        denominator)
    observed = hypergeom(a)
    return float(sum((hypergeom(x) for x in range(lo,hi+1)
                      if hypergeom(x) <= observed), Fraction(0)))

def load(path=DATA):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 8:
        raise ValueError("Expected 6 site composition rows and 2 pooled outcomes")
    for r in rows:
        if r["taxon"] != "Daphne jezoensis" or r["source_url"] != URL:
            raise ValueError("Source/taxon mismatch: do not mix Onishibari and Naniwazu")
    key = {(r["population"],r["pollination_type"]):r for r in rows}
    if len(key) != 8:
        raise ValueError("Duplicate or missing unique populations/treatments")
    for type_ in ("self","outcross"):
        pooled=key[("pooled",type_)]
        n=int(pooled["seeds_sown"])
        by_sites = sum(int(key[(site,type_)]["seeds_sown"])
                       for site in ("CTS","NOP","TOM"))
        if n != by_sites:
            raise ValueError("Pooled seeds differ from site sum")
        germ, survivors = (int(pooled[field]) for field in
                           ("germinated_after_one_year",
                            "survived_three_years_after_germination"))
        if not 0 <= survivors <= germ <= n:
            raise ValueError("Biologically impossible fate table")
    return key

def calculate(key):
    def n(t, fld):
        return int(key[("pooled",t)][fld])
    xo, xs = (n(t,"seeds_sown") for t in ("outcross","self"))
    go, gs = (n(t,"germinated_after_one_year") for t in ("outcross","self"))
    so, ss = (n(t,"survived_three_years_after_germination")
              for t in ("outcross","self"))
    rates = {
        "outcross_n":xo, "self_n":xs,
        "outcross_germinated_n":go, "self_germinated_n":gs,
        "outcross_survivors_n":so, "self_survivors_n":ss,
        "germ_outcross":go/xo,
        "germ_self":gs/xs,
        "survival_cond_outcross":so/go,
        "survival_cond_self":ss/gs,
        "cumulative_outcross":so/xo,
        "cumulative_self":ss/xs,
        "viability_ratio":(ss/xs)/(so/xo),
        "inbreeding_depression":1-(ss/xs)/(so/xo),
        "fisher_germination":exact_fisher(go,xo-go,gs,xs-gs),
        "fisher_conditional_survival":exact_fisher(so,go-so,ss,gs-ss),
        "fisher_cumulative":exact_fisher(so,xo-so,ss,xs-ss),
    }
    # Reconstructed counts reproduce rounded published percentages.
    assert abs(rates["germ_outcross"] - .95)<1e-12
    assert abs(rates["germ_self"] - .474)<.0005
    assert abs(rates["survival_cond_outcross"]-.816)<.0005
    assert abs(rates["survival_cond_self"]-.667)<.0005
    assert abs(rates["inbreeding_depression"]-.59)<.005
    assert abs(rates["fisher_conditional_survival"]-.38)<.02
    return rates

def report(x):
    lines=[
        "# Daphne jezoensis: seed-to-three-year survival pipeline",
        "",
        "**External research reanalysis, NOT an Onishibari (D. pseudomezereum) reproductive experiment.** Source: [Shibata (2020), doctoral thesis, chapter 3]("+URL+").",
        "",
        "The dissertation reports 59 hand-pollinated seeds from CTS, NOP and TOM; 40 were outcrossed and 19 selfed. Percentages reconstruct the following pooled integer outcomes:",
        "",
        "| Stage | Hand outcross | Hand self |",
        "|---|---:|---:|",
        f"| Seeds sown | {x['outcross_n']} | {x['self_n']} |",
        f"| Germinated in year 1 | {x['outcross_germinated_n']} ({x['germ_outcross']*100:.1f}%) | {x['self_germinated_n']} ({x['germ_self']*100:.1f}%) |",
        f"| Survived three years, conditional on germination | {x['outcross_survivors_n']} ({x['survival_cond_outcross']*100:.1f}%) | {x['self_survivors_n']} ({x['survival_cond_self']*100:.1f}%) |",
        f"| Reached end of observation per original seed | {x['outcross_survivors_n']}/{x['outcross_n']} ({x['cumulative_outcross']*100:.1f}%) | {x['self_survivors_n']}/{x['self_n']} ({x['cumulative_self']*100:.1f}%) |",
        "",
        f"**Reproduced cumulative inbreeding depression**: delta = 1 - ({x['self_survivors_n']}/{x['self_n']}) / ({x['outcross_survivors_n']}/{x['outcross_n']}) = **{x['inbreeding_depression']:.3f}** (thesis rounded to 0.59).",
        f"Selfed progeny reached the endpoint at **{x['viability_ratio']*100:.1f}%** the per-seed frequency of outcrossed progeny.",
        "",
        "### Sanity-check of contingency tables (NOT new experimental inference)",
        f"- Two-sided Fisher exact for germination: {x['fisher_germination']:.6g} (thesis: p < 0.001).",
        f"- Two-sided Fisher exact for conditional three-year survival: {x['fisher_conditional_survival']:.3g} (thesis: p = 0.38).",
        f"- Two-sided Fisher exact for cumulative endpoint: {x['fisher_cumulative']:.6g}; exploratory, not reported as primary in the source.",
        "",
        "### Biological and statistical boundaries",
        "- **Fruits/seeds formed after deliberate self-pollination**, hence this nearby species is not completely self-incompatible under those treatments.",
        "- These rates apply **after a seed was produced**; they do not measure success per hand-pollinated flower, pollen-tube compatibility or pre-seed embryo loss.",
        "- The trial pooled 3 source populations because of small individual counts. Seed sibships and maternal identities are unavailable; exact p-values that treat seeds as independent may be anti-conservative.",
        "- Three years of follow-up is not lifetime fitness; further reproduction, male fitness, environmental context and post-transplant survival remain unknown.",
        "- One published cohort cannot estimate Onishibari's seed fertility, its mating system, or the presence of floral self-pollen interference.",
        "",
    ]
    return "\n".join(lines)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--output",type=Path,default=REPORT)
    args=parser.parse_args()
    result=calculate(load())
    if args.check:
        print("SHIBATA_SEED_PIPELINE_CHECK=SUCCESS; delta=%.3f"%result["inbreeding_depression"])
        return
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(report(result),encoding="utf-8")
    print("SHIBATA_SEED_PIPELINE_GENERATED="+str(args.output))

if __name__=="__main__":
    main()
