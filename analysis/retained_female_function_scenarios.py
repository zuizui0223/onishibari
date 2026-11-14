#!/usr/bin/env python3
"""Conditional thought experiment: why retain weak female function?

Inputs: seed survival counts from Shibata (2020) for *Daphne jezoensis*.
Outputs: hypothetical fitness budget for *D. pseudomezereum*, with NO target
species fit, no true colonization probabilities, no sex-system equilibrium.

Compares retention (R) with an otherwise identical, hypothetical fully male
morph (M) in two states: pollen mates locally available vs unavailable.
Only a stripped-down nuclear diploid genetic transmission proxy is used.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "derived" / "shibata2020_seed_pipeline.csv"
DEFAULT_OUT = ROOT / "docs" / "2026-11-14_retained_pistil_thresholds.md"

def source_viability(path=SRC):
    with path.open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    data={(r["population"],r["pollination_type"]):r for r in rows}
    if len(data)!=8:
        raise ValueError("Published source table changed; inspect before running")
    a,b=data[("pooled","outcross")],data[("pooled","self")]
    if {a["taxon"],b["taxon"]}!={"Daphne jezoensis"}:
        raise ValueError("Source is not the intended D. jezoensis study")
    if a["source_url"]!="https://doi.org/10.14943/doctoral.k14152" or b["source_url"]!=a["source_url"]:
        raise ValueError("Thesis source DOI mismatch")
    def as_fraction(row):
        sow=int(row["seeds_sown"])
        germ=int(row["germinated_after_one_year"])
        alive=int(row["survived_three_years_after_germination"])
        if not 0<=alive<=germ<=sow or sow==0:
            raise ValueError("Impossible seed stage counts")
        return Fraction(alive,sow)
    vo,vs=as_fraction(a),as_fraction(b)
    if (vo,vs)!=(Fraction(31,40),Fraction(6,19)):
        raise ValueError("The fixed published seed outcomes were changed")
    return float(vs/vo), float(vo), float(vs)

def transmission_increment(p:float, self_seeds:float, male_cost:float,
                           fixed_cost:float, relative_self_viability:float,
                           gene_copies:float=2.0) -> float:
    """Extra expected transmitted-copy equivalents of R vs fully male M.

    When mates are available (fraction 1-p), retaining female function incurs
    male success loss 'male_cost' in units of viable outcross progeny sired.
    When mates unavailable (fraction p), fully male has zero siring success
    and retaining morph may produce expected 'self_seeds' already *matured*.
    For an isolated self seed the parent transmits two nuclear copies; the
    three-year survival discount is relative to outcross progeny.
    A fixed resource/future fitness cost is paid in both states.

    This excludes all female outcrossing in normal environments, frequency
    dependence, self pollen and seed discounting, genetic sex determination,
    and stochastic population establishment. It cannot establish stability.
    """
    inputs=(p,self_seeds,male_cost,fixed_cost,relative_self_viability,gene_copies)
    if not all(math.isfinite(x) for x in inputs):
        raise ValueError("All inputs must be finite")
    if not 0<=p<=1 or any(x<0 for x in inputs[1:]) or not 0<relative_self_viability<=1:
        raise ValueError("Inputs outside defined scenario range")
    return p*gene_copies*relative_self_viability*self_seeds - (1-p)*male_cost - fixed_cost

def threshold(self_seeds, male_cost, fixed_cost, r, copies=2.0):
    """Smallest mate-scarcity fraction above which R has positive return.

    None means R never has strictly positive advantage in p in [0,1].
    """
    if any(not math.isfinite(x) for x in (self_seeds,male_cost,fixed_cost,r,copies)):
        raise ValueError("Nonfinite parameter")
    if min(self_seeds,male_cost,fixed_cost)<0 or not 0<r<=1 or copies<=0:
        raise ValueError("Invalid parameter")
    if self_seeds==0:
        return None
    denom = copies*r*self_seeds + male_cost
    p=(male_cost+fixed_cost)/denom
    return p if p<1 else None

def check():
    r,vo,vs=source_viability()
    assert abs(r - 240/589)<1e-14
    assert abs(1-r - .593)<.001
    assert 2*r<1  # δ exceeds the 0.5 replacement threshold
    assert threshold(0,0.05,0,r) is None
    assert threshold(0.5,0.05,0,r) < threshold(0.5,0.10,0,r)
    assert threshold(0.5,0.05,0,r) > threshold(1,0.05,0,r)
    t=threshold(.5,.05,0,r)
    assert 0.10 < t < 0.12
    assert abs(transmission_increment(t,.5,.05,0,r))<1e-12
    assert transmission_increment(.05,.5,.05,0,r)<0
    assert transmission_increment(.20,.5,.05,0,r)>0
    assert threshold(.5,.05,100,r) is None
    assert threshold(.5,.05,0,r,1)>threshold(.5,.05,0,r,2)
    print("PISTIL_RETENTION_MODEL_CHECK=SUCCESS",
          f"r={r:.6f} inbreeding_depression={1-r:.6f}",
          f"baseline_break_even_p={t:.6f}",flush=True)

def report():
    r,vo,vs=source_viability()
    s,c,k=.5,.05,0.0
    t=threshold(s,c,k,r)
    td=threshold(s,c,k,r,copies=1)
    lines=[
        "# Conditional maintenance of residual female function: toy threshold model",
        "",
        "**THEORETICAL SCENARIOS ONLY. No *D. pseudomezereum* mating-system data, no estimated colonization frequency, no true sex-system fitness or equilibrium.**",
        "",
        "**REVISION 2026-11-14:** The older 10.9% threshold assumes mature selfed seeds are produced whenever compatible mates are absent; it does NOT model whether effective pollinators can transfer self pollen. This can fail biologically even in self-compatible plants. See [the two-channel model](../analysis/two_channel_pollen_mate_sensitivity.py) for a mate × pollinator sensitivity test; neither model is fitted to Onishibari.",
        "",
        "### Only empirical input: viability after a seed already exists, from a different species",
        "",
        f"- *Daphne jezoensis* Shibata (2020): outcrossed offspring retained {31}/{40} = {vo:.3f} and selfed {6}/{19} = {vs:.3f} through the three-year stage.",
        f"- Selfed survival relative to outcrossed survival **r = {r:.5f}**, implying inbreeding depression **δ = {1-r:.5f}** at this measured horizon.",
        f"- A diploid selfed child receives two parental nuclear gene copies; a child produced via one outcross parental route receives one. The crude relative return of a selfed versus outcrossed seed is **2r = {2*r:.5f}** (below 1 when δ > 0.5).",
        "- These are **post-seed** values, not rates of autonomous selfing, fruiting, pollen transfer, complete lifetime fitness or the target species' parameters.",
        "",
        "### Two-state conditional budget (simplified)",
        "",
        "- p: fraction of reproductive opportunities under severe absence of compatible mates (not measured and NOT equivalent to geographic isolation distance alone).",
        "- S: expected mature selfed seeds produced **per mate-limited opportunity** by a residual-female-function morph (including whether self pollen can be delivered).",
        "- c: loss of viable outcross offspring sired in a mate-rich opportunity, comparing a residual-female morph with an otherwise identical fully male morph.",
        "- k: extra fixed maintenance/future-fitness cost, in the same normalized gene-copy-equivalent units per reproductive opportunity.",
        "- r: post-seed self viability relative to outcross viability, borrowed from *D. jezoensis* solely for illustration.",
        "",
        "For nuclear genetic transmission accounting (gene-copy proxy), an assumed difference between retaining a pistil (R) and full male loss (M) is:",
        "",
        "    ΔG = p × (2 r S) − (1 − p) × c − k",
        "",
        "Female function is predicted to have a positive *direct return in this toy budget* for p above:",
        "",
        "    p* = (c + k)/(2 r S + c), if p* < 1 and S > 0.",
        "",
        "This is NOT a gynodioecy-equilibrium equation and ignores selfing and outcross female output when mates are available, pollen/seed discounting, frequency dependence, nuclear–cytoplasmic sex determination, and variation in paternal mating success. It must not be treated as an adaptation test.",
        "",
        "### Illustrative outcomes, fixed maintenance cost k = 0 (entirely assumed)",
        "",
        "| S: mature self seeds when mates scarce | c: male return lost when mates present | Gene-copy break-even p* | If counting surviving offspring only, p* |",
        "|---:|---:|---:|---:|",
    ]
    for S in (0.25,0.5,1.0):
        for C in (.02,.05,.10):
            tg=threshold(S,C,0,r,copies=2)
            ts=threshold(S,C,0,r,copies=1)
            lines.append(f"| {S:.2f} | {C:.2f} | {100*tg:.1f}% | {100*ts:.1f}% |")
    lines += [
        "",
        f"**Example**: if S={s}, c={c} and k=0, a residual-female morph has positive genetic-transmission return only when severe mate scarcity occurs in more than **{t*100:.1f}%** of opportunities. If counting only surviving descendants, the analogous threshold is **{td*100:.1f}%**. These are distinct units and NOT empirical rates.",
        f"If fixed cost k=0.02 in that same gene-copy scenario, p* becomes **{threshold(s,c,.02,r)*100:.1f}%**.",
        "- If S=0 (no viable selfed mature seeds when mates absent), retaining this function never yields a strict positive direct benefit in this minimal model.",
        "- If c=k=0 and S>0, any p>0 yields a positive direct return but that does not demonstrate an evolutionary constraint or persistence.",
        "",
        "### What to measure in *D. pseudomezereum*",
        "",
        "1. Autonomy: bagged flower self-seed formation and viability in the absence of visitors; if there is no autonomous selfing, mate absence is not necessarily a selfing opportunity.",
        "2. Viability: genetic self versus outcross assignment, mature seeds per manipulated flower, three-year and preferably lifetime survival; estimate r with maternal/genet-level uncertainty.",
        "3. Costs: repeated whole-plant flower production, pollen viability, actual genetic siring and next-year performance, with reproductive manipulation and a sham; histological pistil formation costs are not the same as fruit production costs.",
        "4. Frequency of mate scarcity: mapped, synchronous flowering neighbors, actual compatible pollen deposition, pollinator activity, and colonization history; ordinary geographic distance is not a direct p measurement.",
        "5. Sex-system genetics and frequency-dependent selection require a separate, full population-genetic model after phenotypic and genetic data exist.",
        "",
        "### References",
        "- Shibata (2020) PhD thesis, data source: https://doi.org/10.14943/doctoral.k14152",
        "- Shibata et al. (2018) near-functional dioecy and no assured population selfing: https://doi.org/10.1007/s10265-017-0978-5",
        "- Busch and Delph (2012) reproductive assurance vs automatic transmission advantage: https://pmc.ncbi.nlm.nih.gov/articles/PMC3278291/",
        "- Broader theoretical review of mating-system transitions: https://pmc.ncbi.nlm.nih.gov/articles/PMC3278288/",
        "",
    ]
    return "\n".join(lines)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--output",type=Path,default=DEFAULT_OUT)
    args=parser.parse_args()
    check()
    if not args.check:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(report(),encoding="utf-8")
        print("PISTIL_RETENTION_REPORT="+str(args.output),flush=True)
