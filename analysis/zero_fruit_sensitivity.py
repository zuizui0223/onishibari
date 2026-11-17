#!/usr/bin/env python3
"""Sensitivity of zero fruit observations to plant-level clustering.

Synthetic scenarios only. NEVER interpret as actual Onishibari fruit counts.
Under a beta-binomial thought experiment, each genet's individual fruit
probability varies; flowers from one genet are not independent genets.
"""
import argparse
import math
from pathlib import Path

REPORT = Path(__file__).resolve().parents[1] / "docs" / "zero_fruit_interpretation.md"

def zero_fruit_probability(genets, flowers_per_genet, mean_fruit_prob, icc):
    if not isinstance(genets,int) or not isinstance(flowers_per_genet,int):
        raise TypeError("Independent genet count and flowers must be integers")
    if genets<1 or flowers_per_genet<1:
        raise ValueError("At least one genet and one flower per genet")
    if not (0<=mean_fruit_prob<=1) or not (0<=icc<1):
        raise ValueError("Probabilities or ICC outside valid range")
    if mean_fruit_prob == 0:
        return 1.0
    if mean_fruit_prob == 1:
        return 0.0
    if icc == 0:
        return (1-mean_fruit_prob)**(genets*flowers_per_genet)
    concentration = 1/icc - 1
    beta = (1-mean_fruit_prob)*concentration
    # E[(1 - p_genet) ^ n_flower] for p ~ Beta(alpha,beta)
    p_no_fruit_in_one_genet = math.prod(
        (beta + j)/(concentration + j) for j in range(flowers_per_genet)
    )
    return p_no_fruit_in_one_genet**genets

def upper_limit_zero_independent_genets(n_genets, alpha=.05):
    """One-sided 1-alpha exact binomial upper bound if 0 of n genets succeed.

    It is NOT valid for n flowers sharing plants. Multiple flowers/genet can
    only use this if the endpoint is 'any success per genet' for independent
    genets and comparable exposure per genet.
    """
    if not isinstance(n_genets,int) or n_genets<=0:
        raise ValueError("Need positive number of independent genets")
    if not 0<alpha<1:
        raise ValueError("alpha must be 0..1")
    return 1-alpha**(1/n_genets)

def check():
    assert abs(upper_limit_zero_independent_genets(20)-0.139108)<1e-5
    assert abs(upper_limit_zero_independent_genets(60)-0.048703)<1e-5
    for plants in (5,10,20):
        for flowers in (1,2,5):
            for probability in (.01,.05,.20):
                independent=zero_fruit_probability(plants,flowers,probability,0)
                assert abs(independent-(1-probability)**(plants*flowers))<1e-12
                if flowers==1:
                    grouped=zero_fruit_probability(plants,flowers,probability,.5)
                    assert abs(grouped-independent)<1e-12
                else:
                    grouped=zero_fruit_probability(plants,flowers,probability,.2)
                    assert grouped>=independent-1e-12
    assert zero_fruit_probability(10,5,0,.2)==1.0
    assert zero_fruit_probability(10,5,1,.2)==0.0
    print("ZERO_FRUIT_CLUSTER_CHECK=SUCCESS; independent_genet_upper20="
          f"{upper_limit_zero_independent_genets(20):.4f}",flush=True)

def report():
    out=[
        "# Why zero observed fruit does not prove female-function loss",
        "",
        "**This entire table is hypothetical. No Onishibari (D. pseudomezereum) fruit observations were collected or fitted.**",
        "",
        "A negative result is interpretable only after checking pollen viability, a successful hand-outcross female positive control, bag leakage, sample size, independent genet identity and embryo/seed endpoints.",
        "",
        "## No successes among independent plants",
        "",
        "If zero of n **independent genets** yields at least one mature seed after a specified, standardized treatment, the one-sided exact 95% upper bound for the proportion of successful genets is 1 - 0.05^(1/n). It is not a per-flower limit, and sites or cloned plants can violate independence.",
        "",
        "| Independent genets, each with comparable exposure | Observed successful genets | 95% upper bound for proportion of successful genets |",
        "|---:|---:|---:|",
    ]
    for n in (5,10,20,30,40,60):
        out.append(f"| {n} | 0 | {100*upper_limit_zero_independent_genets(n):.1f}% |")
    out += [
        "",
        "## Hypothetical chance of obtaining zero fruits at a low but nonzero fertility",
        "",
        "Assume each genet has a true fruit probability drawn from Beta(alpha,beta), with mean mu and intragenet correlation ICC. Flowers on each genet are conditionally independent given that genets's unknown propensity. Independent genets are assumed exchangeable. Zero probability = (E[(1-p_genet)^m])^n.",
        "",
        "| True mean fruit probability (hypothetical) | Genets n | Flowers per genet m | ICC | Probability of seeing zero mature fruits |",
        "|---:|---:|---:|---:|---:|",
    ]
    for mu in (.01,.03,.10):
        for genets,m in ((10,1),(20,1),(20,5),(40,5)):
            for icc in (0,.2):
                pr=zero_fruit_probability(genets,m,mu,icc)
                out.append(f"| {100*mu:.0f}% | {genets} | {m} | {icc:.1f} | {100*pr:.1f}% |")
    out += [
        "",
        "## Interpretation limits",
        "- Zero fruit can be common even if true success probability is low but nonzero, especially with few independently tagged genets. Correlated flowers are not independent biological replicates.",
        "- A positive fruit after hand-outcross proves at least some female function under those conditions, but paternity/seed identity, contamination and viability require confirmation.",
        "- An apparent zero self-fruit probability may reflect prezygotic incompatibility, late-acting incompatibility, early inbreeding depression, donor pollen loss, insufficient resource allocation, or inadequate pollination.",
        "- N for a separate pollen-tube/embryo assay cannot be counted as mature fruit trials; they are destructive sampling endpoints.",
        "- This mathematical sensitivity is not a study power analysis; it lacks real intraclass correlations, morph proportions, donor variation, populations and implementation failure.",
        "",
        "## Relevant published contrasts (not target estimates)",
        "- Shibata et al. (2018), *D. jezoensis*, many hermaphrodite genets produced no fruit even after hand outcrossing: https://doi.org/10.1007/s10265-017-0978-5",
        "- Alonso et al. (2016), *D. laureola*, flowers may accumulate pollen but show low pollen-tube conversion, emphasizing quality: https://pubmed.ncbi.nlm.nih.gov/26928007/",
        "",
    ]
    return "\n".join(out)

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--output",type=Path,default=REPORT)
    args=parser.parse_args()
    check()
    if not args.check:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(report(),encoding="utf-8")
        print("ZERO_FRUIT_SENSITIVITY_REPORT="+str(args.output),flush=True)
