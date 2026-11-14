#!/usr/bin/env python3
"""Factorial pollination × mate-availability toy model, not target-species inference.

Repairs an earlier model that treated mate scarcity and pollinator absence as
interchangeable. One empirical INPUT ONLY: post-seed viability of a DIFFERENT
species D. jezoensis (6/19 selfed vs 31/40 outcrossed 3-year survivors).
All other parameters are hypothetical, especially autonomous selfing ability
of D. pseudomezereum which has NOT been experimentally verified.

Unit: normalized expected parental nuclear-transmission equivalents.
This is NOT a full gynodioecy population-genetic stability model.
"""
import argparse
import math
from pathlib import Path
import random
from statistics import mean
from retained_female_function_scenarios import source_viability

DEFAULT = Path(__file__).resolve().parents[1]/"docs/2026-11-14_two_channel_retention_sensitivity.md"
SEED = 20261114

def delta(p, v, s_visited, s_autonomous, f_outcross, c_male, k_fixed, r):
    """Expected retention minus complete male-function-only morph in a 2x2 table.

    p: probability compatible mate is absent.
    v: probability EFFECTIVE pollinator service is present (same in both mate
       states for toy example; visitation != service). Values assumed independent.
    When compatible mates AND effective pollination present, residual morph
       produces f_outcross extra matured outcross seeds (normalized one gene copy)
       but loses c_male outcross offspring sired compared with hypothetical male.
    When mates absent and effective visitor present, residual morph produces
       s_visited selfed mature seeds, discounted by r and doubled for 2 copies.
    Whenever effective visitors absent, produces s_autonomous selfed seeds
       with the same assumed amount regardless of whether mates were present.
    Fixed maintenance/future cost k_fixed applies in all states.
    Does not include selfing under mate-rich visitor-present conditions, pollen
       discounting, clone identity, late-age effects or genetic sex determination.
    """
    pars = (p,v,s_visited,s_autonomous,f_outcross,c_male,k_fixed,r)
    if not all(math.isfinite(x) for x in pars):
        raise ValueError("All parameters must be finite")
    if not 0<=p<=1 or not 0<=v<=1 or min(s_visited,s_autonomous,
                                            f_outcross,c_male,k_fixed)<0 or r<=0:
        raise ValueError("Parameters outside allowed range")
    return ((1-p)*v*(f_outcross-c_male)
            + 2*r*(p*v*s_visited+(1-v)*s_autonomous)
            - k_fixed)

def threshold(v, s_visited, s_autonomous, f_outcross, c_male, k_fixed,r):
    """Classify sign of delta and crossing point within p=0..1."""
    at_zero = delta(0,v,s_visited,s_autonomous,f_outcross,c_male,k_fixed,r)
    at_one = delta(1,v,s_visited,s_autonomous,f_outcross,c_male,k_fixed,r)
    slope=at_one-at_zero
    eps=1e-12
    if abs(slope)<eps:
        return ("positive_everywhere" if at_zero>eps else
                "neutral_everywhere" if abs(at_zero)<=eps else
                "negative_everywhere"),None
    root = -at_zero/slope
    if at_zero>eps and at_one>eps:
        return "positive_everywhere",None
    if at_zero<-eps and at_one<-eps:
        return "negative_everywhere",None
    if slope>0 and root>=0 and root<=1:
        return "positive_if_p_above",max(0.0,min(1.0,root))
    if slope<0 and root>=0 and root<=1:
        return "positive_if_p_below",max(0.0,min(1.0,root))
    # Handles boundary contact and floating-point rounding.
    if at_zero>=-eps and at_one>=-eps:
        return "nonnegative_everywhere",None
    if at_zero<=eps and at_one<=eps:
        return "nonpositive_everywhere",None
    raise AssertionError("Unexpected threshold geometry")

def quantile(vals,q):
    values=sorted(vals)
    at=(len(values)-1)*q
    lo=math.floor(at); hi=math.ceil(at)
    return values[lo]*(hi-at)+values[hi]*(at-lo) if hi!=lo else values[lo]

def propagate_seed_count_uncertainty(n=10000,seed=SEED,params=None):
    """Jeffreys-Beta exploratory sensitivity; treats seeds as independent.

    For sensitivity ONLY, NOT a frequentist CI/validated hierarchical posterior.
    Mothers are shared/unknown and individual genets not available.
    """
    rng=random.Random(seed)
    if params is None:
        params=dict(v=.8,s_visited=.5,s_autonomous=0.,
                    f_outcross=.01,c_male=.08,k_fixed=.01)
    ratios=[]
    breaks=[]
    for _ in range(n):
        vs=rng.betavariate(6+.5,19-6+.5)
        vo=rng.betavariate(31+.5,40-31+.5)
        r=vs/vo
        ratios.append(r)
        status,t=threshold(r=r,**params)
        if status!="positive_if_p_above":
            raise ValueError("Scenario changes threshold geometry; revise sensitivity")
        breaks.append(t)
    return (tuple(quantile(ratios,q) for q in (.05,.5,.95)),
            tuple(quantile(breaks,q) for q in (.05,.5,.95)))

def check():
    r,_,_=source_viability()
    assert abs(r-240/589)<1e-14
    assert abs(delta(p=.6,v=0,s_visited=3,s_autonomous=0,
                     f_outcross=.2,c_male=.3,k_fixed=.01,r=r)+.01)<1e-12
    assert threshold(v=0,s_visited=.5,s_autonomous=0,f_outcross=.01,
                     c_male=.08,k_fixed=.01,r=r)[0]=="negative_everywhere"
    assert threshold(v=0,s_visited=.5,s_autonomous=0,f_outcross=.01,
                     c_male=.08,k_fixed=0,r=r)[0]=="neutral_everywhere"
    t=threshold(v=.8,s_visited=.5,s_autonomous=0,f_outcross=.01,
                c_male=.08,k_fixed=.01,r=r)
    assert t[0]=="positive_if_p_above" and .16<t[1]<.19
    assert abs(delta(p=t[1],v=.8,s_visited=.5,s_autonomous=0,
                     f_outcross=.01,c_male=.08,k_fixed=.01,r=r))<1e-10
    assert threshold(v=.8,s_visited=0,s_autonomous=0,f_outcross=.01,
                     c_male=.08,k_fixed=.01,r=r)[0]=="negative_everywhere"
    assert threshold(v=.8,s_visited=.5,s_autonomous=0,f_outcross=.12,
                     c_male=.08,k_fixed=.01,r=r)[0]=="positive_everywhere"
    x=propagate_seed_count_uncertainty(500,123)
    assert x==propagate_seed_count_uncertainty(500,123)
    print(f"TWO_CHANNEL_CHECK=SUCCESS r={r:.6f} break_even_p={t[1]:.6f}")

def write_report():
    r,vo,vs=source_viability()
    scenarios=[
        ("Baseline: effective pollination 80%", dict(v=.8,s_visited=.5,s_autonomous=0,f_outcross=.01,c_male=.08,k_fixed=.01)),
        ("Low effective pollination 20%", dict(v=.2,s_visited=.5,s_autonomous=0,f_outcross=.01,c_male=.08,k_fixed=.01)),
        ("No effective pollination", dict(v=0,s_visited=.5,s_autonomous=0,f_outcross=.01,c_male=.08,k_fixed=.01)),
        ("Selfing requires service but visitors absent", dict(v=0,s_visited=1,s_autonomous=0,f_outcross=0,c_male=0,k_fixed=0)),
        ("Visitor selfing absent",dict(v=.8,s_visited=0,s_autonomous=0,f_outcross=.01,c_male=.08,k_fixed=.01)),
        ("Possible autonomous selfing (hypothesis)", dict(v=.2,s_visited=.5,s_autonomous=.2,f_outcross=.01,c_male=.08,k_fixed=.01)),
        ("More mate-rich outcross female output", dict(v=.8,s_visited=.5,s_autonomous=0,f_outcross=.10,c_male=.08,k_fixed=.01)),
        ("Higher maintenance cost",dict(v=.8,s_visited=.5,s_autonomous=0,f_outcross=.01,c_male=.08,k_fixed=.05)),
    ]
    lines=[
        "# Residual female function: separate mate shortage from pollinator-service failure",
        "",
        "**Toy model; NOT a fitted Onishibari model, evolutionary equilibrium or an estimate of mate scarcity.**",
        "",
        "### Why the original 10.9% threshold is NOT a biological prediction",
        "",
        "Shibata et al. (2018) reported that even naturally selfed offspring of the related *D. jezoensis* required pollinator service, and usual populations did NOT gain reproductive assurance. Thus it is invalid to assume mate scarcity automatically generates seed in the absence of effective visitors. The older one-state model equated these two failures. Here they are independent hypothetical axes.",
        "",
        f"One empirical viability input ONLY from a different species: Shibata (2020) *D. jezoensis*, outcrossed 31/40 three-year survivors, selfed 6/19. Relative post-seed viability **r={r:.5f}**. No one has measured the target *D. pseudomezereum* value.",
        "",
        "### Four conditional states",
        "",
        "| Compatible mates? | Effective pollinator service? | Hypothetical extra seed/father contribution of pistil retention |",
        "|---|---|---|",
        "| Available | Available | Extra outcross female output F; costs c in extra male progeny sired relative to a hypothetical fully male morph |",
        "| Unavailable | Available | Selfed mature seeds S_visit, only if actual reproductive pollen transfer occurs |",
        "| Available | Unavailable | Selfed mature seeds S_auto, **only** if autonomous selfing genuinely exists |",
        "| Unavailable | Unavailable | Selfed mature seeds S_auto, **only** if autonomous selfing genuinely exists |",
        "",
        "Let p=P(compatible mates unavailable), v=P(effective pollinator service available), taken **independent only for this toy**; F=extra matured outcross seeds in mate+visitor state; c=loss of male-sired viable-outcross progeny in mate+visitor state; k=fixed future/structural costs; S_visit & S_auto=mature selfed seed counts in the appropriate state. Normalize outcross seed offspring survival to one gene-copy equivalent. The crude two-copy selfed nuclear contribution receives a viability discount r.",
        "",
        "    ΔG = (1-p)v(F-c) + 2r[p v S_visit + (1-v) S_auto] - k",
        "",
        "The two-sex nuclear-copy convention is illustrative, not a model of allele-frequency change. It omits pollen/seed discounting, frequency dependence, maternal and paternal survival covariance, environmental correlations, and cytoplasmic male sterility.",
        "",
        "### Numerical scenarios: ALL parameters beyond r are invented",
        "",
        "| Scenario | v | S_visit | S_auto | F | c | k | Region with ΔG positive |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for name,params in scenarios:
        st,root=threshold(r=r,**params)
        label={
           "positive_if_p_above":f"p > {100*root:.1f}%",
           "positive_if_p_below":f"p < {100*root:.1f}%",
           "negative_everywhere":"never (all p)",
           "positive_everywhere":"all p",
           "neutral_everywhere":"never positive; neutral",
           "nonnegative_everywhere":"nonnegative incl. boundary",
           "nonpositive_everywhere":"never positive; boundary zero",
        }[st]
        lines.append("| "+name+" | "+" | ".join(f"{params[x]:.2f}" for x in
                     ("v","s_visited","s_autonomous","f_outcross","c_male","k_fixed"))+" | "+label+" |")
    ratio_q,threshold_q=propagate_seed_count_uncertainty()
    lines += [
        "",
        "### Uncertainty in the *related-species* seed viability input",
        "",
        f"From 6/19 selfed and 31/40 outcrossed survivors, an exploratory independent-seed Jeffreys-Beta simulation (10,000 seeded draws) gives r quantiles (5/50/95%): **{ratio_q[0]:.3f}/{ratio_q[1]:.3f}/{ratio_q[2]:.3f}**. When passed through the hypothetical baseline scenario, the 5/50/95% threshold draws are **{100*threshold_q[0]:.1f}% / {100*threshold_q[1]:.1f}% / {100*threshold_q[2]:.1f}%**.",
        "",
        "**This is not a valid cross-species posterior or confidence interval for Onishibari.** It assumes seed independence even though maternal families were pooled and ignores seed-formation and lifelong fitness. It demonstrates that treating r=0.407 as exact is misleading.",
        "",
        "### Biologically decisive counterexamples",
        "",
        "- **If autonomous selfing S_auto=0 and effective visitors v=0**, neither sexual mate absence nor the presence of a pistil creates any viable self seed. With k>0, the model favors the fully male comparator; with k=0, neutral. This is why 'pollinator scarcity' must be distinguished from 'absence of compatible outcross mates'.",
        "- **If S_visit=0 and S_auto=0**, mate scarcity cannot generate a selfing benefit. Residual pistils might still contribute outcross female success when mates and pollinators are available.",
        "- **If F>c and k is sufficiently small**, residual female function can pay off even without any mate scarcity. Thus failure to find effective selfing would not prove the pistil has no advantage.",
        "- **If neither selfing nor outcross seed benefits and k≈0**, near-neutral developmental carryover could maintain morphology without invoking a benefit, but the model cannot establish developmental constraint or selective history.",
        "",
        "### Empirical discriminating measurements in *D. pseudomezereum*",
        "",
        "1. In independently tagged morphs, cross-pollinate both forms with saturated compatible pollen, and track pollen tubes, embryos, mature seeds, germination (female outcross F).",
        "2. Bag without manipulation; bag and hand-self; bag and hand-cross; compare viable seed outcomes. Bag-no-manipulation isolates **autonomous** selfing from self-compatibility that still needs pollinator transfer.",
        "3. Record local co-flowering partner abundance and effective pollen arrival separately from fly/moth visitation; do not equate visits with pollen deposition. To estimate state probabilities, repeat across several weather and flowering windows and years.",
        "4. Distinguish pistil structure cost from cost of ripening fruit: fruit removal/supplemental resource manipulation measures one component, not the potential lifetime fitness difference of a pure male.",
        "5. Genotype progeny to separate selfed and outcrossed embryos; pollen tube growth alone does not identify pollen donor identity.",
        "",
        "### Literature",
        "- Shibata, Kameyama & Kudo (2018) *J Plant Res*. Pollinator service required despite natural selfing; no reproductive assurance: https://doi.org/10.1007/s10265-017-0978-5",
        "- Shibata (2020) PhD thesis, related-species seed viability: https://doi.org/10.14943/doctoral.k14152",
        "- *Silene nutans*: selfing can change female advantage depending on pollen limitation: https://doi.org/10.1002/ece3.1410",
        "- European chestnut: self pollen may prevent outcross fertilization via late-acting mechanisms, illustrating need to check stages: https://doi.org/10.1002/ajb2.16329",
        "",
    ]
    return "\n".join(lines)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--output",type=Path,default=DEFAULT)
    args=parser.parse_args()
    check()
    if args.check:
        return
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(write_report(),encoding="utf-8")
    print("TWO_CHANNEL_REPORT="+str(args.output),flush=True)

if __name__=="__main__":
    main()
