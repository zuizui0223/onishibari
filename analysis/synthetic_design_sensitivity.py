#!/usr/bin/env python3
"""Illustrative clustered paired-design sensitivity, not empirical Onishibari evidence.

Within each plant, compare emasculated X+X2 versus emasculated X+S flowers:
X+X2: pollen from two unrelated donors; X+S: one unrelated donor plus self.
Both have equal total pollen loads, but X+X2 has twice the viable outcross pollen of X+S, so the simulated contrast does NOT isolate self-pollen clogging. Pilot rates from the target species are unknown.
"""
import argparse
import math
import random
from statistics import mean, stdev

def t_critical(df):
    # Cornish-Fisher approximation to the two-sided 5% Student-t quantile; df>=11.
    if df < 11:
        raise ValueError("At least 12 plants needed")
    z = 1.959963984540054
    return (z + (z**3 + z)/(4*df)
            + (5*z**5 + 16*z**3 + 3*z)/(96*df**2)
            + (3*z**7 + 19*z**5 + 17*z**3 - 15*z)/(384*df**3))

def trial(rng, plants, flowers, baseline, icc, odds_ratio):
    # Beta-distributed true plant-specific fruiting probabilities.
    concentration = 1 / icc - 1
    a, b = baseline * concentration, (1-baseline) * concentration
    differences = []
    for _ in range(plants):
        p_control = rng.betavariate(a, b)
        p_self = odds_ratio*p_control / (1-p_control+odds_ratio*p_control)
        control = sum(rng.random() < p_control for _ in range(flowers))
        treated = sum(rng.random() < p_self for _ in range(flowers))
        differences.append((treated - control) / flowers)
    effect = mean(differences)
    sd = stdev(differences)
    t = effect / (sd / math.sqrt(plants)) if sd else 0.0
    return effect, abs(t) > t_critical(plants-1)

def estimate(plants, flowers, baseline, icc, odds_ratio, reps, seed):
    rng = random.Random(seed)
    effects = []
    rejects = 0
    for _ in range(reps):
        effect, reject = trial(rng, plants, flowers, baseline, icc, odds_ratio)
        effects.append(effect)
        rejects += reject
    rate = rejects/reps
    mc_se = math.sqrt(rate*(1-rate)/reps)
    return mean(effects), rate, mc_se

def report(reps, seed):
    scenarios = [
        (0.20, 0.20, 0.50, "large_contrast"),
        (0.20, 0.20, 0.70, "moderate_contrast"),
        (0.20, 0.20, 1.00, "null"),
        (0.20, 0.40, 0.50, "large_contrast_high_heterogeneity"),
    ]
    output = [
        "# Matched-pollen, plant-paired sensitivity scenarios",
        "",
        "**Synthetic scenarios, not observations of D. pseudomezereum.**",
        "All flowers emasculated. X+X2 (two unrelated donors) and X+S (unrelated donor + self) receive equal total pollen grains. Each plant gets flowers in both groups. Different donor quality and flower-age effects must be counterbalanced in the actual experiment.",
        "",
        "Control fruit probability 20% is **hypothetical**, not measured under X+X2. For context only, hand-outcross fruit set in *D. jezoensis* was 21.1% in Shibata & Kudo (2017); pollen compositions and species are different.",
        "The model draws plant-specific baseline probabilities from a Beta distribution (ICC listed), with constant X+S-versus-X+X2 odds ratio conditional on baseline. Testing uses an approximate two-sided plant-level paired t-test (nominal alpha=0.05), not a mixed-model confirmatory analysis.",
        "",
        "| ICC | OR X+S/X+X2 | Scenario | Plants | Flowers/arm/plant | Simulated difference (pp) | Rejection % | MC SE (pp) |",
        "|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for i, (p, icc, effect_or, label) in enumerate(scenarios):
        for plants in (20, 40, 60):
            for flowers in (4, 8):
                d, power, error = estimate(plants, flowers, p, icc, effect_or, reps,
                                           seed + 1000*i + 17*plants + flowers)
                output.append(f"| {icc:.2f} | {effect_or:.2f} | {label} | {plants} | {flowers} | {100*d:+.1f} | {100*power:.1f} | {100*error:.1f} |")
    output += [
        "",
        "## Limitations and interpretation",
        "- Simulated rejection is an approximate sensitivity indicator, not observed power or a guarantee of significance.",
        "- The null scenario should reject in about 5% of simulations; Monte Carlo error is reported in the last column.",
        "- CRITICAL: X+S vs X+X2 also differs in the amount of active outcross pollen, despite equal total grains. These odds ratios are hypothetical mixture-treatment effects, NOT isolated self-pollen interference.",
        "- Even with equal total grains, X+S vs X+X2 does not separate stigma clogging from self-fertilization followed by embryo abortion. Add fluorescent pollen-tube observation and genetically assigned offspring.",
        "- Compare against X-only to measure pollen dose, and sham/intact versus emasculated to test the naturally present self-pollen pathway; latter contrast remains vulnerable to emasculation injury.",
        "- Multiple populations require explicit population-level heterogeneity; these simulations assume one exchangeable plant pool and can be optimistic.",
        "- Formal pilot-informed preregistration should use hierarchical beta-binomial/logistic models and confidence/credible intervals, plus multiple comparisons control.",
        "- Onishibari has **no data** in this file. A morphology record or close-relative fruit rate is not a measured experimental effect.",
        "",
    ]
    return "\n".join(output)

def check():
    assert 1.96 < t_critical(39) < 2.1
    first = estimate(20, 4, 0.20, 0.20, 1.0, 100, 12345)
    second = estimate(20, 4, 0.20, 0.20, 1.0, 100, 12345)
    assert first == second
    assert 0 <= first[1] <= 1
    assert estimate(20, 8, 0.20, 0.20, 0.5, 100, 67890)[0] < 0
    print("Synthetic sensitivity checks passed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    parser.add_argument("--simulations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        check()
    if args.simulations < 100:
        parser.error("--simulations must be >=100")
    result = report(args.simulations, args.seed)
    if args.output:
        from pathlib import Path
        p = Path(args.output)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(result, encoding="utf-8")
    elif not args.check:
        print(result)
