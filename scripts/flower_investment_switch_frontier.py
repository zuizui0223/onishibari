#!/usr/bin/env python3
"""Reproductive-investment gate model with explicit nonempirical assumptions.

Builds on Harder & Hamnett-Mount (2026), doi:10.1002/ajb2.70204, Eq. 1a/b.
Extension: deduct floral production/maintenance cost c_display * flowers from
a fixed available reproductive resource budget BEFORE fruit/seed development.
Optional separate paternal return uses Karron & Mitchell (2012) doi:
10.1093/aob/mcr193's single-study exponent 0.46; scale is arbitrarily chosen.
This is a sensitivity analysis, NOT a fitted population/trait model.

Only public paper constants are used as named *illustrative* parameters.
No simulated row is evidence for Daphne pseudomezereum, its mating system,
resource pool, pollination service, seed production, or evolutionary optimum.
"""
from __future__ import annotations
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

OUT = Path("analysis_results_allocation_gate")
OUT.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Params:
    budget: float = 120.0             # arbitrary units; no field estimate
    floral_cost: float = 0.4           # arbitrary units per flower
    fruit_cost: float = 2.0            # arbitrary units per fruit
    seed_cost: float = 0.5             # arbitrary units per developing seed
    ovules: float = 10.0              # arbitrary ovules per flower
    z: float = 0.61                   # *illustrative* 2026 Brassica fruit conditional mean
    male_exponent: float = 0.46       # published Mimulus total siring ~ F^0.46
    male_scale: float = 8.0           # arbitrary comparative fitness scale
    male_opportunity: float = 0.6     # arbitrary, distinct from female p

    def check(self):
        assert self.budget > 0 and self.floral_cost > 0
        assert self.fruit_cost > 0 and self.seed_cost >= 0 and self.ovules > 0
        assert 0 < self.z <= 1 and 0 < self.male_exponent <= 1
        assert self.male_scale > 0 and 0 <= self.male_opportunity <= 1

    @property
    def per_fruit_dev_cost(self):
        return self.fruit_cost + self.z * self.ovules * self.seed_cost

    @property
    def max_flowers(self):
        return int(math.floor(self.budget / self.floral_cost))


def female(F: float, p: float, t: Params) -> float:
    """Original paper's min(pF, resource-limited fruits), with floral cost.

    z is treated as an *effective* constant fertilization fraction for the
    schematic model. The 0.61 source datum is conditional on sampled developing
    B. rapa fruits and cannot be extrapolated to the plant or Daphne.
    """
    if not (0 <= p <= 1):
        raise ValueError("female pollination fraction must lie in [0,1]")
    if F < 0:
        raise ValueError("negative flowers")
    remaining = max(0.0, t.budget - t.floral_cost * F)
    fruits = min(p * F, remaining / t.per_fruit_dev_cost)
    return t.z * t.ovules * fruits


def paternal(F: float, t: Params) -> float:
    """Independent male opportunity; NOT measured from the same plants."""
    if F < 0:raise ValueError("negative flowers")
    return t.male_opportunity * t.male_scale * F ** t.male_exponent


def frontier(p: float, t: Params) -> float | None:
    if p <= 0:return None  # no useful female pollen and no female fitness gain
    return t.budget / (t.floral_cost + p * t.per_fruit_dev_cost)


def expected_female(F: float, mixture: tuple[tuple[float,float],...], t: Params) -> float:
    if not math.isclose(sum(w for p,w in mixture),1.0,abs_tol=1e-10):
        raise ValueError("year-state probabilities must sum to 1")
    return sum(w * female(F,p,t) for p,w in mixture)


def optimize(mix: tuple[tuple[float,float],...], t: Params, male_weight: float = 0.0) -> dict:
    if male_weight < 0:raise ValueError("negative relative paternal weight")
    vals=[(F,expected_female(F,mix,t),paternal(F,t)) for F in range(t.max_flowers+1)]
    k=max(vals,key=lambda row:(row[1]+male_weight*row[2],-row[0]))
    return {
        "integer_F":k[0],
        "female_seeds_model_units":round(k[1],5),
        "male_index_model_units":round(k[2],5),
        "combined_index_arbitrary_weighted_units":round(k[1]+male_weight*k[2],5),
        "male_weight_not_a_real_fitness_conversion":male_weight
    }


def test(t: Params):
    t.check()
    assert female(0,.4,t)==0 and female(50,0,t)==0
    assert female(t.max_flowers,1,t)==0
    for p in (.1,.25,.4,.75,1.0):
        Fstar=frontier(p,t)
        # Analytic equality at crossover
        supply=p*Fstar
        capacity=(t.budget-t.floral_cost*Fstar)/t.per_fruit_dev_cost
        assert abs(supply-capacity)<1e-8
        # Analytic breakpoint is integer optimum to rounding tolerance
        best=optimize(((p,1.0),),t)["integer_F"]
        assert abs(best-Fstar)<=1.0, (p,best,Fstar)
        # Increasing before crossing; decreasing after crossing
        assert female(Fstar*.95,p,t)<female(Fstar,p,t)
        assert female(Fstar*1.05,p,t)<female(Fstar,p,t)
    # Jensen concavity in p: mean of capped annual outputs <= capped mean
    mixed=((.05,.5),(.75,.5))
    for F in range(1,t.max_flowers,3):
        assert expected_female(F,mixed,t)<=female(F,.4,t)+1e-9
    baseline=optimize(((.4,1.0),),t)["integer_F"]
    low_high=optimize(mixed,t)["integer_F"]
    skew=optimize(((.05,.25),(.3,.5),(.95,.25)),t)["integer_F"]
    assert low_high<baseline<skew,(low_high,baseline,skew)
    higher_male=optimize(((.4,1.0),),t,3.0)["integer_F"]
    assert higher_male>baseline
    # Separate real structural display from paternity outcome.
    assert paternal(16,t)/paternal(2,t)==(16/2)**t.male_exponent
    print("SELF_TEST_OK",{
        "analytic_discrete_optimum_close":True,
        "jensen_verified_grid":True,
        "variance_optimum_direction_not_universal":[low_high,baseline,skew],
        "male_gain_can_offset_female_loss":True
    },flush=True)


def main():
    t=Params()
    test(t)
    runs=[]
    for p in (.1,.25,.4,.75,1.0):
        r=optimize(((p,1.0),),t)
        r.update({"regime":"fixed_p","p_mean":p,"p_year_states":str(p),
                   "frontier_continuous_F":round(frontier(p,t),5),
                   "female_pollen_fraction_is_not_visitor_rate":True})
        runs.append(r)
    for label,mix in [
        ("constant_mean_0.4",((.4,1.0),)),
        ("equal_0.05_0.75",((.05,.5),(.75,.5))),
        ("equal_0.0_0.8",((0.0,.5),(.8,.5))),
        ("skew_0.05_0.3_0.95",((.05,.25),(.3,.5),(.95,.25))),
    ]:
        r=optimize(mix,t)
        r.update({"regime":label,"p_mean":round(sum(p*w for p,w in mix),5),
                "p_year_states":str(mix),
                "frontier_continuous_F":None,
                "female_pollen_fraction_is_not_visitor_rate":True})
        runs.append(r)
    for w in (0.0,1.0,3.0):
        r=optimize(((.4,1.0),),t,male_weight=w)
        r.update({"regime":"p_0.4_male_weight_"+str(w),
                "p_mean":.4,"p_year_states":"0.4",
                "frontier_continuous_F":round(frontier(.4,t),5),
                "female_pollen_fraction_is_not_visitor_rate":True})
        runs.append(r)
    fields=list(runs[0].keys())
    with (OUT/"investment_switch_model_scenarios.csv").open("w",newline="") as f:
        dw=csv.DictWriter(f,fieldnames=fields)
        dw.writeheader()
        dw.writerows(runs)
    tracks=[]
    for p in (.1,.4,.75):
        for F in (0,10,20,30,40,50,60,80,100,130,160,200,300):
            if F>t.max_flowers:continue
            tracks.append({"p_illustrative":p,
                "flowers":F,
                "seeds_model_expected":round(female(F,p,t),5),
                "paternal_index_model":round(paternal(F,t),5),
                "regime":"pollination_flower_limit" if p*F <
                    max(0,t.budget-t.floral_cost*F)/t.per_fruit_dev_cost else "development_capacity_limit",
                "p_is_not_measured_onishibari":True})
    with (OUT/"investment_switch_model_curves.csv").open("w",newline="") as f:
        dw=csv.DictWriter(f,fieldnames=list(tracks[0]))
        dw.writeheader()
        dw.writerows(tracks)
    report={
        "status":"success_theoretical_only",
        "original_reference":"Harder & Hamnett-Mount 2026 doi:10.1002/ajb2.70204 Eq. 1a/b",
        "model_extension":"flower expense drawn before seed-development capacity (not calibrated)",
        "independent_male_exponent_reference":"Karron & Mitchell 2012 doi:10.1093/aob/mcr193",
        "parameters":vars(t),
        "fruit_capacity_cost_units":t.per_fruit_dev_cost,
        "threshold":"Fstar=B/(c_flower + p*(c_fruit+z*ovules*c_seed)) for constant p>0",
        "derivative_left":"dSeeds/dF=z*ovules*p (>0 if p>0)",
        "derivative_right":"dSeeds/dF=-(z*ovules*c_flower)/(c_fruit+z*ovules*c_seed) (<0 if c_flower>0)",
        "distribution_effect":"year-state variability with same mean changes optimal flower counts in either direction, depending on full distribution; cap is concave in p",
        "no_field_estimates":True,
        "not_a_new_brassica_dataset_fit":True,
        "not_a_lifetime_model":True,
        "resource_units_arbitrary":True,
        "ovule_fertilization_source_is_conditioned_on_fruits":True,
        "male_gain_exponent_does_not_generalize_from_mimulus":True,
        "scenario_results":runs
    }
    (OUT/"investment_switch_model_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    for r in runs:print("MODEL_SCENARIO",r,flush=True)


if __name__=="__main__":
    main()
