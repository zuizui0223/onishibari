#!/usr/bin/env python3
"""Test whether *per-flower visitation dilution* implies fitness decline.

Pure mathematical stress test, not fitted to any plant. Decompose:
  all_visits(F) = a F^beta                    (attraction)
  quality_per_visit(F) = q0 F^(-gamma)        (effective compatible transfer)
  effective_arrivals(F) = lambda F^(beta-gamma) [lambda = a*q0]
  per-flower effective arrival hazard = effective_arrivals/F
  expected flowers with >=1 effective event = F[1-exp(-lambda*F^(beta-gamma-1))]

Random-uniform independently allocated effective contacts are assumed.
Positive observation-level beta from UK PoMS must NOT be taken to mean a
within-genet causal attraction exponent. beta=.17 below is a hypothetical
scenario with approximately the PoMS value; q0, gamma are invented.

Maternal seed output optionally also obeys a finite resource budget and
an extra per-flower construction cost, following the author's existing v0.4
model. This is NOT the measured Brassica Eq1 or an independent experiment.
"""
from __future__ import annotations
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

OUT=Path("analysis_results_visitation_quality_frontier")
OUT.mkdir(parents=True,exist_ok=True)

@dataclass(frozen=True)
class Assumptions:
    beta: float = .17             # scenario, NOT within-plant field estimate
    base_visitors: float = 4.0    # visitors expected with F=1, invented
    q0: float = .5                # compatible effective events per visit F=1, invented
    budget: float = 120.0        # arbitrary resource units
    flower_cost: float = .4      # per additional flower
    fruit_cost: float = 2.0
    ovules_per_flower: float = 10.0
    fertilization_conditional: float = .6  # illustrative not measured pollen conversion
    seed_cost: float = .5

    def validate(self):
        assert self.beta >= 0 and self.base_visitors>0 and 0<self.q0<=1
        assert self.budget>0 and self.flower_cost>0
        assert self.fruit_cost>0 and self.seed_cost>=0
        assert self.ovules_per_flower>0 and 0<self.fertilization_conditional<=1

    @property
    def lam(self):
        return self.base_visitors*self.q0

    @property
    def seed_cost_per_fruit(self):
        return self.fruit_cost+self.ovules_per_flower*self.fertilization_conditional*self.seed_cost

    @property
    def max_F(self):
        return int(self.budget/self.flower_cost)


def visit_count(F:float, t:Assumptions)->float:
    if F<=0:return 0.
    return t.base_visitors*F**t.beta


def effective_count(F:float,gamma:float,t:Assumptions)->float:
    if F<=0:return 0.
    return t.lam*F**(t.beta-gamma)


def p_effective_per_flower(F:float,gamma:float,t:Assumptions)->float:
    if F<=0:return 0.
    # Mathematically exact under Poisson and independent, uniform allocation.
    return -math.expm1(-effective_count(F,gamma,t)/F)


def serviced_flower_count(F:float,gamma:float,t:Assumptions)->float:
    return F*p_effective_per_flower(F,gamma,t)


def d_serviced_dF(F:float,gamma:float,t:Assumptions)->float:
    if F<=0:raise ValueError("F must be positive for derivative")
    k=t.beta-gamma
    x=t.lam*F**(k-1.0)
    return -math.expm1(-x)+(k-1)*x*math.exp(-x)


def capacity_in_fruits(F:float,t:Assumptions)->float:
    return max(0.0,t.budget-t.flower_cost*F)/t.seed_cost_per_fruit


def maternal_seed_index(F:float,gamma:float,t:Assumptions)->float:
    """No paternity or offspring survival; units are not observed seeds."""
    return t.ovules_per_flower*t.fertilization_conditional * min(
        serviced_flower_count(F,gamma,t),capacity_in_fruits(F,t))


def opt(gamma:float,t:Assumptions):
    valid=list(range(1,t.max_F+1))
    f_poll=max(valid,key=lambda x:(serviced_flower_count(x,gamma,t),-x))
    f_seed=max(valid,key=lambda x:(maternal_seed_index(x,gamma,t),-x))
    return {
      "gamma_hypothetical":gamma,
      "effective_elasticity_k_beta_minus_gamma":round(t.beta-gamma,4),
      "F_peak_effective_serviced_flowers_in_grid":f_poll,
      "F_peak_female_seeds_with_cost_in_grid":f_seed,
      "peak_serviced_flower_count":round(serviced_flower_count(f_poll,gamma,t),6),
      "female_seeds_at_model_optimum":round(maternal_seed_index(f_seed,gamma,t),6),
      "binding_at_female_optimum":"compatible_pollen_coverage" if
          serviced_flower_count(f_seed,gamma,t)<capacity_in_fruits(f_seed,t) else "fruit_development_budget",
      "female_optimum_not_a_species_estimate":True
    }


def checks(t:Assumptions):
    t.validate()
    # For k>=0, per-flower dilution can coexist with monotone INCREASE in
    # flowers with an event. Even k=0 is nondecreasing (asymptotes at lam).
    for gamma in (0.0,t.beta):
        prev=0.
        for F in range(1,10001):
            val=serviced_flower_count(F,gamma,t)
            assert val>=prev-1e-9,(gamma,F,val,prev)
            assert d_serviced_dF(F,gamma,t)>=-1e-12
            prev=val
    # At k=0, expected serviced flowers cannot exceed # effective events.
    for F in (1,2,9,16,64,128,1000,10000):
        assert serviced_flower_count(F,t.beta,t)<t.lam+1e-9
    # Even with increasing raw visits, sufficiently falling effectiveness
    # can produce a peak and subsequent DECREASE in plant-wide coverage.
    gamma=t.beta+.10
    assert visit_count(64,t)>visit_count(9,t)
    assert serviced_flower_count(64,gamma,t)<serviced_flower_count(9,gamma,t)
    # Check analytic derivative against finite differences over regimes.
    for gamma in (0.0,t.beta,t.beta+.10,t.beta+.40):
        for F in (1.5,4.0,9.0,64.0,200.0):
            h=1e-4*F
            numeric=(serviced_flower_count(F+h,gamma,t)-
                     serviced_flower_count(F-h,gamma,t))/(2*h)
            assert math.isclose(d_serviced_dF(F,gamma,t),numeric,rel_tol=5e-7,abs_tol=1e-9)
    assert maternal_seed_index(t.max_F,.0,t)==0
    print("TESTS_OK: k>=0 monotone, k<0 decline possible, derivative verified, finite budget",flush=True)


def main():
    t=Assumptions()
    checks(t)
    gammas=(0.0,round(t.beta,3),round(t.beta+.10,3),round(t.beta+.30,3))
    cases=[opt(g,t) for g in gammas]
    curves=[]
    for gamma in gammas:
        for F in (1,2,4,9,16,32,64,128,256,300):
            curves.append({
                "F_count_hypothetical":F,
                "gamma_unmeasured_effective_quality_loss":gamma,
                "beta_hypothetical_attraction":t.beta,
                "observed_visitors_mean_model":round(visit_count(F,t),6),
                "visitors_per_flower_model":round(visit_count(F,t)/F,6),
                "effective_events_mean_model":round(effective_count(F,gamma,t),6),
                "q_compatible_per_visit_model":round(t.q0*F**(-gamma),6),
                "p_flower_with_one_or_more_compatible_events_model":round(p_effective_per_flower(F,gamma,t),6),
                "n_flowers_with_compatible_events_model":round(serviced_flower_count(F,gamma,t),6),
                "maternal_seed_model_with_budget":round(maternal_seed_index(F,gamma,t),6),
                "provenance":"SYNTHETIC_ONLY; none of these observations were made on real plants"
            })
    with (OUT/"effective_visits_flower_coverage_synthetic_curves.csv").open("w",newline="") as f:
        dw=csv.DictWriter(f,fieldnames=list(curves[0]))
        dw.writeheader();dw.writerows(curves)
    result={
      "status":"validated_mathematical_model_not_empirical",
      "question":"Can lower visits per flower alone drive a decline in a genet's total compatible pollen coverage?",
      "answer_under_assumed_poisson_allocation":"No if total effective visit intensity is nondecreasing (k=beta-gamma>=0). Falling real-visit effectiveness (gamma>beta), or a development resource cost, can create total declines.",
      "beta_positive_cannot_be_assumed_within_genet_from_PoMS":True,
      "gamma_and_lam_unmeasured":True,
      "analytical_derivative":"1-exp(-x)+(k-1)x*exp(-x), x=lambda*F^(k-1), k=beta-gamma",
      "no_causal_inference_from_public_visitation_counts":True,
      "per_flower_effective_poisson_independence_is_an_assumption":True,
      "resource_costs_not_measured":True,
      "parameters":vars(t),
      "scenarios":cases,
      "monotonicity_check_F1_to_F10000":True,
      "do_not_use_as_fitness_estimate":True
    }
    (OUT/"quality_bottleneck_theoretical_validation.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
