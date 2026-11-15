#!/usr/bin/env python3
"""A strictly illustrative lifetime rare-interaction encounter model.

Annual survival after reproduction s, independent Poisson effective encounter
events of mean lambda per reproductive year, with no age/stage structure.
P(>=one encounter before death) = (1-exp(-lambda))/(1-s*exp(-lambda)).
This is not seed set, realized fitness, or demographic persistence.
Compare constant lambda with a persistent two-regime population: 50% of
individuals lambda=0, 50% lambda=2*lambda, equal overall mean intensity.
A correlated bad regime depresses lifetime success even with same mean.
"""
from pathlib import Path
import csv,json,math
OUT=Path("analysis_results_theory")
OUT.mkdir(exist_ok=True)

def p_encounter(s,lam):
    if not (0 <= s <1 and lam>=0):raise ValueError("invalid model parameters")
    if lam==0:return 0.0
    q=-math.expm1(-lam)
    return q/(1-s*(1-q))

def main():
    results=[]
    for s in [0.5,0.75,0.90,0.95]:
        for lam in [0.02,0.1,0.5]:
            stable=p_encounter(s,lam)
            persistent=.5*p_encounter(s,0.0)+.5*p_encounter(s,2*lam)
            results.append(dict(survival_s=s,effective_events_mean_lambda=lam,
               expected_reproductive_seasons=1/(1-s),
               p_at_least_one_event_constant_rate=round(stable,6),
               p_at_least_one_event_persistent_mixture=round(persistent,6),
               difference=round(stable-persistent,6)))
            if persistent>stable+1e-12:raise AssertionError("unexpected mixture bound")
    with (OUT/"rare_interaction_lifetime_sensitivity.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    notes={"model":"survival-geometric annual Poisson encounters, effective successful events",
      "formula":"[1-exp(-lambda)]/[1-s*exp(-lambda)]",
      "meaning":"probability of at least one effective interaction before death, not seed production, self-incompatibility, or allele flow",
      "assumptions":["survival independent of reproduction and encounter history",
       "Poisson events within each year","reproduction possible in every year",
       "constant case independent year encounters",
       "heterogeneous case persistent high/zero individual-specific interaction regimes with same mean λ"],
      "theory_only":True,
      "not_observed_in_Daphne_or_GloPL":True,
      "primary_prediction":"The temporal organization of rare interactions may matter as much as their mean rate for whether a strategy persists",
      "next_empirical_requirement":"Individual × year visit-effective-pollination measures, vital rates, mate availability and seed recruitment; distinguish mean from serial correlation."}
    (OUT/"model_scope.json").write_text(json.dumps(notes,indent=2)+"\n")
    print(json.dumps({"results":results,"scope":notes},indent=2))
if __name__=="__main__":main()
