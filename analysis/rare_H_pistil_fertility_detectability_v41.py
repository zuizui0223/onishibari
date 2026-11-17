#!/usr/bin/env python3
"""v0.41 | Prospective rare functional pistil detectability for Onishibari, NOT data.

For D. pseudomezereum H-morph pistils, a successful trial is defined *strictly*
as at least one independently verified viable, genetically assigned maternal
seed from a bud under bagged hand-outcross pollination. A mature fruit, pollen
tube, unfathered seed or H-morph label alone is NOT a confirmed success.

No source studies supply Onishibari p, so p=.01/.03/.05/.1 and within-genet
intra-class correlations 0/.3/.6 are HYPOTHETICAL design sensitivities.
The independent unit is a genet, not an individual flower; beta-binomial
mixture approximates exchangeable within-genet outcomes.

Under the existing prospective 2027 five-arm conservation-safe pilot, first
hand-cross attempt uses one of five manipulated buds. n outcross flowers/genet
needs 4+n manipulated buds, <=10 and <=50% verified buds, requiring at least
2*(4+n) original verified buds. Not a simulation of actual source observations,
not power for overall sexual-system evolution or multi-year viability.
"""
from __future__ import annotations
import csv,json,math
from pathlib import Path
from itertools import product

OUT=Path("analysis_results_onishibari_rare_fertility_v41")
PS=(.01,.03,.05,.10)
ICS=(0.0,.3,.6)
GS=(10,20,40,60,100,160)
NS=(1,2,3)
ALPHA=.05

def no_success_one_genet(p:float,n:int,rho:float)->float:
    assert 0<=p<=1 and n>=1 and 0<=rho<1
    if p==0:return 1.0
    if p==1:return 0.0
    if rho==0:return (1-p)**n
    a=p*(1-rho)/rho
    b=(1-p)*(1-rho)/rho
    result=1.0
    for k in range(n):
        result*=(b+k)/(a+b+k)
    return result

def probabilities(p,n,rho,G):
    no=no_success_one_genet(p,n,rho)
    yes=1-no
    p0=no**G
    p1=G*yes*(no**(G-1))
    return {"p_zero_positive_genets":p0,
      "p_at_least_one_positive_genet":1-p0,
      "p_at_least_two_positive_genets":max(0.,1-p0-p1),
      "expected_positive_genets":G*yes,
      "per_genet_positive_probability":yes}

def required_genets(p,n,rho,n_successes,confidence=.95,max_G=20000):
    for G in range(1,max_G+1):
        m=probabilities(p,n,rho,G)
        if m["p_at_least_one_positive_genet" if n_successes==1 else
             "p_at_least_two_positive_genets"]>=confidence-1e-12:return G
    return None

def upper_p_if_zero(G,n,rho):
    # One-sided alpha=0.05 bound solving Pr(zero positive genets | p)=.05.
    low,high=0.,1.
    for _ in range(90):
        mid=(low+high)/2
        if no_success_one_genet(mid,n,rho)**G>ALPHA:
            low=mid
        else:high=mid
    return (low+high)/2

def record(p,n,rho,G):
    z=probabilities(p,n,rho,G)
    return {"probability_success_per_flower_hypothetical":p,
       "genets":G,"hand_cross_flowers_per_genet":n,
       "within_genet_ICC_hypothetical":rho,
       "n_total_hand_cross_flowers":n*G,
       "assigned_arms_per_genet_5arm_plus_cross":4+n,
       "minimum_verified_buds_per_genet_for_half_cap":2*(4+n),
       "cumulative_manipulated_buds_per_genet_le10":(4+n)<=10,
       **{k:round(v,8) for k,v in z.items()}}

def main():
    OUT.mkdir(exist_ok=True)
    assert no_success_one_genet(.03,1,0.0)==.97
    assert abs(no_success_one_genet(.03,2,.3)-((.97**2)+.3*.03*.97))<1e-12
    assert abs(probabilities(.03,1,0,20)["p_at_least_one_positive_genet"]-(1-.97**20))<1e-12
    assert required_genets(.03,1,0,1)==99
    assert required_genets(.03,1,0,2)==157
    assert required_genets(.05,1,0,1)==59
    assert required_genets(.05,1,0,2)==93
    assert required_genets(.03,2,0,1)==50
    assert required_genets(.03,2,0,2)==79
    assert required_genets(.03,2,.3,1)>=50
    assert required_genets(.03,2,.6,1)>=required_genets(.03,2,.3,1)
    assert abs(upper_p_if_zero(20,1,0)-(1-ALPHA**(1/20)))<1e-10
    rows=[record(p,n,rho,G) for p,n,rho,G in product(PS,NS,ICS,GS)]
    path=OUT/"onishibari_2027_rare_H_outcross_viable_seed_detectability_v41.csv"
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    req=[]
    for p,n,rho in product(PS,NS,ICS):
        req.append({"hypothetical_p_per_flower":p,
          "hand_cross_flowers_per_independent_genet":n,
          "genet_ICC":rho,
          "min_buds_required_for_five_arm_safe_pilot":2*(4+n),
          "genets_for_95pct_detect_at_least_one":required_genets(p,n,rho,1),
          "genets_for_95pct_detect_at_least_two_independent_genets":required_genets(p,n,rho,2)})
    upper=[]
    for G,n,rho in product((10,20,40,60),(1,2),ICS):
        upper.append({"n_genets_with_zero_successes":G,
          "hand_cross_flowers_per_genet":n,
          "assumed_genet_ICC":rho,
          "conditional_one_sided_95pct_upper_bound_p_per_flower":round(upper_p_if_zero(G,n,rho),7)})
    fpc={p:probabilities(p,1,0,20) for p in (.05,.20,.50)}
    output={
      "status":"PASS_ONISHIBARI_PROSPECTIVE_RARE_H_FERTILITY_DETECTABILITY_V41",
      "not_a_measured_onishibari_sex_fertility_probability":True,
      "endpoint":"genet_with_at_least_one_verified_viable_maternal_seed_after_bagged_hand_outcross",
      "success_requires_genetic_maternity_viability_fruit_alone_insufficient":True,
      "independent_units_are_genets":True,
      "hypothetical_flower_probabilities":list(PS),
      "hypothetical_within_genet_ICC":list(ICS),
      "sample_size_grid":list(GS),
      "flowers_per_genet_grid":list(NS),
      "current_five_arm_budget":{
         "first_four_non_cross_arms":"open,bag_no_hand,bag_sham,bag_hand_self",
         "one_extra_or_more_bag_hand_cross":"n=1,2,3",
         "max_cumulative_manipulations_per_genet":10,
         "max_manipulated_fraction_of_verified_buds":.5,
         "bud_minimum_per_genet_for_n1":10,
         "bud_minimum_per_genet_for_n2":12,
         "bud_minimum_per_genet_for_n3":14,
         "independence_of_genet_may_require_marker_validation":True},
      "required_independent_genets_for_detection":req,
      "zero_success_upper_bound_scenarios":upper,
      "not_observed_female_F_positive_control_sensitivity":{
        str(k):{n:round(v,8) for n,v in z.items()}
        for k,z in fpc.items()},
      "highlight_example_p03_G20_n1_rho0":record(.03,1,0,20),
      "highlight_example_p03_G20_n2_rho03":record(.03,2,.3,20),
      "limitations":{
        "hand_cross_success_probability_unknown_in_focal":True,
        "no_translated_jezoensis_p_or_cross_species_trait_substitution":True,
        "mother_and_donor_compatibility_and_pollen_viability_are_essential_controls":True,
        "no_automatic_selfing_or_retained_pistil_adaptive_fitness_inferred":True,
        "small_number_of_genets_with_zero_success_does_not_prove_sterility":True,
        "even_verified_rare_maternal_success_does_not_prove_selection_for_retention":True,
        "genet_level_responses_may_share_population_and_year_environment":True,
        "same_young_seed_not_final_established_adult_fitness":True,
        "budget_examples_presume_sufficient_verified_predehiscence_buds":True,
        "this_is_an_illustrative_power_and_detection_model_not_observed_results":True}}
    (OUT/"onishibari_2027_focal_H_rare_maternal_detectability_v41.json").write_text(
       json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    compact={"status":output["status"],
       "p03_G20_n1":output["highlight_example_p03_G20_n1_rho0"],
       "p03_G20_n2_rho03":output["highlight_example_p03_G20_n2_rho03"],
       "minimum_required_for_p03":
         [v for v in req if v["hypothetical_p_per_flower"]==.03],
       "zero_success_p_upper_G20n1":next(v for v in upper if v["n_genets_with_zero_successes"]==20
           and v["hand_cross_flowers_per_genet"]==1 and v["assumed_genet_ICC"]==0),
       "n_sensitivity_grid_rows":len(rows),
       "not_observed_Onishibari_outcomes":True}
    print("ONISHIBARI_RARE_FERTILITY_DESIGN_V41",json.dumps(compact,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
