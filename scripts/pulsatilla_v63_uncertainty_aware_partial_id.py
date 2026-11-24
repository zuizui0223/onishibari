#!/usr/bin/env python3
"""v0.63 *Pulsatilla alpina* uncertainty-aware original cohort identification.

Source: Chen & Pannell (2024) J Ecol DOI 10.1111/1365-2745.14378;
the overlapping 2025 Current Biology endpoints and original single-flower
2022 census, exact SHA-pinned by v61 (not independent biological replicates).

For 81 original randomized bisexual single flowers (54 intact;15 half;12 all)
only 59 typed-seed flowers have EXACT original census ID and mature seed
counts (42/11/6). Another 5 have mature outcomes with untyped selfing
fractions (3/1/1); 17 no source mature outcome matches (9/3/5).

This script JOINTLY audits:
 A) deterministic v62 "sharp" bounds at delta x assumed K, conditional on
    observed selfed fraction at typed flowers (NOT a sampling CI);
 B) Beta-Jeffreys posterior seed-origin uncertainty conditional on original
    observed typed flowers (NO flowery resampling);
 C) sensitivity combining flower bootstrap + seed-origin Beta draws at
    typed flowers (not the fixed 81 enrolled cohort's design-based CI).
The untyped known mature yields remain bounded [(1-delta)M,M]; lost
outcomes each bounded [0,K], with K an explicit HYPOTHETICAL missing seed
cap, NOT a measured pistil number. Seed paternity/call uncertainty, δ
inference error, missing-not-at-random mechanism and actual recruitment
are not identified. "Posterior" below conditions on an arbitrary
Jeffreys Beta prior and exchangeable typed seed draws.

No source original individual IDs, SSR alleles, coordinates, or primary
flower-level data are printed or saved; only aggregates.
"""
from __future__ import annotations
import argparse,json,tempfile
from pathlib import Path
import numpy as np
from pulsatilla_v61_selection_tipping import source
from pulsatilla_v62_partial_identification import main as v62_main,sharp

OUT=Path("analysis_results_pulsatilla_uncertainty_partial_id_v63")
SEED=20261124
N_DRAWS=20000
DELTA=(0.,.5,.75,.9,.95,.99,1.)
CAPS=(0,5,10,20,25,50,100)
ARMS=("intact","half","all")
INIT_N={"intact":54,"half":15,"all":12}
TYPED_N={"intact":42,"half":11,"all":6}
KNOWN_UNTYPED_MATURE={"intact":126.,"half":43.,"all":113.}
UNMATCHED_N={"intact":9,"half":3,"all":5}

def get_source():
    p,t,cb,_=source(Path(tempfile.gettempdir())/"pulsatilla_v63_original_pinned")
    t=t[t.Ind_ID.isin(p.Ind_ID)].copy()
    assert len(p)==81 and len(t)==59
    res={}
    for arm in ARMS:
        g=t[t.arm==arm].sort_values("Ind_ID")
        assert len(g)==TYPED_N[arm]
        m=np.asarray(g.Mature_seed_number,float)
        typed=np.asarray(g.Total_typed_seed,int)
        selfed=np.asarray(g.Total_selfed_seed,int)
        assert len(m)>0 and np.all(m>=0) and (typed>0).all()
        assert (selfed>=0).all() and (selfed<=typed).all()
        res[arm]={"m":m,"n":typed,"k":selfed}
    assert res["all"]["k"].sum()==0 and res["all"]["n"].sum()==49
    return res

def summaries(x):
    return {
      "mean":float(np.mean(x)),
      "median":float(np.median(x)),
      "mc_percentile95":np.quantile(x,[.025,.975]).tolist(),
      "prob_positive":float(np.mean(x>0)),
      "prob_nonpositive":float(np.mean(x<=0))
    }

def pairs_by_arm(src,rng):
    out={}
    for arm in ARMS:
        z=src[arm]
        m=z["m"];n=z["n"];k=z["k"]
        a=k+.5;b=n-k+.5
        # Each typed flower is kept with its original measured seed count.
        seed=rng.beta(a[None,:],b[None,:],size=(N_DRAWS,len(m)))
        idx=rng.integers(len(m),size=(N_DRAWS,len(m)))
        # Additional flower-resampling sensitivity, independent seed origins.
        seed_b=rng.beta(a[idx],b[idx])
        out[arm]={
          "flower_mature_seed_total_observed":float(m.sum()),
          "flower_selfed_seed_posterior_rate_sum":None,
          "source_typed_genetic_seeds":int(n.sum()),
          "source_selfed_typed_seeds":int(k.sum()),
          "seed_only_mature_total":np.full(N_DRAWS,m.sum()),
          "seed_only_weighted_self_sum":(seed*m[None,:]).sum(axis=1),
          "nested_mature_total":m[idx].sum(axis=1),
          "nested_weighted_self_sum":(seed_b*m[idx]).sum(axis=1),
        }
    return out

def bound_samples(arr,delta,K,mode):
    assert mode in ("seed_only","nested")
    if mode=="seed_only":
        nameM="seed_only_mature_total";nameS="seed_only_weighted_self_sum"
    else:
        nameM="nested_mature_total";nameS="nested_weighted_self_sum"
    def typed_total(a):
        return arr[a][nameM]-delta*arr[a][nameS]
    sa=typed_total("all");si=typed_total("intact")
    la=(sa+(1-delta)*KNOWN_UNTYPED_MATURE["all"])/INIT_N["all"]
    ua=(sa+KNOWN_UNTYPED_MATURE["all"]+
         UNMATCHED_N["all"]*K)/INIT_N["all"]
    li=(si+(1-delta)*KNOWN_UNTYPED_MATURE["intact"])/INIT_N["intact"]
    ui=(si+KNOWN_UNTYPED_MATURE["intact"]+
         UNMATCHED_N["intact"]*K)/INIT_N["intact"]
    return la-ui,ua-li

def critical_delta(arr,K,mode):
    if mode=="seed_only":
        nameM="seed_only_mature_total";nameS="seed_only_weighted_self_sum"
    else:
        nameM="nested_mature_total";nameS="nested_weighted_self_sum"
    # The worst-case lower bound is linear in δ:
    # L = [(Ma+113)/12-(Mi+126+9K)/54]
    #   +δ[-(Qa+113)/12+Qi/54].
    Ma=arr["all"][nameM];Mi=arr["intact"][nameM]
    Qa=arr["all"][nameS];Qi=arr["intact"][nameS]
    intercept=(Ma+113.)/12-(Mi+126.+9*K)/54
    slope=-(Qa+113.)/12+Qi/54
    threshold=np.full(N_DRAWS,np.nan)
    valid=slope>1e-10
    threshold[valid]=-intercept[valid]/slope[valid]
    return {
       "fraction_positive_delta_slope":float(valid.mean()),
       "fraction_delta_threshold_not_in_0_1":float(np.mean(~valid|((threshold<0)|(threshold>1)))),
       "fraction_no_feasible_delta_0_1_for_guaranteed_positive":
          float(np.mean((~valid)|((threshold>1)&valid))),
       "threshold_when_positive_slope":summaries(threshold[valid]),
       "fraction_delta95_exceeds_threshold":float(np.mean(valid & (threshold<.95)))
    }

def main():
    src=get_source()
    v62=v62_main()
    rng=np.random.default_rng(SEED)
    draws=pairs_by_arm(src,rng)
    expected={a:v62["enrollment_and_outcome_link"][a] for a in ARMS}
    assert [expected[a]["genotyped_mature_outcome"] for a in ARMS]==[42,11,6]
    assert [expected[a]["mature_outcome_known_but_selfing_untyped"] for a in ARMS]==[3,1,1]
    assert [expected[a]["no_matching_mature_outcome_in_source"] for a in ARMS]==[9,3,5]
    deltas={}
    for d in DELTA:
        b={}
        for K in CAPS:
            key=str(int(K))
            # v62 only archives a finite coarse δ grid. For additional
            # δ=0.9/0.99, reconstruct the same sharp bound from source
            # typed flower data (never interpolate an observed outcome).
            sum_all=float(np.sum(src["all"]["m"]*(1-d*src["all"]["k"]/src["all"]["n"])))
            sum_intact=float(np.sum(src["intact"]["m"]*(1-d*src["intact"]["k"]/src["intact"]["n"])))
            low=(sum_all+(1-d)*113.)/12-(sum_intact+126.+9*K)/54
            high=(sum_all+113.+5*K)/12-(sum_intact+(1-d)*126.)/54
            original={"sharp_treatment_mean_difference_lower":low,
                      "sharp_treatment_mean_difference_upper":high}
            if str(float(d)) in v62["δ_sensitivity"]:
                archived=v62["δ_sensitivity"][str(float(d))]["symmetric_missing_outcome_hypothetical_cap_scenarios"]
                verified=next(x for x in archived if x["hypothetical_cap_effective_seeds_per_untracked_all_flower"]==K)
                assert abs(low-verified["sharp_treatment_mean_difference_lower"])<1e-8
                assert abs(high-verified["sharp_treatment_mean_difference_upper"])<1e-8
            modes={}
            for mode in ("seed_only","nested"):
                lo,hi=bound_samples(draws,d,K,mode)
                assert np.all(lo<=hi+1e-8)
                modes[mode]={
                  "lower_bound_samples":summaries(lo),
                  "upper_bound_samples":summaries(hi),
                  "joint_probability_bounds_strictly_positive":float((lo>0).mean()),
                  "joint_probability_bounds_strictly_negative":float((hi<0).mean()),
                  "partial_identification_disagreement_fraction":float(((lo<=0)&(hi>=0)).mean()),
                  "source_MC_draws":N_DRAWS,
                }
            b[key]={"K_unmatched_effective_seeds_max_each":int(K),
                    "v62_point_id_lower":original["sharp_treatment_mean_difference_lower"],
                    "v62_point_id_upper":original["sharp_treatment_mean_difference_upper"],
                    "uncertainty":modes,
                    "threshold_delta":{mode:critical_delta(draws,K,mode) for mode in ("seed_only","nested")} if K in (0,10,20,25) else None}
        deltas[str(float(d))]=b
    # v62 means and v61 random seed sample source audited independently.
    at=deltas["0.95"]["0"]
    assert abs(at["v62_point_id_lower"]-3.6173834653001293)<1e-9
    assert at["uncertainty"]["seed_only"]["joint_probability_bounds_strictly_positive"]<1
    assert (deltas["0.95"]["25"]["v62_point_id_lower"]<0<
            deltas["0.95"]["25"]["v62_point_id_upper"])
    return {
      "status":"v63_Pulsatilla_seed_origin_uncertainty_on_sharp_original_treatment_bounds",
      "original_study":"10.1111/1365-2745.14378",
      "overlapping_study":"10.1016/j.cub.2025.04.022",
      "reused_v61_source_sha256":{k:v62["original_sources"][k]["sha256"] for k in ("census","phenotypes","fitness")},
      "source_first_81_enrolled_bisexual":81,
      "original_treatment_n":INIT_N,
      "exact_enrollment_linked_typed_n":TYPED_N,
      "known_untyped_mature_seed_totals":KNOWN_UNTYPED_MATURE,
      "unmatched_mature_outcome_n":UNMATCHED_N,
      "typed_seed_counts_by_arm":{a:{
          "genotyped_seeds":draws[a]["source_typed_genetic_seeds"],
          "selfed_genotyped_seeds":draws[a]["source_selfed_typed_seeds"]}
          for a in ARMS},
      "rng_seed":SEED,
      "monte_carlo_draws":N_DRAWS,
      "delta_x_missing_cap_uncertainty":deltas,
      "limits":{
        "seed_origin_Jeffreys_beta_is_prior_dependent":True,
        "seed_only_has_fixed_observed_typed_flower_mature_seed_totals":True,
        "nested_bootstrap_adds_sampling_variation_as_sensitivity_not_design_based_trial_CI":True,
        "per_seed_sampling_exchangeability_assumed":True,
        "known_untyped_mature_seeds_have_unidentified_self_fraction":True,
        "unmatched_original_flowers_have_no_individual_pistil_count_to_set_physical_K":True,
        "this_is_NOT_a_posterior_for_treatment_causal_effect":True,
        "not_true_seedling_survival_or_lifetime_fitness":True,
        "inbreeding_depression_delta_is_assumed_fixed_at_each_grid_point":True,
        "male_parentage_success_and_biomass_reallocation_not_measured":True,
        "source_original_IDs_and_genotypes_not_exported":True,
        "not_Daphne_pseudomezereum_data":True}}
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=OUT/"pulsatilla_seed_origin_missing_robustness_v63.json")
    v=ap.parse_args()
    v.output.parent.mkdir(parents=True,exist_ok=True)
    results=main()
    v.output.write_text(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    subset={}
    for delta in ("0.0","0.75","0.9","0.95","0.99","1.0"):
        q=results["delta_x_missing_cap_uncertainty"][delta]
        subset[delta]={
           k:{"v62_point":[q[k]["v62_point_id_lower"],q[k]["v62_point_id_upper"]],
              "seed_only":{a:q[k]["uncertainty"]["seed_only"][a] for a in
                     ("lower_bound_samples","joint_probability_bounds_strictly_positive")},
              "nested":{a:q[k]["uncertainty"]["nested"][a] for a in
                     ("lower_bound_samples","joint_probability_bounds_strictly_positive")}}
           for k in ("0","10","20","25","50")}
        subset[delta]["threshold_K0"]=q["0"]["threshold_delta"]
    print("PULSATILLA_V63_UNCERTAINTY_SUMMARY",json.dumps(subset,ensure_ascii=False,allow_nan=False),flush=True)
