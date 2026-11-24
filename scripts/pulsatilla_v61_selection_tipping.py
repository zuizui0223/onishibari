#!/usr/bin/env python3
"""v0.61 *Pulsatilla alpina* published 2022 stamen-removal experiment.

Question: Is the conditional positive "all stamens removed -> more outcross
maternal seed equivalents" robust to typed-seed uncertainty and selection
into genotyped/recovered endpoints, when the original single-flower bisexual
manipulation cohorts were larger?

Primary originals:
Chen & Pannell (2024), Journal of Ecology 10.1111/1365-2745.14378,
Dryad 10.5061/dryad.hmgqnk9sd, Zenodo 12748590;
Current Biology (2025) 10.1016/j.cub.2025.04.022, Zenodo 15183584;
2022 single-flower plot inventory, Zenodo 10121671.
Three source data files and ZIP member SHA256 pinned by v51. Use original
case-sensitive unique Ind_ID and original treatment labels; do not mistake
the 2 outcome records missing the original 2022 census for enrolled flowers.
Do NOT back-solve offspring selfing from companion paper's fitness estimates.
Outcross mature seed equivalents are NOT observed surviving offspring.

Audit assigned->mature outcome->genotyped stages and estimate the missing
outcross contribution required to erase the conditional complete-case gain.
A reverse tipping point is a HYPOTHETICAL missing-not-at-random sensitivity,
not a claim that losses actually had the stipulated seeds.

Two-level seed-origin and flower bootstrap in the genotyped group:
Jeffreys Beta(k+0.5, n-k+0.5) within each sampled flower, then resample
whole genotyped flowers within treatment. Does not propagate SSR parentage
genotyping errors, unknown missing fertility, or δ=.95 uncertainty.
Outputs ONLY aggregates; original IDs not saved or printed.
"""
from __future__ import annotations
import argparse,io,json,tempfile
from pathlib import Path
from zipfile import ZipFile
from hashlib import sha256
import numpy as np
import pandas as pd
from pulsatilla_selfing_refit_v51 import (
    FILES, CB_MEMBER, CB_SHA256, ARMS, read_frame, download)

OUT=Path("analysis_results_pulsatilla_selection_v61")
DELTA=(0.0,0.25,0.5,0.75,0.95,1.0)
SEED=20261124
N_BOOT=16000
N_SEED_MC=16000

def source(cache):
    cache.mkdir(parents=True,exist_ok=True)
    census=read_frame(download(cache,"census"))
    typed=read_frame(download(cache,"phenotypes"))
    typed=typed[typed.Inbreeding_depression==0].copy()
    selfing=read_frame(download(cache,"selfing"))
    with ZipFile(download(cache,"fitness")) as z:
        buf=z.read(CB_MEMBER)
    assert sha256(buf).hexdigest()==CB_SHA256
    cb=read_frame(io.BytesIO(buf))
    p=census[(census.Flower_number==1)&(census.Floral_gender=="H")].copy()
    assert len(census)==175 and len(p)==81 and len(cb)==87 and len(typed)==60 and len(selfing)==60
    assert all(x.Ind_ID.is_unique for x in (p,typed,cb,selfing))
    assert dict(p.arm.value_counts())=={"intact":54,"half":15,"all":12}
    assert dict(typed.arm.value_counts())=={"intact":43,"half":11,"all":6}
    assert (typed.Total_typed_seed>0).all()
    assert np.array_equal(typed.F_RS_ID0.to_numpy(),typed.Mature_seed_number.to_numpy()) if "F_RS_ID0" in typed.columns else True
    assert np.all(typed.Total_typed_seed>=typed.Total_selfed_seed)
    # The one mismatched ID between phenotype and selfing was already
    # cross-audited by v51 using unique organ and date descriptors. Use
    # the 2024 phenotypes file for typed seed counts / mature seed count.
    # All outcome and genotype records must have treatment consistent
    # whenever they carry the same source ID.
    joined=p[["Ind_ID","arm"]].merge(cb[["Ind_ID","arm"]],on="Ind_ID",how="inner",validate="one_to_one")
    assert (joined.arm_x==joined.arm_y).all()
    genofollow=typed[["Ind_ID","arm"]].merge(cb[["Ind_ID","arm"]],on="Ind_ID",validate="one_to_one")
    assert len(genofollow)==60 and (genofollow.arm_x==genofollow.arm_y).all()
    return p,typed,cb,joined

def arm_table(p,typed,cb):
    result={}
    lost_total=0
    for arm in ARMS:
        a=p[p.arm==arm].copy()
        g=typed[typed.arm==arm].copy()
        f=cb[(cb.arm==arm)&(cb.Floral_gender=="H")].copy()
        exact=f[f.Ind_ID.isin(a.Ind_ID)].copy()
        g_in_census=g[g.Ind_ID.isin(a.Ind_ID)]
        followed_untyped=exact[~exact.Ind_ID.isin(g.Ind_ID)].copy()
        lost=a[~a.Ind_ID.isin(exact.Ind_ID)].copy()
        assert len(g_in_census)==len(g), f"typed original enrollment mismatch in arm {arm}"
        assert all(pd.to_numeric(exact.F_RS_ID0,errors="coerce").ge(0))
        assert len(g)+len(followed_untyped)+len(lost)==len(a)
        assert len(f)==(46 if arm=="intact" else 12 if arm=="half" else 8)
        maxfollow=float(pd.to_numeric(followed_untyped.F_RS_ID0).max()) if len(followed_untyped) else None
        rs={
          "original_assigned_single_flower_bisexual":len(a),
          "genotyped_mature_outcome":len(g),
          "mature_outcome_known_but_selfing_untyped":len(followed_untyped),
          "no_matching_mature_outcome_in_source":len(lost),
          "extra_outcome_not_in_original_104_single_flower_census":len(f)-len(exact),
          "genotyped_fraction_of_original":len(g)/len(a),
          "genotyped_fraction_of_matureoutcome_matched":len(g)/len(exact),
          "known_mature_outcome_untyped_harvest_seeds_sum":
              float(pd.to_numeric(followed_untyped.F_RS_ID0).sum()),
          "known_mature_outcome_untyped_harvest_seeds_mean":
              float(pd.to_numeric(followed_untyped.F_RS_ID0).mean()) if len(followed_untyped) else None,
          "known_untyped_mature_seed_minmax":
              [float(pd.to_numeric(followed_untyped.F_RS_ID0).min()),maxfollow] if len(followed_untyped) else None,
          "typed_harvested_seeds_sum":float(g.Mature_seed_number.sum()),
          "typed_selfed_genotyped_seeds_sum":int(g.Total_selfed_seed.sum()),
          "typed_genotyped_seeds_sum":int(g.Total_typed_seed.sum()),
          "not_genotyped_original_flower_total":len(a)-len(g),
          "census_and_followup_coding_verified":True,
        }
        result[arm]=rs
        lost_total+=len(lost)
    assert sum(result[a]["original_assigned_single_flower_bisexual"] for a in ARMS)==81
    assert sum(result[a]["genotyped_mature_outcome"] for a in ARMS)==60
    assert sum(result[a]["mature_outcome_known_but_selfing_untyped"] for a in ARMS)==4
    assert sum(result[a]["no_matching_mature_outcome_in_source"] for a in ARMS)==17
    assert lost_total==17
    return result

def observed_scores(typed,delta):
    ans={}
    for arm in ARMS:
        g=typed[typed.arm==arm]
        s=np.asarray(g.Total_selfed_seed,float)/np.asarray(g.Total_typed_seed,float)
        y=g.Mature_seed_number.to_numpy(float)
        ans[arm]=y*(1-delta*s)
    return ans

def tipping(original_census,typed):
    total={arm:int((original_census.arm==arm).sum()) for arm in ARMS}
    typed_n={arm:int((typed.arm==arm).sum()) for arm in ARMS}
    rs={}
    for delta in DELTA:
        scores=observed_scores(typed,delta)
        sum_all=float(scores["all"].sum())
        sum_intact=float(scores["intact"].sum())
        nm_all=total["all"]-typed_n["all"]
        nm_intact=total["intact"]-typed_n["intact"]
        assert nm_all==6 and nm_intact==11
        # If missing all-removal flowers have mean mA adjusted maternal
        # seed equivalents and missing intact have mean mI, the original
        # assigned-arm contrast is:
        # (Sall+6*mA)/12 - (Sintact+11*mI)/54.
        # Tipping mI at that mA equals 0:
        # mI = [54/12*(Sall+6*mA)-Sintact]/11.
        threshold={}
        for miss_all in [0,10,25,50,100]:
            thr=((total["intact"]/total["all"])*(sum_all+nm_all*miss_all)-sum_intact)/nm_intact
            threshold[str(miss_all)]=float(thr)
        rs[str(delta)]={
          "typed_n":typed_n,
          "genotype_complete_case_means":
              {arm:float(np.mean(scores[arm])) for arm in ARMS},
          "typed_all_minus_intact":float(scores["all"].mean()-scores["intact"].mean()),
          "original_assigned_outcross_contrast_formula":
              f"(S_all+{nm_all}*m_all)/{total['all']} -(S_intact+{nm_intact}*m_intact)/{total['intact']}",
          "typed_sum_all":sum_all,
          "typed_sum_intact":sum_intact,
          "non_genotyped_all":nm_all,
          "non_genotyped_intact":nm_intact,
          "intact_missing_mean_for_original_contrast_zero_conditional_on_all_missing_mean":threshold,
          "adversarial_all_missing_zero_intact_missing_zero_contrast":
              float(sum_all/total["all"]-sum_intact/total["intact"]),
          "not_a_recovered_original_assigned_cohort_mean":True
        }
    return rs

def partial_outcome_bounds(p,typed,cb,arm_info):
    # Known untyped mature-seed counts can rule out some completely
    # implausible imputations, but their selfed fractions are untyped.
    by={}
    for arm in ["intact","all"]:
        g=typed[typed.arm==arm]
        q=p[p.arm==arm]
        got=cb[(cb.Floral_gender=="H")&(cb.arm==arm)&cb.Ind_ID.isin(q.Ind_ID)]
        follow_untyped=got[~got.Ind_ID.isin(g.Ind_ID)]
        vals=follow_untyped.F_RS_ID0.to_numpy(float)
        assert np.all(vals>=0)
        # δ=.95, selfing fraction could be 0..1, so a measured mature seed
        # sum contributes between 5% and 100% of its recorded value.
        by[arm]={
          "known_followed_untyped_mature_harvest_sum":float(vals.sum()),
          "known_followed_untyped_delta95_lower_if_all_selfed":float(.05*vals.sum()),
          "known_followed_untyped_delta95_upper_if_all_outcrossed":float(vals.sum()),
          "unmatched_original_census_flowers":arm_info[arm]["no_matching_mature_outcome_in_source"],
          "not_assigned_zero_unmatched_outcome_by_default":True
        }
    return by

def nested_complete_case_boot(typed):
    # For each original typed flower, self-origin counts from random typed
    # seeds are Beta sampled using Jeffreys prior; whole flowers are
    # resampled independently by treatment in the same simulation.
    rng=np.random.default_rng(SEED)
    scenarios={}
    for delta in [0.0,0.5,0.95,1.0]:
        results={}
        for arm in ["intact","half","all"]:
            d=typed[typed.arm==arm]
            n=len(d)
            k=d.Total_selfed_seed.to_numpy(int)
            nn=d.Total_typed_seed.to_numpy(int)
            mat=d.Mature_seed_number.to_numpy(float)
            sel=rng.integers(0,n,size=(N_BOOT,n))
            self_p=rng.beta(k[sel]+.5,nn[sel]-k[sel]+.5)
            sc=np.mean(mat[sel]*(1-delta*self_p),axis=1)
            results[arm]={
              "n":n,
              "model_free_assigned_genotype_outcross_point":
                  float(np.mean(mat*(1-delta*k/nn))),
              "Jeffreys_beta_and_whole_flower_bootstrap_mean":
                  float(sc.mean()),
              "Jeffreys_beta_and_whole_flower_bootstrap95":
                  np.quantile(sc,[.025,.975]).tolist(),
              "draws":N_BOOT
            }
            # carry only anonymous array in memory
            results[arm]["_draws"]=sc
        contrast=results["all"]["_draws"]-results["intact"]["_draws"]
        result={
            "delta":delta,
            "all_minus_intact_genotype_complete_case_Beta_plant_bootstrap95":
                np.quantile(contrast,[.025,.975]).tolist(),
            "all_minus_intact_genotype_complete_case_Beta_plant_bootstrap_mean":
                float(contrast.mean()),
            "all_minus_intact_bootstrap_probability_positive":float(np.mean(contrast>0)),
            "source_female_mature_seeds_not_recruitment_or_seedling_survival":True,
            "genotyped_flower_selected_subset_only":True,
            "source_hermaphrodite_flower_resample_independent_not_typed_seed_independent":True,
            "arms":{arm:{key:value for key,value in x.items() if key!="_draws"} for arm,x in results.items()}
        }
        scenarios[str(delta)]=result
    return scenarios

def main():
    cache=Path(tempfile.gettempdir())/"onishibari_pulsatilla_v61_pinned"
    p,typed,cb,joined=source(cache)
    arms=arm_table(p,typed,cb)
    scenarios=tipping(p,typed)
    # tipping expects numeric info for N assigned (the underlying function
    # only needs source enrolled arm counts); keep original assigned census.
    return {
      "status":"published_Pulsatilla_stamen_removal_selection_tipping_v61",
      "source_dois":["10.1111/1365-2745.14378","10.1016/j.cub.2025.04.022"],
      "provenance":{k:{"record":x[0],"filename":x[1],"sha256":x[2]} for k,x in FILES.items()},
      "Fitness_data_RS_zip_member_sha256":CB_SHA256,
      "source_single_flower_bisexual_original_assigned":len(p),
      "source_mature_outcomes_all_genders":len(cb),
      "source_genotyped_bisexual":len(typed),
      "source_original_assigned_outcome_overlap":len(joined),
      "cohort_flow_by_arm":arms,
      "tipping_by_delta_original_assigned":scenarios,
      "followed_untyped_mature_outcomes_delta95_partial_bounds":partial_outcome_bounds(p,typed,cb,arms),
      "nested_typed_flower_and_typed_seed_bootstrap":nested_complete_case_boot(typed),
      "limitations":{
         "one_published_population_one_season_same_as_v50_v51":True,
         "all_originally_assigned_mature_progeny_not_observed":True,
         "flower_outcomes_complete_case_not_intention_to_treat":True,
         "original_2024_and_2025_scenarios_have_nonidentical_selfing_rates":True,
         "typed_seed_origin_not_true_seedling_survival":True,
         "selfed_seed_Beta_Jeffreys_uncertainty_conditions_on_assumed_independent_genotyped_seeds":True,
         "no_SSR_genotyping_call_uncertainty_or_paternal_assignments_propagated":True,
         "missing_all_and_intact_mean_values_are_hypothetical_MNAR_scenarios":True,
         "inbreeding_depression_delta_values_not_new_2022_offspring_survival_measurements":True,
         "random_removal_is_treatment_not_random_genotyping_recovery":True,
         "not_Daphne_pseudomezereum_or_Daphne_jezoensis":True,
         "individual_plant_IDs_and_raw_primary_records_never_output":True
      }}
if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--output",type=Path,default=OUT/"pulsatilla_stamen_removal_selection_tipping_v61.json")
    opt=args.parse_args()
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    result=main()
    opt.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("PULSATILLA_V61_SELECTION_TIPPING",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
