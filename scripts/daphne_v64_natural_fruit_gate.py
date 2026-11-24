#!/usr/bin/env python3
"""v0.64 — Is residual maternal reproduction in H-morph Daphne jezoensis
a diffuse low-probability per-flower process or concentrated in rare genets?

This is a secondary reanalysis of the SAME Shibata & Kudo 2023 (AJB),
Zenodo 8253921 adult genet 2018/2020 observed natural flower and fruit
counts, pinned raw ZIP+CSV SHA by v28/v31. It is NOT an independent field
study or an experiment on focal Daphne pseudomezereum.

Nulls are exact flower/fruit marginal randomizations (multivariate
hypergeometric) within sex/year, also conditioned on original observed
fruit totals in flower-number quartile strata and coarse 10/20 m spatial
blocks with 2 origin shifts. Per flower chance is exchangeable WITHIN
stratum, NOT evidence about individual flower-level independence.

Outcome: per-genet observed 0-vs-positive mature fruits; separate from
fertilization, viable seed, recruited offspring, pollen donor fitness.
Original individual IDs, genotype, pairwise fruit edges, source coordinates
NEVER written to JSON or logs. Quantile/aggregate counts only.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,ZIP_SHA,stringify

OUT=Path("analysis_results_daphne_fruit_activation_v64")
SEED=20261124
REPS=6000
YEARS=(2018,2020)
SEXES=("H","F")
EXPECTED_H={2018:{"genets":105,"fruits":53,"flowers":4765},
            2020:{"genets":71,"fruits":28,"flowers":2573}}
EXPECTED_F={2018:{"fruits":255,"flowers":619},2020:{"fruits":205,"flowers":719}}

def make(df,year,sex):
    sx=f"Sex_{year}.g"
    cf=f"NaturalFlowerNo_{year}.g"
    cr=f"NaturalFruitNo_{year}.g"
    q=df.loc[df[sx].astype(str).eq(sex)].copy()
    q["id"]=stringify(q.GenetID2)
    q["n"]=pd.to_numeric(q[cf],errors="coerce")
    q["y"]=pd.to_numeric(q[cr],errors="coerce")
    valid=q.n.gt(0)&q.y.ge(0)&q.y.le(q.n)&q.y.notna()
    q=q.loc[valid].copy().sort_values("id")
    assert len(q)>=15 and q.id.is_unique
    n=q.n.to_numpy(np.int64)
    y=q.y.to_numpy(np.int64)
    assert np.all(n>=y) and n.sum()>y.sum()>0
    xy=q[["East.g","North.g"]].apply(pd.to_numeric,errors="coerce").to_numpy()
    assert len(xy)==len(n)
    if sex=="H":
        exp=EXPECTED_H[year]
        assert (len(q),int(n.sum()),int(y.sum()))==(exp["genets"],exp["flowers"],exp["fruits"]),("H original source changed",year,len(q),int(n.sum()),int(y.sum()))
    if sex=="F":
        exp=EXPECTED_F[year]
        assert (int(n.sum()),int(y.sum()))==(exp["flowers"],exp["fruits"])
    return {"n":n,"y":y,"xy":xy,"ids":q.id.tolist()}

def obs_metrics(y):
    y=np.asarray(y,dtype=np.int64)
    assert y.ndim==1 and y.sum()>0
    yp=y[y>0]
    return {
     "fruit_positive_genets":int((y>0).sum()),
     "zero_fruit_genets":int((y==0).sum()),
     "fruit_count_per_fruiting_genet":float(yp.mean()),
     "effective_number_of_fruit_mothers_herfindahl":float(y.sum()**2 / (y@y)),
     "share_all_fruits_top_tenth_genets":float(np.sort(y)[-max(1,int(np.ceil(len(y)*.1))):].sum()/y.sum()),
    }

def sample_within(n,y,labels,rng,draws=REPS):
    """Exact within-stratum fruit quotas conditional on observed source n,y.
    Simulated fruit allocations never exceed available per-genet flower counts.
    """
    lab=np.asarray(labels,dtype=np.int64)
    uniq=np.unique(lab)
    assert len(n)==len(y)==len(lab)
    sims=np.zeros((draws,len(n)),dtype=np.int32)
    for group in uniq:
        ix=np.flatnonzero(lab==group)
        count=int(y[ix].sum())
        if count==0:continue
        if count==int(n[ix].sum()):
            sims[:,ix]=np.tile(n[ix],(draws,1))
        elif len(ix)==1:
            sims[:,ix[0]]=count
        else:
            sims[:,ix]=rng.multivariate_hypergeometric(n[ix],count,size=draws,method="marginals")
    assert np.all(sims.sum(axis=1)==y.sum())
    assert np.all(sims<=n[None,:])
    return sims

def simulation_summary(obs,sim):
    active=(sim>0).sum(axis=1).astype(int)
    sumy=sim.sum(axis=1)
    hhi=np.square(sumy[:,None]).ravel()/(sim.astype(float)**2).sum(axis=1)
    topn=max(1,int(np.ceil(sim.shape[1]*.1)))
    share=np.sort(sim,axis=1)[:,-topn:].sum(axis=1)/sumy
    a0=int(obs["fruit_positive_genets"])
    e0=float(obs["effective_number_of_fruit_mothers_herfindahl"])
    s0=float(obs["share_all_fruits_top_tenth_genets"])
    def d(z,observed,tail):
        return {"null_mean":float(z.mean()),
                "null_95percentile":np.quantile(z,[.025,.975]).tolist(),
                "observed":float(observed),
                "null_tail_p_one_sided":float((1+np.count_nonzero(z<=observed if tail=="low" else z>=observed))/(len(z)+1)),
                "tail":tail}
    return {
      "fruit_positive":d(active,a0,"low"),
      "effective_fruit_mothers":d(hhi,e0,"low"),
      "top_10pct_fruit_share":d(share,s0,"high"),
      "draws":len(sim)
    }

def spatial(q,side,offset):
    xy=q["xy"]
    valid=np.isfinite(xy).all(axis=1)
    if not valid.all():
        return {"available":False,"fraction_coord_valid":float(valid.mean())}
    # Origin anchored within the same source sample; shifted grid tests
    # coarse space conditioned fruit allocation. Absolute coords never output.
    origin=xy.min(axis=0)
    cells=np.floor((xy-origin+offset)/side).astype(np.int64)
    _,strata=np.unique(cells,axis=0,return_inverse=True)
    return {"available":True,"strata":strata,"blocks":int(strata.max()+1),
            "fruit_blocks":int(np.unique(strata[np.asarray(q["y"])>0]).size),
            "single_genet_blocks":int((np.bincount(strata)==1).sum())}

def flower_quartiles(n):
    # rank gives equal-number count partitions despite nonunique n and
    # explicit tie-break per genet. Do NOT let outcomes decide stratum.
    rank=np.argsort(np.argsort(n,kind="stable"),kind="stable")
    return np.minimum(3,(rank*4)//len(n))

def analyze(q,year,sex,rng):
    n,y=q["n"],q["y"]
    observed=obs_metrics(y)
    n_active=observed["fruit_positive_genets"]
    naive=float(y.sum()/n.sum())
    by_sex={
     "year":year,"morph_exact_label":sex,
     "n_original_genets_valid":len(n),
     "flowers_total":int(n.sum()),"fruits_total":int(y.sum()),
     "natural_fruit_per_flower":naive,
     "fruit_positive_genets":n_active,
     "zero_fruit_genets":int((y==0).sum()),
     "mean_genet_fruit_rate":float(np.mean(y/n)),
     "observed":observed,
     "null_tests":{}
    }
    designs={"all_flowers_within_morph_year":np.zeros(len(n),dtype=int),
             "flower_count_quartile_conditioned":flower_quartiles(n)}
    for name,lab in designs.items():
        sims=sample_within(n,y,lab,rng)
        by_sex["null_tests"][name]={"strata":int(len(np.unique(lab))),
                    "summary":simulation_summary(observed,sims)}
    # Source-independent outlier guard: re-run a flower-budget conditional
    # null after withholding the largest observed maternal fruit producer.
    # This changes the population being described, and is NOT a genotype
    # dropout/assignment correction.
    largest=int(np.argmax(y))
    keep=np.arange(len(y))!=largest
    ry=y[keep];rn=n[keep]
    if ry.sum()>0:
        rob=obs_metrics(ry)
        by_sex["single_largest_fruit_genet_excluded"]={
          "n_retained_genets":len(rn),
          "n_retained_fruits":int(ry.sum()),
          "n_original_fruits_removed_by_one_genet":int(y[largest]),
          "observed":rob,
          "flower_budget_global_null":simulation_summary(
               rob,sample_within(rn,ry,np.zeros(len(rn),dtype=int),rng)),
          "flower_quartile_conditional_null":simulation_summary(
               rob,sample_within(rn,ry,flower_quartiles(rn),rng)),
          "not_outlier_independent_population_replication":True
        }
    for side in (10,20):
        for origin_shift in (0,side/2):
            s=spatial(q,side,origin_shift)
            name=f"spatial_{side}m_shift_{int(origin_shift)}m"
            if not s["available"]:
                by_sex["null_tests"][name]={"available":False,"valid_coord_fraction":s["fraction_coord_valid"]}
                continue
            sim=sample_within(n,y,s["strata"],rng)
            by_sex["null_tests"][name]={
              "available":True,"spatial_block_size_m":side,
              "shift_m":origin_shift,
              "blocks":s["blocks"],"fruit_positive_blocks":s["fruit_blocks"],
              "single_genet_blocks":s["single_genet_blocks"],
              "summary":simulation_summary(observed,sim)}
            # Simultaneous control for coarse location AND within-year
            # flower crop rank, without assigning arbitrary biological
            # neighborhoods to the raw coordinates. Many single-genet
            # strata reduce testing opportunity; report loss of support.
            halves=np.minimum(1,np.argsort(np.argsort(n,kind="stable"),kind="stable")*2//len(n))
            joint=s["strata"]*2+halves
            _,joint=np.unique(joint,return_inverse=True)
            group_sizes=np.bincount(joint)
            sampleable=np.isin(joint,np.flatnonzero(group_sizes>1))
            # Fruits fixed in singleton groups contribute no randomness.
            fruit_in_sampleable=int(y[sampleable].sum())
            simjoint=sample_within(n,y,joint,rng)
            by_sex["null_tests"][name+"_plus_flower_half"]={
              "available":True,"spatial_block_size_m":side,
              "shift_m":origin_shift,
              "blocks_after_joint_flower_size":len(group_sizes),
              "singleton_genet_strata":int((group_sizes==1).sum()),
              "n_genets_in_multigenet_strata":int(sampleable.sum()),
              "n_fruits_in_multigenet_strata":fruit_in_sampleable,
              "null_test_reduced_randomization_due_to_fine_strata":True,
              "summary":simulation_summary(observed,simjoint)}
    return by_sex

def temporal(d):
    a,b=d[2018]["H"],d[2020]["H"]
    a_ids={z:int(f>0) for z,f in zip(a["ids"],a["y"])}
    b_ids={z:int(f>0) for z,f in zip(b["ids"],b["y"])}
    # Main stable H-to-H census 45 was previously source-confirmed,
    # restrict both years to same morphological sex H, positive flowers,
    # valid fruit totals in both years.
    keep=sorted(set(a_ids)&set(b_ids))
    n=len(keep)
    x=np.array([a_ids[k] for k in keep],int)
    y=np.array([b_ids[k] for k in keep],int)
    both=int((x*y).sum())
    if n==0:return {"both_years_valid_H":0}
    expected=float(x.sum()*y.sum()/n)
    return {
      "same_H_genets_with_valid_natural_fruit_both_years":n,
      "fruit_positive_H_2018":int(x.sum()),
      "fruit_positive_H_2020":int(y.sum()),
      "positive_both_years":both,
      "null_expected_repeat_if_years_exchangeable":expected,
      "hypergeometric_one_sided_repeat_excess_p":float(hypergeom.sf(both-1,n,int(x.sum()),int(y.sum()))),
      "null_does_not_preserve_year_specific_flower_counts_or_spatial_network":True,
      "no_causal_temporal_residual_pistil_evolution":True
    }

def main():
    original=data()
    adult=original["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(adult)==492
    rng=np.random.default_rng(SEED)
    private={year:{sex:make(adult,year,sex) for sex in SEXES} for year in YEARS}
    by=[analyze(private[year][sex],year,sex,rng) for year in YEARS for sex in SEXES]
    return {"status":"original_Daphne_jezoensis_residual_female_activation_v64",
            "source_article":"10.1002/ajb2.16239",
            "relevant_separate_hand_outcross_study":"10.1007/s10265-017-0978-5",
            "original_zenodo":"10.5281/zenodo.8253921",
            "original_zip_sha256":ZIP_SHA,"rng_seed":SEED,
            "conditioned_draws_per_null":REPS,
            "sex_year":by,"strict_H_temporal_recurrence":temporal(private),
            "limits":{
              "fruit_count_mature_not_viable_seed_or_genetic_offspring":True,
              "rare_fruits_do_not_prove_pistil_intrinsic_viability":True,
              "null_assumes_exchangeable_flower_fate_within_given_strata":True,
              "source_flower_counts_not_ovule_or_pollination_day_match":True,
              "spatial_grid_conditions_coarse_coordinates_not_complete_resource_and_pollen_exposure":True,
              "source_prior_related_parentage_not_used_in_natural_fruit_counts":True,
              "observed_fruit_clustering_not_equivalent_to_mechanistic_two_stage_reproductive_failure":True,
              "Shibata2018_hand_outcross_nonfruiting_is_previously_published":True,
              "no_focal_D_pseudomezereum_flower_or_pollen_tube_data":True,
              "one_source_population_repeated_years_not_two_independent_replicates":True,
              "no_original_genetIDs_alleles_coords_or_parentage_rows_exported":True}}
if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--output",type=Path,default=OUT/"daphne_natural_fruit_activation_conditioned_v64.json")
    x=arg.parse_args()
    x.output.parent.mkdir(parents=True,exist_ok=True)
    result=main()
    x.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("DAPHNE_V64_FRUIT_GATE_AGGREGATES",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
