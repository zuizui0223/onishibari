#!/usr/bin/env python3
"""v0.65: Do zero-fruiting H genets require a discrete sterile subpopulation?

Original Shibata & Kudo (2023) AJB, DOI 10.1002/ajb2.16239,
author Zenodo 10.5281/zenodo.8253921, source adult data raw ZIP SHA pinned
by v28. Same genets as v64; NOT another source/field experiment and NOT
Daphne pseudomezereum.

Independent primary unit = adult genet in each morph x year (H/F, 2018/2020).
Observed (natural fruits y_i) <= (natural flowers n_i). Compare:
  Binomial: one common per-flower fruit probability p;
  Beta-binomial: genotype-specific latent fruit chance varying continuously;
  Zero-inflated binomial: a point-mass of structurally zero genets pi;
  Zero-inflated beta-binomial: mixture of both (full-fit only diagnostic).
The 2-parameter BB and 2-parameter ZIB can both produce many zeros.
10-fold repeated 5-fold heldout-GENET CV distinguishes predictive fit from
in-sample zero description. Report same-source bootstrap uncertainty of
paired cross-validated per-genet predictive log score differences.
A positive ZIB pi never proves permanent anatomical sterility.

Models do not include donor pollen compatibility, phenology, ovule counts,
plant resource or precise position: v64 conditioned hypergeometric nulls
separately test crop and 10/20m spatial cells. Repeated years share some
genets and must not be counted as independent populations.
Do not export original IDs, coordinates, genotype/individual source rows.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit,logit,betaln,gammaln,logsumexp
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data,ZIP_SHA
from daphne_v64_natural_fruit_gate import make,YEARS,SEXES

OUT=Path("analysis_results_daphne_zero_activation_v65")
MODELS=("binomial","beta_binomial","zero_inflated_binomial")
FULL_DIAGNOSTIC="zero_inflated_beta_binomial"
SEED=20261124
CV_REPEATS=2
CV_FOLDS=5
BOOT=8000
EXPECTED={(2018,"H"):(105,4765,53,14),
          (2018,"F"):(63,619,255,37),
          (2020,"H"):(71,2573,28,9),
          (2020,"F"):(57,719,205,32)}

def logpmf(n,y,model,params):
    n=np.asarray(n,float)
    y=np.asarray(y,float)
    p=float(expit(params[0]))
    base=(gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1))
    if model in ("binomial","zero_inflated_binomial"):
        base+=y*np.log(p)+(n-y)*np.log1p(-p)
    else:
        conc=float(np.exp(params[1]))
        alpha=p*conc;beta=(1-p)*conc
        base+=betaln(y+alpha,n-y+beta)-betaln(alpha,beta)
    if model in ("zero_inflated_binomial","zero_inflated_beta_binomial"):
        pi=float(expit(params[-1]))
        base=np.where(y==0,
             np.logaddexp(np.log(pi),np.log1p(-pi)+base),
             np.log1p(-pi)+base)
    return base

def summarize_params(model,param):
    r={"conditional_nonstructural_fruit_probability":float(expit(param[0]))}
    if "beta_binomial" in model:
        k=float(np.exp(param[1]))
        r["latent_beta_concentration"]=k
        r["latent_beta_probability_variance"]=float(expit(param[0])*(1-expit(param[0]))/(1+k))
        r["latent_beta_intragenet_correlation_rho"]=float(1/(1+k))
    if model.startswith("zero_inflated"):
        r["hypothetical_point_zero_genet_fraction"]=float(expit(param[-1]))
    return r

def fit(n,y,model):
    n=np.asarray(n,float);y=np.asarray(y,float)
    assert np.all(n>=y) and n.sum()>0
    total_p=float(np.clip(y.sum()/n.sum(),1e-7,.99))
    p0=float(logit(total_p))
    if model=="binomial":
        return np.array([p0]),float(logpmf(n,y,model,[p0]).sum()),True
    # Limits prevent numerical overflow but deliberately allow BB dispersion
    # / zero fraction near boundaries to reveal nonidentifiability.
    starts=[]
    if model=="beta_binomial":
        bounds=[(-16,8),(-9,13)]
        starts=[[p0,np.log(k)] for k in (.025,.25,2,25,250)]
    elif model=="zero_inflated_binomial":
        bounds=[(-16,8),(-12,7)]
        starts=[[logit(min(.9,total_p/(1-pi))),logit(pi)] for pi in (.02,.3,.6,.85)]
    elif model=="zero_inflated_beta_binomial":
        bounds=[(-16,8),(-9,13),(-12,7)]
        starts=[[logit(min(.9,total_p/(1-pi))),np.log(k),logit(pi)]
                for k in (.1,5,50) for pi in (.02,.4,.8)]
    else:raise ValueError(model)
    best=None
    for initial in starts:
        res=minimize(lambda z:-float(logpmf(n,y,model,z).sum()),
                     x0=np.array(initial,float),method="L-BFGS-B",
                     bounds=bounds,options={"maxiter":350,"ftol":1e-10})
        if np.isfinite(res.fun) and (best is None or res.fun<best.fun):
            best=res
    assert best is not None and np.isfinite(best.fun)
    return np.asarray(best.x),-float(best.fun),bool(best.success)

def crossvalidate(n,y,rng):
    N=len(n)
    indices=np.arange(N)
    # Each original adult is held out in exactly one fold per repeat.
    cvlog={name:np.zeros((CV_REPEATS,N),float) for name in MODELS}
    cvzero={name:np.zeros((CV_REPEATS,N),float) for name in MODELS}
    cvfit={name:0 for name in MODELS}
    for rep in range(CV_REPEATS):
        shuffled=rng.permutation(indices)
        folds=np.array_split(shuffled,CV_FOLDS)
        assert sum(len(x) for x in folds)==N
        for held in folds:
            train=np.ones(N,bool);train[held]=False
            for model in MODELS:
                param,_,success=fit(n[train],y[train],model)
                cvfit[model]+=int(success)
                pred=logpmf(n[held],y[held],model,param)
                zero=logpmf(n[held],np.zeros(len(held),int),model,param)
                cvlog[model][rep,held]=pred
                cvzero[model][rep,held]=np.exp(zero)
    avg={z:x.mean(axis=0) for z,x in cvlog.items()}
    obszer=(y==0).astype(float)
    result={}
    for model in MODELS:
        l=avg[model]
        predzero=cvzero[model].mean(axis=0)
        result[model]={
         "heldout_genets":N,
         "avg_log_predictive_score_per_genet":float(l.mean()),
         "sum_log_predictive_score_per_genet":float(l.sum()),
         "mean_predicted_zero_genet_fraction":float(predzero.mean()),
         "actual_zero_genet_fraction":float(obszer.mean()),
         "predictive_zero_probability_brier":float(np.mean((predzero-obszer)**2)),
         "optimization_successful_fits":cvfit[model],
         "total_fits":CV_REPEATS*CV_FOLDS,
        }
    # Difference scores are paired at adult genet level. Bootstrap reflects
    # within-cohort finite genet variation but not new-year/site generality.
    def diff(a,b):
        dz=avg[a]-avg[b]
        ix=rng.integers(N,size=(BOOT,N))
        mc=dz[ix].mean(axis=1)
        return {
         "comparison":a+"_minus_"+b,
         "mean_additional_heldout_log_score_per_genet":float(dz.mean()),
         "adult_genet_bootstrap95":np.quantile(mc,[.025,.975]).tolist(),
         "positive_genets":int(np.count_nonzero(dz>0)),
         "negative_genets":int(np.count_nonzero(dz<0)),
         "bootstrap_draws":BOOT,
         "one_original_population_no_external_validation":True
        }
    result["comparisons"]={
      "beta_binomial_minus_binomial":diff("beta_binomial","binomial"),
      "zero_inflated_binomial_minus_binomial":diff("zero_inflated_binomial","binomial"),
      "zero_inflated_binomial_minus_beta_binomial":diff("zero_inflated_binomial","beta_binomial")
    }
    result["folds"]=CV_FOLDS
    result["repeats"]=CV_REPEATS
    return result

def analyze(record,year,sex,rng):
    n=record["n"];y=record["y"]
    exp=EXPECTED[(year,sex)]
    assert (len(n),int(n.sum()),int(y.sum()),int((y>0).sum()))==exp
    fits={}
    for model in list(MODELS)+[FULL_DIAGNOSTIC]:
        param,ll,ok=fit(n,y,model)
        k=len(param)
        fits[model]={
           "log_likelihood":ll,
           "AIC":float(2*k-2*ll),
           "BIC_genets_unit":float(np.log(len(n))*k-2*ll),
           "number_parameters":k,
           "optimizer_success":ok,
           "parameters":summarize_params(model,param),
           "probability_zero_implied_by_fitted_params_avg":float(np.mean(np.exp(logpmf(n,np.zeros(len(n)),model,param)))),
           "conditional_on_same_original_fruit_and_flower_sample":True
        }
    return {
     "year":year,"morph_exact":sex,
     "genets":len(n),"flowers":int(n.sum()),"fruits":int(y.sum()),
     "fruit_positive_genets":int((y>0).sum()),
     "zero_fruit_genets":int((y==0).sum()),
     "total_fruit_rate":float(y.sum()/n.sum()),
     "model_fit":fits,
     "genet_level_crossvalidation":crossvalidate(n,y,rng),
    }

def main():
    source=data()
    adults=source["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(adults)==492
    rng=np.random.default_rng(SEED)
    out=[analyze(make(adults,year,sex),year,sex,rng)
         for year in YEARS for sex in SEXES]
    return {
      "status":"original_Daphne_jezoensis_zero_activation_vs_continuous_heterogeneity_v65",
      "source_article":"10.1002/ajb2.16239",
      "source_doctoral_thesis":"10.14943/doctoral.k14152",
      "original_author_zenodo":"10.5281/zenodo.8253921",
      "original_zip_sha256":ZIP_SHA,
      "source_original_2018_2020_adult_genets":len(adults),
      "seed":SEED,"fixed_model_designs":list(MODELS)+[FULL_DIAGNOSTIC],
      "sex_year":out,
      "previous_microscopy_published_thesis_stage_evidence":{
        "original_2016_hand_outcross_microscopy_pistils_H":55,
        "genets_H":6,
        "pistils_F":11,
        "genets_F":11,
        "thesis_reported_at_least_one_tube_reaching_ovule_more_than_90pct_each_morph":True,
        "pollen_tubes_arriving_at_ovule_does_not_establish_fertilization":True,
        "not_same_genets_years_flowers_as_2018_2020_natural_fruit_data":True,
        "stage_after_tube_arrival_and_before_mature_fruit_is_not_directly_measured":True,
        "original_resource_fertilization_leaf_cut_shade_nulls_previously_published":True,
        "no_new_discovery_of_early_pollen_tube_success_claimed":True
      },
      "limitations":{
        "zero_inflation_point_mass_is_model_component_not_observed_sterility":True,
        "beta_binomial_latent_heterogeneity_can_generate_most_zeros":True,
        "mixture_parameter_near_boundary_is_nonidentified":True,
        "multilocus_parentage_spatial_prior_not_used":True,
        "one_originating_source_per_two_repeated_years":True,
        "crossval_holds_adult_genets_within_sex_year_not_external_location":True,
        "2018_2020_overlap_not_between_year_external_independence":True,
        "true_seedling_survival_and_pollen_tube_to_fruit_transition_not_observed":True,
        "doctorate_2020_pollen_tube_micro_data_not_available_as_child_adult_join":True,
        "no_original_sensitive_raw_ID_SSR_alleles_or_coords_output":True,
        "focus_D_pseudomezereum_not_directly_observed":True,
      }}
if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--output",type=Path,default=OUT/"daphne_fruit_zero_state_model_comparison_v65.json")
    args=arg.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result=main()
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("DAPHNE_V65_HURDLE_CV",json.dumps(result,ensure_ascii=False,allow_nan=False),flush=True)
