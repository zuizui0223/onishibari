#!/usr/bin/env python3
"""Test 2018/2020 repeatability of residual assigned father success in Daphne jezoensis.

This is a NEW secondary diagnostic of the same Shibata & Kudo (2023) AJB
source used by v28/v30/v38 (Zenodo 8253921). Not an additional study.

Primary null: mothers independently allocate their fixed observed assigned
offspring counts to eligible H sires with probabilities proportional to
exp(-distance/10m)*same-year source natural-flower-count^fitted_beta.
For genets H-coded and eligible in BOTH years, correlate excess siring
log((observed+.5)/(model expected+.5)) across 2018 and 2020.
10k conditional Monte Carlo nulls include correlated competing fathers
within maternal family and heterogeneous sire opportunity; they condition
on fitted beta, maternal sampling frames, father eligibility and positions.
This does NOT identify stable heritable quality or true total reproductive
fitness. Sensitive individual IDs and coordinates never leave memory.

Also report same-genet H/H natural-fruit recurrence as a descriptive
crossyear contrast with a fixed-marginal Fisher exact test (small n).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import softmax
from scipy.stats import rankdata, spearmanr, fisher_exact
from daphne_jezoensis_2018_2020_HF_fathers_v28 import data, stringify, ZIP_SHA
from daphne_jezoensis_2018_2020_flower_to_paternity_choice_v30 import build, estimate

OUT=Path("analysis_results_daphne_stable_sire_v53")
LAMBDA=10.0
N_SIM=10000
SEED=20261122

def donor_keys(adult,year):
    sex=f"Sex_{year}.g"
    flowers=f"NaturalFlowerNo_{year}.g"
    d=adult[["GenetID2","East.g","North.g",sex,flowers]].copy()
    d["key"]=stringify(d.GenetID2)
    d["x"]=pd.to_numeric(d["East.g"],errors="coerce")
    d["y"]=pd.to_numeric(d["North.g"],errors="coerce")
    d["flower"]=pd.to_numeric(d[flowers],errors="coerce")
    d=d[(d[sex].astype(str)=="H")&d.flower.gt(0)&d.x.notna()&d.y.notna()]
    assert not d.key.duplicated().any()
    return d.key.tolist()

def siring(q):
    b=float(estimate(q,LAMBDA)["beta_hat"])
    lf=q["logflowers"]
    n=q["counts"].sum(axis=1)
    prob=softmax(-q["distance"]/LAMBDA+b*lf[None,:],axis=1)
    expected=np.sum(n[:,None]*prob,axis=0)
    observed=q["counts"].sum(axis=0)
    assert len(prob)==len(n) and np.isclose(expected.sum(),observed.sum())
    return dict(beta=b,n=n,prob=prob,expected=expected,observed=observed,
                residual=np.log((observed+.5)/(expected+.5)))

def simulate(parent, rng):
    rows=np.zeros((N_SIM,len(parent["observed"])),dtype=np.int32)
    for n,prob in zip(parent["n"],parent["prob"]):
        rows+=rng.multinomial(int(n),prob,size=N_SIM).astype(np.int32)
    assert np.all(rows.sum(axis=1)==int(parent["observed"].sum()))
    return rows

def corr_rows(a,b):
    # Stable tie-aware per-row Spearman using scipy column-rank transform.
    ra=rankdata(a,axis=1,method="average")
    rb=rankdata(b,axis=1,method="average")
    ra-=ra.mean(axis=1,keepdims=True)
    rb-=rb.mean(axis=1,keepdims=True)
    den=np.linalg.norm(ra,axis=1)*np.linalg.norm(rb,axis=1)
    return np.divide((ra*rb).sum(axis=1),den,out=np.full(len(den),np.nan),where=den>0)

def stat(v):
    return float(v) if np.isfinite(v) else None

def fruit_recurrence(adult):
    keep=pd.Series(True,index=adult.index)
    for y in [2018,2020]:
        sx=f"Sex_{y}.g"
        fl=f"NaturalFlowerNo_{y}.g"
        fr=f"NaturalFruitNo_{y}.g"
        keep &= adult[sx].astype(str).eq("H")
        keep &= pd.to_numeric(adult[fl],errors="coerce").gt(0)
        v=pd.to_numeric(adult[fr],errors="coerce")
        keep &= v.notna()&v.ge(0)&v.le(pd.to_numeric(adult[fl],errors="coerce"))
    d=adult.loc[keep].copy()
    if len(d)<10:raise AssertionError("Too few repeated strict-H source adults")
    a=pd.to_numeric(d["NaturalFruitNo_2018.g"])>0
    b=pd.to_numeric(d["NaturalFruitNo_2020.g"])>0
    tab=np.array([[int((~a&~b).sum()),int((~a&b).sum())],
                  [int((a&~b).sum()),int((a&b).sum())]])
    assert tab.sum()==len(d)
    OR,p=fisher_exact(tab)
    return {"n_eligible_both_years":len(d),"2018_fruit_positive":int(a.sum()),
      "2020_fruit_positive":int(b.sum()),
      "same_genet_fruiting_both_years":int((a&b).sum()),
      "fruit_binary_2x2_2018_rows_2020_cols":tab.tolist(),
      "independent_fixed_margins_expected_repeated":float(a.sum()*b.sum()/len(d)),
      "fisher_two_sided_p":float(p),"fisher_odds_ratio":stat(OR),
      "limitations":"conditional on flowering H both years; no 2019, true viable seeds unknown, plant size/fruit probability not adjusted"}

def main():
    src=data()
    ad=src["GenetID_geno_phenotype_1820Nopadult_211109.csv"]
    assert len(ad)==492 and ad.GenetID2.nunique()==492
    y={year:build(year,src,include_HF=False) for year in (2018,2020)}
    stats={yr:siring(q) for yr,q in y.items()}
    keys={yr:donor_keys(ad,yr) for yr in (2018,2020)}
    for yr in (2018,2020):
        assert len(keys[yr])==len(stats[yr]["observed"])
        assert len(keys[yr])==(105 if yr==2018 else 71)
        assert (y[yr]["counts"].sum(axis=1)>0).all()
    # Filter strictly H in both years, in the exact order of original adult table.
    common=[key for key in keys[2018] if key in set(keys[2020])]
    idx={yr:np.asarray([keys[yr].index(k) for k in common],int) for yr in (2018,2020)}
    assert len(common)>=35 and len(common)<=71
    s1=stats[2018];s2=stats[2020]
    i1=idx[2018];i2=idx[2020]
    def subset(name,ids1=i1,ids2=i2):
        return getattr(spearmanr(s1[name][ids1],s2[name][ids2]),"statistic")
    primary=subset("residual")
    o1=s1["observed"][i1];o2=s2["observed"][i2]
    exp1=s1["expected"][i1];exp2=s2["expected"][i2]
    overlap=int(((o1>0)&(o2>0)).sum())
    rng=np.random.default_rng(SEED)
    sim1=simulate(s1,rng)[:,i1]
    sim2=simulate(s2,rng)[:,i2]
    r1=np.log((sim1+.5)/(exp1+.5))
    r2=np.log((sim2+.5)/(exp2+.5))
    null=corr_rows(r1,r2)
    valid=null[np.isfinite(null)]
    assert len(valid)>N_SIM*.99
    p_up=float((1+np.sum(valid>=primary))/(1+len(valid)))
    p_two=float((1+np.sum(np.abs(valid)>=abs(primary)))/(1+len(valid)))
    overlap_null=((sim1>0)&(sim2>0)).sum(axis=1)
    p_overlap_up=float((1+np.sum(overlap_null>=overlap))/(N_SIM+1))
    # One-donor exclusion affects the biological question (a post-hoc
    # influence probe), not the prespecified Monte Carlo p-value.
    loo=[]
    for a in range(len(common)):
        keep=np.arange(len(common))!=a
        loo.append(float(spearmanr(s1["residual"][i1[keep]],s2["residual"][i2[keep]]).statistic))
    def support_threshold(v):
        mask=(exp1>=v)&(exp2>=v)
        return {"n":int(mask.sum()),"spearman":stat(spearmanr(s1["residual"][i1[mask]],s2["residual"][i2[mask]]).statistic)
                if mask.sum()>=6 else None}
    fruit=fruit_recurrence(ad)
    return {
      "status":"original_Daphne_same_genet_multiyear_sire_residual_repeatability_v53",
      "source_doi":"10.1002/ajb2.16239","zenodo_record":"10.5281/zenodo.8253921",
      "original_zip_sha256":ZIP_SHA,"adult_genets":len(ad),
      "years":[{"year":yr,"candidate_strict_H":len(keys[yr]),
               "offspring_assigned_to_eligible_H":int(stats[yr]["observed"].sum()),
               "sampled_mother_genets":len(stats[yr]["n"]),
               "distance_conditional_beta":float(stats[yr]["beta"])}
                for yr in (2018,2020)],
      "same_strict_H_genet_in_both_years":len(common),
      "repeatability":{
        "primary_score":"Spearman correlation of per-donor log((assigned offspring +0.5)/(model-expected +0.5)); same eligible H sire genets 2018 and 2020",
        "primary_spearman":float(primary),
        "conditional_independent_maternal_multinomial_null_mean":float(valid.mean()),
        "null_95percent":np.quantile(valid,[.025,.975]).tolist(),
        "upper_tail_p":p_up,"two_sided_p":p_two,
        "null_draws":len(valid),"seed":SEED,
        "raw_assigned_sire_count_spearman":stat(spearmanr(o1,o2).statistic),
        "expected_sire_count_spearman":stat(spearmanr(exp1,exp2).statistic),
        "overlapping_sires_with_at_least_one_assigned_both_years":overlap,
        "null_overlap_mean":float(overlap_null.mean()),
        "null_overlap_95percent":np.quantile(overlap_null,[.025,.975]).tolist(),
        "null_overlap_upper_tail_p":p_overlap_up,
        "leave_one_genet_out_primary_range":[float(min(loo)),float(max(loo))],
        "expected_count_support_both_at_least_0p5":support_threshold(.5),
        "expected_count_support_both_at_least_1":support_threshold(1),
        "expected_count_support_both_at_least_2":support_threshold(2),
      },
      "same_H_fruit_recurrence_descriptive":fruit,
      "strict_limits":{
        "same_population_two_years":True,
        "source_95percent_SSR_father_assignment_not_whole_population_success":True,
        "mother_offspring_counts_fixed_within_year":True,
        "mean_model_conditions_on_fitted_year_beta_and_10m_kernel":True,
        "sire_flower_count_observational_and_spatial_opportunity_model_misspec_possible":True,
        "not_stable_genetic_quality_or_gene_flow_causal_identification":True,
        "original_parent_IDs_and_coordinates_never_output":True,
        "target_Daphne_pseudomezereum_unmeasured":True}}
if __name__=="__main__":
    OUT.mkdir(exist_ok=True,parents=True)
    x=main()
    path=OUT/"daphne_same_H_sire_residual_repeatability_2018_2020_v53.json"
    path.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("DAPHNE_SIRE_REPEATABILITY_V53",json.dumps(x,ensure_ascii=False,allow_nan=False),flush=True)
