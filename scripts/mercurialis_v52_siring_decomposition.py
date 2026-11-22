#!/usr/bin/env python3
"""Original-record v0.52 paternal-self vs paternal-outcross allocation audit.
Data: Chen & Pannell (2025) New Phytol DOI 10.1111/nph.70356,
Zenodo fixed 10.5281/zenodo.15644955, 181 parents x 2 assumed delta.
All output is aggregates; no IDs, parentage rows or genotype data emitted.
Male and female fitness in original ID=1 scenario are FORMULA-SCENARIO
outcomes, not measured offspring survival or independently observed children.
"""
import argparse, hashlib, io, json, urllib.request, zipfile, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

URL="https://zenodo.org/api/records/15644955/files/Geito_NewPhyto_v2.zip/content"
SHA="a18b7d296110be1e9817d6e1a3994365a4e7ae0aab50435e7e6883e1f3de4042"
SEED=70357
DRAWS=1599
OUTCOMES=["outcross_siring","self_siring","all_siring","outcross_maternal"]
def archive(path):
    if path:
        b=path.read_bytes()
    else:
        q=urllib.request.Request(URL,headers={"User-Agent":"Source-backed-fitness-components-0.52"})
        b=urllib.request.urlopen(q,timeout=90).read()
    if hashlib.sha256(b).hexdigest()!=SHA: raise ValueError("Pinned source archive checksum mismatch")
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        a=pd.read_csv(io.BytesIO(z.read("Geito_NewPhyto_v2/3G.csv")))
        f=pd.read_csv(io.BytesIO(z.read("Geito_NewPhyto_v2/3G_ID01.csv")))
    return a,f

def prep(a,f):
    assert len(a)==181 and len(f)==362
    assert not a.duplicated(["Population","ID"]).any()
    f=f.copy()
    f["Population"]=f.Population.str.replace("3G_P","Population ",regex=False)
    assert not f.duplicated(["Population","ID","Inbreeding"]).any()
    assert set(f.Inbreeding)=={"ID=0","ID=1"}
    p=f.pivot(index=["Population","ID"],columns="Inbreeding",
              values=["Male_fitness","Female_fitness","Total_fitness"])
    p.columns=["%s_%s"%(x,y.replace("ID=","d")) for x,y in p.columns]
    p=p.reset_index()
    m=a.merge(p,on=["Population","ID"],how="inner",validate="one_to_one")
    assert len(m)==181
    for d in (0,1):
        x=m[f"Total_fitness_d{d}"]
        y=m[f"Male_fitness_d{d}"]+m[f"Female_fitness_d{d}"]
        q=(x-y).abs().dropna()
        assert len(q)>=180 and q.max()<2e-5, f"sex-fitness sum mismatch d{d} {q.max()}"
    m["outcross_siring"]=m.Male_fitness_d1
    m["all_siring"]=m.Male_fitness_d0
    m["self_siring"]=m.Male_fitness_d0-m.Male_fitness_d1
    m["outcross_maternal"]=m.Female_fitness_d1
    m["self_maternal"]=m.Female_fitness_d0-m.Female_fitness_d1
    q=(m.self_siring-m.self_maternal).abs().dropna()
    consistency={"paired_valid":int(len(q)),"self_male_female_max_abs_difference":float(q.max()),
                 "self_male_female_abs_difference_median":float(q.median())}
    # The original construction attributes each selfed seed equally
    # to both sex functions. Distinct partition is only accepted if
    # formula identities and non-negativity are source-verified.
    assert len(q)>=180 and q.max()<2e-5, f"paternal and maternal self components disagree {q.max()}"
    for name in ["self_siring","outcross_siring","outcross_maternal"]:
        v=m[name].dropna()
        assert v.min()>-1e-5, (name,v.min())
        m[name]=m[name].clip(lower=0)
    s=m.loc[m.Selfing_rate.notna() & m.Female_fitness_d0.gt(0)].copy()
    err=((s.Female_fitness_d0-s.Female_fitness_d1)/s.Female_fitness_d0-s.Selfing_rate).abs()
    assert len(err)>=172 and err.max()<1e-5
    consistency["selfing_rate_discount_max_error"]=float(err.max())
    return m,consistency

def Xcols(df,center=None):
    pop=pd.get_dummies(df.Population)
    assert set(pop)=={"Population 1","Population 2","Population 3"}
    x=np.column_stack([pop[["Population 1","Population 2","Population 3"]].to_numpy(float),
         np.log2(df.Male_flower_total.to_numpy(float)),
         np.log2(df.Female_flower_total.to_numpy(float)),
         np.log2(df.Biomass_whole.to_numpy(float))])
    if center is None: center=x[:,3:].mean(axis=0)
    x[:,3:]-=center
    return x,center

def ppml(X,y):
    # Nonnegative fractional shares from the authors' paternity expansion:
    # Poisson pseudo-ML with an HC3 plant sandwich, NOT a literal Poisson seed model.
    assert y.min()>=0 and (y>0).any()
    popid=np.argmax(X[:,:3],axis=1)
    init=np.zeros(X.shape[1],float)
    for j in range(3):
        init[j]=np.log(max(.001,y[popid==j].mean()))
    def objective(b):
        eta=np.clip(X@b,-25,25)
        mu=np.exp(eta)
        return float(np.sum(mu-y*eta)),X.T@(mu-y)
    opt=minimize(objective,init,jac=True,method="L-BFGS-B",
                 options={"maxiter":1000,"ftol":1e-12,"gtol":1e-6})
    if not opt.success and np.max(abs(opt.jac))>1e-3:
        raise RuntimeError(f"PPML: {opt.message}, gradient max {np.max(abs(opt.jac)):.6g}")
    b=opt.x
    mu=np.exp(X@b)
    bread=np.linalg.inv(X.T@(mu[:,None]*X))
    leverage=mu*np.einsum("ij,jk,ik->i",X,bread,X)
    score=X*((y-mu)/np.maximum(.05,1-leverage))[:,None]
    cov=bread@(score.T@score)@bread
    se=np.sqrt(np.diag(cov))
    return b,se,float(((y-mu)**2/np.maximum(mu,1e-9)).sum()/(len(y)-len(b)))

def estimate(df):
    X,cent=Xcols(df)
    n=len(df)
    results={}
    for outcome in OUTCOMES:
        y=df[outcome].to_numpy(float)
        b,se,disp=ppml(X,y)
        # β is the elasticity of expected component per doubling of male
        # flowers, conditional on female flowers, final biomass, and garden.
        point=np.exp(b[3])
        interval=np.exp([b[3]-1.959963984540054*se[3],b[3]+1.959963984540054*se[3]])
        results[outcome]={"per_male_flower_doubling_ratio":float(point),
              "hc3_95":interval.tolist(),"ppml_log2_male_b":float(b[3]),
              "ppml_male_hc3_se":float(se[3]),"pearson_dispersion":disp,
              "mean_expanded_source_success":float(y.mean()),"nonzero":int((y>0).sum())}
    # Resampling at the father plant level *within* fixed arrays excludes
    # genetic parentage uncertainty, interference and new population effects.
    rng=np.random.default_rng(SEED)
    group=[np.flatnonzero(df.Population.to_numpy()==x) for x in
           ["Population 1","Population 2","Population 3"]]
    B={outcome:[] for outcome in OUTCOMES}
    for i in range(DRAWS):
        ix=np.concatenate([rng.choice(g,len(g),replace=True) for g in group])
        Xi=X[ix]
        for outcome in OUTCOMES:
            y=df[outcome].to_numpy(float)[ix]
            if not (y>0).any():continue
            bb,_,_=ppml(Xi,y)
            B[outcome].append(float(np.exp(bb[3])))
    for outcome in OUTCOMES:
        results[outcome]["resample_success"]=len(B[outcome])
        results[outcome]["within_array_father_bootstrap95"]=np.quantile(B[outcome],[.025,.975]).tolist()
    return results

def main(a,f):
    m,consistency=prep(a,f)
    d=m.loc[m.Male_flower_total.gt(0)&m.Female_flower_total.gt(0)&m.Biomass_whole.notna()].copy()
    d=d.dropna(subset=OUTCOMES).reset_index(drop=True)
    assert len(d)>=160
    result={"status":"original_data_secondary_v52","source_doi":"10.1111/nph.70356",
        "zenodo_version":"10.5281/zenodo.15644955","archive_sha256":SHA,
        "source_parents":len(m),"paired_scenario_rows":len(f),
        "formula_identity_audit":consistency,
        "complete_bisexual_size_subset_n":len(d),
        "n_both_source_sex_function_outputs":int(m.Female_fitness_d1.notna().sum()),
        "estimation":"PPML mean ratios for expanded genetically estimated successes, male+female flower counts and aboveground biomass log2, three population fixed intercepts; HC3 father robust and 1599 within-population plant bootstraps",
        "outcome_definitions":{
            "outcross_siring":"Male_fitness under assumed delta=1: paternity shares expanded by each maternal seed total, selfed progeny discounted to zero",
            "self_siring":"Male_fitness(delta=0)-Male_fitness(delta=1), same seed selfing contribution also accounted through female fitness",
            "all_siring":"Male_fitness under assumed delta=0",
            "outcross_maternal":"Female_fitness under assumed delta=1"},
        "effects":estimate(d),
        "limitations":["observational male flower allocation not experimentally randomized",
        "paired delta endpoints are arithmetic scenarios, not measured viable offspring",
        "paternity-attributed paternal totals are expanded from a sample, not raw count",
        "three garden arrays not independent species and not enough for population-level cluster inference",
        "fixed network, related genets and paternity assignment errors not in plant HC3 or plant bootstrap",
        "male flower number not floral dry biomass"] }
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    a,f=archive(args.archive)
    result=main(a,f)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("FITNESS_DECOMPOSITION_V52",json.dumps(result,ensure_ascii=False,allow_nan=False))
