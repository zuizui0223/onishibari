#!/usr/bin/env python3
"""v66: annual persistence of latent fruiting propensity, conditional on flower counts.
Same original Daphne jezoensis study (Shibata & Kudo 2023, DOI
10.1002/ajb2.16239), not independent replication or D. pseudomezereum.
Use exact source adult IDs only in memory. No original individual data exported.

Models for same-morph genets in both years:
D: independent annual beta-binomial latent probabilities with common beta prior
P: fully persistent probability p_i ~ Beta(a,b) in both years
M: mixture (1-w) D + w P, with w between 0 and 1.
A shared p_i gives beta-binomial posterior prediction of second-year count
conditional on previous-year n,y. This is a diagnostic of statistical persistence,
NOT an anatomical fertility test or causal model. Common year distribution is
motivated by v65 for H; F is comparison/sensitivity, not assumed stationary.
Primary likelihood is joint count PMF on each same-genet pair. Compare leave-one-
genet-out predictive gain and 2018<->2020 new-genet-only transfer from v65.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln,gammaln,expit,logit
from scipy.stats import spearmanr
from daphne_v65_zero_state_predictive_comparison import pinned_adults_only,fit,logpmf,MODELS
from daphne_v64_natural_fruit_gate import make

OUT=Path("analysis_results_daphne_temporal_persistence_v66")
SEED=20261126
BOOT=5000
SIM=10000
def comb(n,y):
    return gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)
def lp_bb(n,y,a,b):
    return comb(n,y)+betaln(a+y,b+n-y)-betaln(a,b)
def components(n1,y1,n2,y2,a,b):
    d=lp_bb(n1,y1,a,b)+lp_bb(n2,y2,a,b)
    p=comb(n1,y1)+comb(n2,y2)+betaln(a+y1+y2,b+n1+n2-y1-y2)-betaln(a,b)
    return d,p
def joint(n1,y1,n2,y2,theta,kind):
    p=expit(theta[0]);k=np.exp(theta[1]);a=p*k;b=(1-p)*k
    d,s=components(n1,y1,n2,y2,a,b)
    if kind=="independent":return d
    if kind=="persistent":return s
    w=expit(theta[2])
    return np.logaddexp(np.log1p(-w)+d,np.log(w)+s)
def fit_joint(n1,y1,n2,y2,kind):
    rate=np.clip((y1.sum()+y2.sum())/(n1.sum()+n2.sum()),1e-7,.99)
    starts=[[logit(rate),np.log(k)] for k in (.3,3,30,200)]
    bounds=[(-16,8),(-9,11)]
    if kind=="mixture":
        starts=[x+[logit(w)] for x in starts for w in (.1,.5,.9)]
        bounds+=[(-14,14)]
    best=None
    for x in starts:
        z=minimize(lambda t:-joint(n1,y1,n2,y2,t,kind).sum(),x,method="L-BFGS-B",bounds=bounds,
                   options={"maxiter":300,"ftol":1e-10})
        if np.isfinite(z.fun) and (best is None or z.fun<best.fun):best=z
    if best is None:raise RuntimeError("no finite fit")
    return best.x,float(-best.fun),bool(best.success)
def p_summary(x,kind):
    p=expit(x[0]);k=np.exp(x[1])
    return {"mean":float(p),"concentration":float(k),"weight_persistent":float(expit(x[2])) if kind=="mixture" else float(kind=="persistent")}
def sim_overlap(n1,n2,a,b,rng):
    # Conditional on observed flower opportunities and fitted marginal distribution,
    # simulate whether each genet fruits at least once. Shared vs independent.
    # No forced same total fruits; includes count stochasticity and latent variation.
    p1=rng.beta(a,b,size=(SIM,len(n1)))
    p2=rng.beta(a,b,size=(SIM,len(n2)))
    z1=rng.binomial(n1,p1)>0
    z2d=rng.binomial(n2,p2)>0
    z2p=rng.binomial(n2,p1)>0
    return {k:np.sum(z1&z,axis=1) for k,z in [("independent",z2d),("persistent",z2p)]}
def analyze(cohorts,sex,rng):
    g1,g2=[cohorts[(y,sex)] for y in (2018,2020)]
    i1={k:i for i,k in enumerate(g1["ids"])};i2={k:i for i,k in enumerate(g2["ids"])}
    ids=sorted(set(i1)&set(i2));ix=np.array([i1[k] for k in ids]);jx=np.array([i2[k] for k in ids])
    n1,y1=g1["n"][ix],g1["y"][ix];n2,y2=g2["n"][jx],g2["y"][jx]
    N=len(ids)
    assert N>10
    if sex=="H":assert (N,int((y1>0).sum()),int((y2>0).sum()),int(((y1>0)&(y2>0)).sum()))==(45,6,6,1)
    kinds=("independent","persistent","mixture")
    full={}
    for kind in kinds:
        th,ll,ok=fit_joint(n1,y1,n2,y2,kind)
        full[kind]={"parameters":p_summary(th,kind),"joint_loglikelihood":ll,"AIC":2*len(th)-2*ll,"optimizer_success":ok}
    # Paired leave-one-genet-out, entire two-year record excluded from estimation.
    scores={k:np.zeros(N) for k in kinds}
    conditional={k:np.zeros((2,N)) for k in kinds}
    opt={k:0 for k in kinds}
    for i in range(N):
        keep=np.arange(N)!=i
        for k in kinds:
            th,ll,ok=fit_joint(n1[keep],y1[keep],n2[keep],y2[keep],k)
            opt[k]+=int(ok)
            score=joint(n1[i:i+1],y1[i:i+1],n2[i:i+1],y2[i:i+1],th,k)[0]
            p=expit(th[0]);conc=np.exp(th[1]);a=p*conc;b=(1-p)*conc
            scores[k][i]=score
            conditional[k][0,i]=score-lp_bb(n1[i],y1[i],a,b)
            conditional[k][1,i]=score-lp_bb(n2[i],y2[i],a,b)
    draws=rng.integers(N,size=(BOOT,N))
    def compare(z):
        return {"mean":float(z.mean()),"genet_resample95":np.quantile(z[draws].mean(axis=1),[.025,.975]).tolist(),
                "positive_genets":int((z>0).sum())}
    comparisons={}
    for k in ("persistent","mixture"):
        comparisons[k+"_minus_independent"]={
          "joint":compare(scores[k]-scores["independent"]),
          "conditional_2018_to_2020":compare(conditional[k][0]-conditional["independent"][0]),
          "conditional_2020_to_2018":compare(conditional[k][1]-conditional["independent"][1])}
    par=full["independent"]["parameters"]
    sims=sim_overlap(n1,n2,par["mean"]*par["concentration"],(1-par["mean"])*par["concentration"],rng)
    overlap=int(((y1>0)&(y2>0)).sum())
    simsum={k:{"mean":float(z.mean()),"interval95":np.quantile(z,[.025,.975]).tolist(),"tail_le_observed":float((1+(z<=overlap).sum())/(SIM+1))} for k,z in sims.items()}
    transfer=[]
    for yr1,yr2 in ((2018,2020),(2020,2018)):
        a=cohorts[(yr1,sex)];b=cohorts[(yr2,sex)]
        unseen=np.array([x not in set(a["ids"]) for x in b["ids"]])
        assert unseen.any()
        row={"train_year":yr1,"test_year":yr2,"unseen_genets":int(unseen.sum()),
             "unseen_fruit_positive":int((b["y"][unseen]>0).sum()),"models":{}}
        for model in MODELS:
            th,_,ok=fit(a["n"],a["y"],model)
            n=b["n"][unseen];y=b["y"][unseen]
            row["models"][model]={"mean_logscore":float(logpmf(n,y,model,th).mean()),
                "predicted_zero_fraction":float(np.exp(logpmf(n,np.zeros(len(n),int),model,th)).mean()),
                "observed_zero_fraction":float((y==0).mean()),"optimizer_success":ok}
        transfer.append(row)
    return {"sex":sex,"paired_genets":N,"flower_totals":[int(n1.sum()),int(n2.sum())],
      "fruit_totals":[int(y1.sum()),int(y2.sum())],"positive_genets":[int((y1>0).sum()),int((y2>0).sum())],
      "positive_both":overlap,"flower_count_spearman":float(spearmanr(n1,n2).statistic),
      "fruit_rate_spearman_descriptive":float(spearmanr(y1/n1,y2/n2).statistic),
      "joint_models":full,"leave_one_genet_out":{
        k:{"joint_mean_logscore":float(scores[k].mean()),"conditional_mean_logscore_18_to_20":float(conditional[k][0].mean()),"conditional_mean_logscore_20_to_18":float(conditional[k][1].mean()),"optimizer_success_count":opt[k]} for k in kinds},
      "comparisons":comparisons,"simulated_overlap":simsum,"new_genet_only_crossyear_transfer":transfer}
def main():
    df=pinned_adults_only()
    cohorts={(y,s):make(df,y,s) for y in (2018,2020) for s in ("H","F")}
    rng=np.random.default_rng(SEED)
    return {"source_doi":"10.1002/ajb2.16239","source_record":"10.5281/zenodo.8253921",
      "status":"v66_annual_persistence_paired_genet_models",
      "seed":SEED,"bootstrap_draws":BOOT,"overlap_simulations":SIM,
      "by_sex":[analyze(cohorts,s,rng) for s in ("H","F")],
      "limitations":["common beta prior across years; year-specific environment can differ",
       "mixture weight is statistical persistence component, not actual proportion permanently fertile",
       "shared population and survival/flowering selection, only two observed years",
       "within-cohort leave-one-genet-out not independent population",
       "overlap simulation plug-in uncertainty not propagated",
       "genet resampling intervals descriptive conditional on fitted CV",
       "natural mature fruit counts not fertilization, recruitment or D. pseudomezereum data"]}
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,default=OUT/"daphne_temporal_persistence_v66.json")
    args=p.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    v=main();args.output.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
    print("DAPHNE_V66_TEMPORAL",json.dumps(v,allow_nan=False),flush=True)
