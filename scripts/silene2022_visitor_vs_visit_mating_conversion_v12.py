#!/usr/bin/env python3
"""Re-analyse Barbot et al 2022 individual Silene dioica visitation vs fitness.

Primary original Dryad v251371 doi:10.5061/dryad.0rxwdbs4d,
FLODIS dataID_FloDis.txt pinned SHA. Tests whether **source 'nbVis'
(number of visiting pollinators, unidentified individual identities)** provides
incremental plant-level prediction over 'nbVisTot' (number of observed
visits) for genetic mating diversity, seed production or male siring.

Not a test of whether *fresh uniquely identified* bees versus same-bee
returns causally affect fitness: independent pollinator IDs and sequential
within-genet visiting histories are not present in this data table. Both
predictors are post-randomization observed outcomes, hence adjusting for
them is not a causal mediation analysis of fake flower treatment.

Use exact analytic leave-one-plant-out PRESS predictions for four nested
linear models, **only identical complete-case plants per endpoint**. Report
predictive differences with plant-bootstrap intervals of paired squared
errors; these are descriptive and not a pre-registered superiority test.
No original IDs or raw observation rows exported.
"""
from __future__ import annotations
import hashlib,io,json,math,urllib.request,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

OUT=Path("analysis_results_silene_visitors_vs_repeat_v12");OUT.mkdir(exist_ok=True)
VERSION=251371
URL=f"https://datadryad.org/api/v2/versions/{VERSION}/download"
ORIGINAL_SHA="1d245bb2dd6222bdaa05b0e23521eafd286cfa533b3826f8a5d8a65f715d3f9a"
RNG=np.random.default_rng(20261116);BOOT=7000

def val(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(np.integer,int)):return int(x)
    if isinstance(x,(np.floating,float)):return round(float(x),7) if np.isfinite(x) else None
    if isinstance(x,(list,tuple,np.ndarray)):return [val(z) for z in x]
    if isinstance(x,dict):return {str(k):val(z) for k,z in x.items()}
    return x

def get():
    req=urllib.request.Request(URL,headers={"Accept":"application/zip","User-Agent":"Open-Silene-individual-visitation-and-fitness-audit-v12/0.1"})
    with urllib.request.urlopen(req,timeout=90) as r:b=r.read(750000)
    if not zipfile.is_zipfile(io.BytesIO(b)):raise ValueError("unexpected Dryad archive format")
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        for info in z.infolist():
            if Path(info.filename).name=="dataID_FloDis.txt":
                raw=z.read(info)
                if hashlib.sha256(raw).hexdigest()!=ORIGINAL_SHA:
                    raise ValueError("original FLODIS individual source SHA mismatch")
                s=raw.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
                x=pd.read_csv(io.StringIO(s),sep=None,engine="python")
                x.columns=[str(k).strip() for k in x.columns]
                return x
    raise ValueError("missing original individual table")

def study():
    x=get()
    fields=["ID","nbVis","nbVisTot","mean_nbFlo","prodTOT","mean_nbOvTOT",
            "MedMRS","QttTot","MS_EchqueFem8"]
    if not set(fields)<=set(x):raise ValueError("FLODIS schema changed")
    if len(x)!=192 or x.ID.duplicated().any():raise ValueError("expected 192 unique plants")
    x["group"]=x.ID.astype(str).str.extract(r"^([A-Za-z]+)",expand=False)
    for name in fields[1:]:x[name]=pd.to_numeric(x[name],errors="coerce")
    x=x[x.group.isin(["FC","FT","MC","MT"])].copy()
    if len(x)!=172:raise ValueError("missing one original fake-display group")
    x["male"]=x.group.str.startswith("M").astype(int)
    x["fake"]=x.group.str.endswith("T").astype(int)
    both=x.dropna(subset=["nbVis","nbVisTot"])
    qa={
       "n_original_plants":192,"n_four_main_randomized_groups":len(x),
       "number_of_visiting_pollinators_is_NOT_tagged_individual_ID":True,
       "total_visits_is_observed_number_of_visits":True,
       "n_with_both_visitor_variables":len(both),
       "n_nbVis_greater_than_nbVisTot":int((both.nbVis>both.nbVisTot).sum()),
       "n_visitor_count_zero":int((both.nbVis==0).sum()),
       "n_total_visits_zero":int((both.nbVisTot==0).sum()),
       "n_visitors_gt0_but_total0":int(((both.nbVis>0)&(both.nbVisTot==0)).sum()),
       "n_source_visitors_missing":int(x.nbVis.isna().sum()),
       "n_source_visits_missing":int(x.nbVisTot.isna().sum()),
       "source_attached_repeated_individual_bee_IDs":False,
       "groups":{str(k):int(v) for k,v in x.group.value_counts().items()}
    }
    if (both[["nbVis","nbVisTot"]]<0).any().any():raise ValueError("negative visitation")
    if qa["n_nbVis_greater_than_nbVisTot"]>0:
        qa["repeat_count_decomposition_not_mathematically_valid"]=True
    else:
        qa["repeat_count_decomposition_not_mathematically_valid"]=False
        qa["n_returned_contacts_proxy_positive"]=int(((both.nbVisTot-both.nbVis)>0).sum())
    endpoints={
       "female_viable_seed_index":(0,"prodTOT","mean_nbOvTOT"),
       "female_genetic_mate_number":(0,"MS_EchqueFem8","mean_nbOvTOT"),
       "male_siring_MEMM_index":(1,"MedMRS","QttTot"),
       "male_genetic_mate_number":(1,"MS_EchqueFem8","QttTot")
    }
    models={}
    for title,(male,outcome,gamete) in endpoints.items():
        d=x[x.male==male].copy()
        d=d.dropna(subset=["nbVis","nbVisTot","mean_nbFlo",gamete,outcome])
        d=d[(d[["nbVis","nbVisTot","mean_nbFlo",gamete,outcome]]>=0).all(axis=1)].copy()
        if len(d)<60:raise ValueError("unexpected missingness "+title+" "+str(len(d)))
        d["outcome"]=np.log1p(d[outcome])
        d["lnF"]=np.log1p(d.mean_nbFlo)
        d["lnG"]=np.log1p(d[gamete])
        d["lnVisitors"]=np.log1p(d.nbVis)
        d["lnVisits"]=np.log1p(d.nbVisTot)
        designs={
           "true_flower_gamete_and_treatment":["fake","lnF","lnG"],
           "plus_visits":["fake","lnF","lnG","lnVisits"],
           "plus_visitor_number":["fake","lnF","lnG","lnVisitors"],
           "plus_both_visit_and_visitor_number":["fake","lnF","lnG","lnVisits","lnVisitors"]
        }
        folds={}
        res={"n_independent_plants":len(d),
             "n_observed_positive_pollinator_visitors":int((d.nbVis>0).sum()),
             "n_observed_positive_visits":int((d.nbVisTot>0).sum()),
             "visitors_and_total_visits_spearman":val(spearmanr(d.nbVis,d.nbVisTot).statistic),
             "original_source_outcome":outcome,
             "original_plant_level_selected_treatment_group_counts":{str(k):int(v) for k,v in d.group.value_counts().items()},
             "same_exact_complete_case_genets_across_all_models":True,
             "models":{}}
        for key,cols in designs.items():
            X=sm.add_constant(d[cols],has_constant="add").to_numpy(float)
            Y=d.outcome.to_numpy(float)
            fit=sm.OLS(Y,X).fit(cov_type="HC3")
            # Exact LOO error for ordinary least squares with same fixed X
            h=np.sum(X*(X@np.linalg.pinv(X.T@X)),axis=1)
            if np.any(h>=.97):res["high_leverage_warning"]=int((h>=.97).sum())
            loo=fit.resid/(1-h)
            folds[key]=np.square(loo)
            named={}
            for idx,name in enumerate(["intercept"]+cols):
                ci=fit.conf_int()[idx]
                named[name]={"beta":val(fit.params[idx]),
                    "95_HC3":[val(ci[0]),val(ci[1])],"p":val(fit.pvalues[idx])}
            res["models"][key]={
               "n_complete_case_plants":len(d),"terms":named,
               "in_sample_R2":val(fit.rsquared),
               "LOO_plant_MSE":val(np.mean(loo**2)),
               "LOO_plant_RMSE":val(np.sqrt(np.mean(loo**2))),
               "condition_number":val(np.linalg.cond(X)),
               "additional_covariates_postrandomization_not_causal":True}
        main="plus_both_visit_and_visitor_number"
        for compared in ("plus_visits","plus_visitor_number","true_flower_gamete_and_treatment"):
            a=folds[compared]-folds[main]
            bs=a[RNG.integers(0,len(a),size=(BOOT,len(a)))].mean(axis=1)
            res["models"][main]["LOO_MSE_improvement_vs_"+compared]={
              "MSE_compared_minus_joint":val(a.mean()),
              "paired_plant_bootstrap95":[val(t) for t in np.quantile(bs,[.025,.975])],
              "model_selection_exploratory_not_preactivated":True}
        models[title]=res
    return {"status":"success_original_2022_Silene_same_plant_visitation_fitness_audit",
            "doi":"10.1098/rspb.2022.1987",
            "dryad":"10.5061/dryad.0rxwdbs4d",
            "source_sha256":ORIGINAL_SHA,
            "source_exposure_definitions":qa,
            "comparisons_by_sex_and_endpoint":models,
            "new_arrivals_cannot_be_identified_from_unmarked_visitors":True,
            "study_design_warning":"Observational associations conditional on visual-only randomized assignment; no causal mediation and no joint actual effective pollen deposition",
            "no_original_individual_ids_or_raw_records_exported":True}

if __name__=="__main__":
    r=study()
    path=OUT/"silene2022_visitors_repeat_proxy_vs_genetic_fitness_v12.json"
    path.write_text(json.dumps(val(r),ensure_ascii=False,indent=2)+"\n")
    print("SILENE_VISIT_REVISIT_PRIMARY",json.dumps(val(r),ensure_ascii=False)[:40000],flush=True)
