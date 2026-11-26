#!/usr/bin/env python3
"""Source-verified audit of per-flower biological seed upper caps and attrition in
Pulsatilla alpina 2022 experimental cohort. Does not report plant IDs.
"""
import tempfile,json
from pathlib import Path
import pandas as pd,numpy as np
from pulsatilla_v61_selection_tipping import source
from pulsatilla_selfing_refit_v51 import ARMS

def output():
    p,t,cb,_=source(Path(tempfile.gettempdir())/"onishibari_pulsatilla_v62_pinned")
    statuses={a:[] for a in ARMS}
    result={"source_census_2022_columns":list(p.columns),
            "source_typed_2024_columns":list(t.columns),
            "source_outcomes_2025_columns":list(cb.columns),
            "groups":{}}
    for arm in ARMS:
        x=p[p.arm==arm].copy()
        typed_ids=set(t.loc[t.arm==arm,"Ind_ID"])
        matched=set(cb.loc[(cb.arm==arm)&(cb.Floral_gender=="H"),"Ind_ID"])
        m1=x.Ind_ID.isin(typed_ids)
        m2=x.Ind_ID.isin(matched)
        groups={}
        for label,mask in (
         ("genotyped",m1),
         ("outcome_known_untyped",m2&~m1),
         ("no_matched_mature_outcome",~m2)):
            y=x.loc[mask]
            features={}
            for col in x.columns:
                if col in ["Ind_ID","arm","Stamen_removal","Floral_gender"] or "coord" in col.lower():continue
                z=pd.to_numeric(y[col],errors="coerce")
                full=pd.to_numeric(x[col],errors="coerce")
                if full.notna().sum()<.6*len(x):continue
                features[col]={
                  "count":int(z.notna().sum()),
                  "missing":int(z.isna().sum()),
                  "min":float(z.min()) if z.notna().any() else None,
                  "median":float(z.median()) if z.notna().any() else None,
                  "mean":float(z.mean()) if z.notna().any() else None,
                  "max":float(z.max()) if z.notna().any() else None,
                  "total":float(z.sum()) if z.notna().any() else None
                }
            groups[label]={"flowers":len(y),"fields":features}
        result["groups"][arm]=groups
    return result
if __name__=="__main__":
    print("PULSATILLA_V62_SOURCE_FIELD_AUDIT",json.dumps(output(),allow_nan=False,ensure_ascii=False),flush=True)
