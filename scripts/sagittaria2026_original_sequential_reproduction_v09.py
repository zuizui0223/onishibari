#!/usr/bin/env python3
"""Original 2026 Sagittaria trifolia within-season reproductive compensation audit.

Tang & Dai (2026) Ecol Evol DOI 10.1002/ece3.73955, openly published
Figshare DOI 10.6084/m9.figshare.27610539. Source uses clone genotypes across
three *site-confounded* environmental treatments (net/shaded/open).
The net treatment is NOT a factorial randomized pollen-vs-light contrast:
net and shade are separate sites; cover changes microclimate too. No inference
of universal compensation from between-site effects.

Source plant-level 'TFS' and inflorescence-level fruits recorded on the same
plants. Analyze ordered inflorescences 1 vs >=2 within *reproductive plants*,
then whole plant seed counts, keeping missing vs true zeros separate.
Genotype-level bootstrap within each treatment; all outputs aggregated.
"""
from __future__ import annotations
import hashlib,io,json,math,urllib.request
from collections import Counter
from pathlib import Path
import numpy as np,pandas as pd
OUT=Path("analysis_results_sagittaria2026_sequential_fitness_v09");OUT.mkdir(exist_ok=True)
BASE="https://ndownloader.figshare.com/files/"
FILES={
 "inflorescences.csv":(50204946,"72ff4762ef2b90455bdbb7c39be3c047ef7182e912964b92ae3b5c81f70e71fa"),
 "total.csv":(63676992,"0b24262abe2dfe5520c5f5dd3141eb3c00bbf1acae29716c8ff871fbe5739551"),
 "plant information.csv":(50204961,"962b712015b62b409c4f42f7306d87414c271982a517675cef5a9a605e7baf11"),
}
RNG=np.random.default_rng(20261116);N_BOOT=5000;N_PERM=6000

def original():
    tables={}
    for name,(fid,sha) in FILES.items():
        r=urllib.request.Request(BASE+str(fid),headers={"User-Agent":"Original-reproductive-compensation-audit-v09/0.1","Accept":"*/*"})
        with urllib.request.urlopen(r,timeout=65) as response:blob=response.read(100000)
        if hashlib.sha256(blob).hexdigest()!=sha:
            raise ValueError("source original file checksum mismatch "+name)
        try:
            decoded=blob.decode("utf-8-sig")
            encoding="utf-8-sig"
        except UnicodeDecodeError:
            decoded=blob.decode("latin-1")
            encoding="latin-1"
        d=pd.read_csv(io.StringIO(decoded),na_values=["NA","NaN"])
        print("SOURCE_FILE_ENCODING",name,encoding,flush=True)
        d.columns=[str(c).strip() for c in d.columns]
        tables[name]=d
    return tables
def v(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return float(round(x,6)) if np.isfinite(x) else None
    if isinstance(x,(dict,)):return {str(k):v(val) for k,val in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [v(val) for val in x]
    return x

def bootstrap_by_genotype(frame,response,n=N_BOOT):
    if frame.empty:return {"n":0}
    groups=[g[response].dropna().to_numpy(float) for _,g in frame.groupby("genotype") if g[response].notna().any()]
    groups=[x for x in groups if len(x)]
    if not groups:return {"n":0}
    vals=np.concatenate(groups)
    bs=np.empty(n)
    for i in range(n):
        idx=RNG.integers(0,len(groups),size=len(groups))
        pool=np.concatenate([groups[j] for j in idx])
        bs[i]=pool.mean()
    return {"n_plants":len(vals),"n_genotypes":len(groups),
            "mean":v(vals.mean()),"median":v(np.median(vals)),
            "genotype_bootstrap95":[v(x) for x in np.quantile(bs,[.025,.975])]}

def prepare():
    d=original()
    flo=d["inflorescences.csv"].copy()
    total=d["total.csv"].copy()
    info=d["plant information.csv"].copy()
    for q in [flo,total,info]:
        for c in q.columns:
            if c.lower() in ("genotype","indi","plantid","id1","id"):
                q[c]=q[c].astype("string").str.strip()
    required_f=["ID1","ID2","treatment","xu","male","female","fruit","fruitset"]
    required_t=["genotype","indi","treatment","TFS","inflorno","femaleno","maleno","fitTFS","cormnumber"]
    required_i=["genotype","plantid","liveor","treatment"]
    for req,tab in ((required_f,flo),(required_t,total),(required_i,info)):
        if not set(req)<=set(tab):raise ValueError("unrecognized source file header "+str(set(req)-set(tab)))
    for c in ["treatment","xu","male","female","fruit","fruitset"]:
        flo[c]=pd.to_numeric(flo[c],errors="coerce")
    for c in ["treatment","TFS","inflorno","femaleno","maleno","fitTFS","cormnumber"]:
        total[c]=pd.to_numeric(total[c],errors="coerce")
    info["liveor"]=pd.to_numeric(info.liveor,errors="coerce")
    print("SOURCE_PLANT_ID_GRAIN",{
        "total_n":len(total),"n_total_id_unique":int(total.ID.nunique()),
        "n_total_indi_unique":int(total.indi.nunique()),
        "info_n":len(info),"n_info_id_unique":int(info.id.nunique()),
        "n_info_plantid_unique":int(info.plantid.nunique())},flush=True)
    if total.ID.duplicated().any():
        raise ValueError("Original total seeds table canonical ID unexpectedly duplicated")
    n_ambiguous_by_name=info.groupby("id")["treatment"].nunique()
    ambiguous=set(n_ambiguous_by_name[n_ambiguous_by_name>1].index)
    print("SOURCE_INFO_DUPLICATION_AUDIT",{
        "n_repeated_entries":int(info.id.duplicated().sum()),
        "n_ambiguous_ids_with_conflicting_named_treatment":len(ambiguous),
        "repeated_id_treatment_name_nunique_distribution":{
           str(k):int(v) for k,v in n_ambiguous_by_name.value_counts().items()}},flush=True)
    # Deduplicate only for the provenance mapping. The original file itself
    # remains intact; no survival values are modified or counted as zeros.
    info_for_mapping=info[~info.id.isin(ambiguous)].drop_duplicates("id")
    assert set(total.treatment.dropna().unique())<=set([1,2,5]),set(total.treatment.dropna().unique())
    assert set(flo.treatment.dropna().unique())<=set([1,2,5])
    # First confirm mapping numeric treatment -> site-level Net/Shade/Open
    mapping=info_for_mapping.merge(total[["ID","treatment"]],left_on="id",right_on="ID",how="inner",validate="one_to_one",suffixes=("_name","_code"))
    ct=pd.crosstab(mapping.treatment_code,mapping.treatment_name)
    print("SOURCE_TREATMENT_CROSSWALK",ct.to_dict(),flush=True)
    assignment={int(code):str(ct.columns[np.argmax(row.to_numpy())]) for code,row in ct.iterrows()}
    if any((row>0).sum()!=1 for _,row in ct.iterrows()) or len(set(assignment.values()))!=len(assignment):
        raise ValueError("source plant treatment mapping not one to one")
    return flo,total,info,assignment,ct

def main():
    flo,total,info,treatmap,crosstab=prepare()
    print("TREATMENT_MAPPING",treatmap,"source inflo rows",len(flo),"plant rows",len(total),"planted entries",len(info),flush=True)
    report={
      "status":"failed",
      "source":"Tang & Dai 2026 Ecol Evol doi:10.1002/ece3.73955",
      "original_data_doi":"10.6084/m9.figshare.27610539",
      "source_file_sha256_verified":{k:sha for k,(fid,sha) in FILES.items()},
      "primary_sample_grain":"clonal plant; inflorescences within plants; genotypes may repeat across treatments",
      "treatment_site_confounded_no_factorial_pollen_or_light_test":True,
      "source_tables":{"inflorescence_rows":len(flo),"total_rows":len(total),"planted_rows":len(info)},
      "treatment_code_to_name":treatmap,
      "site_assignment_crosstab":{str(index):{str(col):int(crosstab.loc[index,col]) for col in crosstab.columns} for index in crosstab.index},
      "by_treatment":{}
    }
    if len(flo)<200 or len(total)<70:raise ValueError("source sample unexpectedly low")
    # Meaning of 'xu' is sequential inflorescence order; confirm unique ID1,xu.
    report["duplicate_plant_inflo_order"]=int(flo.duplicated(["treatment","ID1","xu"]).sum())
    report["missing_by_source"]={x:{str(k):int(z) for k,z in d.isna().sum().items() if z>0} for x,d in [("inflorescences",flo),("total",total),("plants",info)]}
    # Preserve missing fruit/female as missing rather than falsely coding failure.
    reproductive=flo.dropna(subset=["ID1","xu","female","fruit","treatment"]).copy()
    reproductive=reproductive[reproductive.female>0]
    if (reproductive.fruit<0).any() or (reproductive.fruit>reproductive.female).any():
        raise ValueError("fruit count cannot exceed female flower count")
    # Total source records should track genetically distinct individual clones, not inflorescences.
    reproductive["plant"]=reproductive["ID1"].astype("string")
    grouped=[]
    for (plant,treatment),g in reproductive.groupby(["plant","treatment"]):
        tval=int(treatment); gname=treatmap.get(tval,"unknown")
        if g.treatment.nunique()!=1:raise ValueError("plant with multiple source treatments")
        genotype=str(g.ID2.iloc[0])
        first=g[g.xu==1];later=g[g.xu>1]
        if first.empty:continue
        f1=first.female.sum();r1=first.fruit.sum()
        flater=later.female.sum();rlater=later.fruit.sum()
        row={"plant":str(plant),"genotype":genotype,"treatment_code":tval,"treatment":gname,
           "n_inflorescences_with_female":int(g.xu.nunique()),
           "first_female":f1,"first_fruit":r1,"first_fruit_fraction":r1/f1,
           "late_female":flater,"late_fruit":rlater,
           "later_fruit_fraction":rlater/flater if flater>0 else np.nan,
           "late_minus_first_fraction":(rlater/flater-r1/f1) if flater>0 else np.nan,
           "fruit_fraction_all":(r1+rlater)/(f1+flater),
           "late_share_of_total_fruit":rlater/(r1+rlater) if (r1+rlater)>0 else np.nan,
           "first_zero_later_positive":bool(r1==0 and rlater>0)}
        grouped.append(row)
    frame=pd.DataFrame(grouped)
    if frame.empty:raise ValueError("no reproductive plant matching")
    valid=total.copy()
    valid["treatment_name"]=valid.treatment.map(treatmap)
    missing_total=valid.TFS.isna().sum()
    report["total_seed_count_missing_not_zero_filled"]=int(missing_total)
    report["n_genotypes_total"]=int(valid.genotype.nunique())
    report["genotypes_shared_between_treatments"]={
       str(k):int(v) for k,v in valid.groupby("genotype").treatment.nunique().value_counts().sort_index().items()}
    plantdiff=frame.merge(valid[["indi","treatment","TFS","inflorno","femaleno","genotype","treatment_name"]],
                          left_on=["plant","treatment_code"],right_on=["indi","treatment"],
                          how="left",validate="one_to_one",suffixes=("","_total"))
    report["repro_plants_to_total_matched"]=int(plantdiff.TFS.notna().sum())
    report["repro_plants_without_matched_total_seeds"]=int(plantdiff.TFS.isna().sum())
    if int(plantdiff.TFS.notna().sum())<60:raise ValueError("too many original plant-to-total mismatches")
    matched=plantdiff[plantdiff.TFS.notna()].copy()
    genotype_mismatch=(matched.genotype.astype(str)!=matched.genotype_total.astype(str))
    report["genotype_inconsistencies_on_matched_individuals"]=int(genotype_mismatch.sum())
    if genotype_mismatch.any():
        print("INFO: inflorescence ID2 label differs from source total.csv genotype in",
               int(genotype_mismatch.sum()),"matched records; using verified total.csv genotype",flush=True)
    matched["genotype"]=matched["genotype_total"]
    report["analysis_uses_only_matched_plant_level_seed_and_inflorescence_records"]=True
    # Genotype-stratified summaries by environment, separately single vs multi inflorescence plants.
    for code,name in sorted(treatmap.items()):
        w=matched[matched.treatment_code==code].copy()
        t=valid[valid.treatment==code].copy()
        p=w[w.late_female>0].copy()
        a={
          "n_original_total_plants":len(t),"n_genotypes_original_total":int(t.genotype.nunique()),
          "n_eligible_reproductive_plants_with_first_female_flower":len(w),
          "n_plant_pairs_with_first_and_later_female_flowers":len(p),
          "n_original_rows_in_this_treatment":int((flo.treatment==code).sum()),
          "plant_level_total_seeds_TFS":bootstrap_by_genotype(t.dropna(subset=["TFS"]),"TFS"),
          "plant_level_inflorescence_count":bootstrap_by_genotype(t.dropna(subset=["inflorno"]),"inflorno"),
          "first_inflorescence_fruit_rate":bootstrap_by_genotype(w,"first_fruit_fraction"),
          "later_inflorescence_fruit_rate_among_multi":bootstrap_by_genotype(p,"later_fruit_fraction"),
          "late_minus_first_fruit_fraction_same_plant":bootstrap_by_genotype(p,"late_minus_first_fraction"),
          "proportion_plants_initially_zero_but_later_fruit":{
             "n_initial_zero":int((p.first_fruit==0).sum()),
             "n_rescued_later":int(p.first_zero_later_positive.sum()),
             "prob_recovery_given_initial_zero":float(p.first_zero_later_positive.sum()/(p.first_fruit==0).sum()) if (p.first_fruit==0).sum()>0 else None},
          "later_share_of_fruits_among_multi":bootstrap_by_genotype(p.dropna(subset=["late_share_of_total_fruit"]),"late_share_of_total_fruit"),
          "genotype_counts":int(w.genotype.nunique())}
        report["by_treatment"][name]=a
    # Compare different sites only within the same SOURCE CLONAL GENOTYPE.
    # Site is nevertheless fully confounded with treatment: differences are
    # descriptive, not a causal isolated pollinator/light effect.
    report["cross_site_common_genotype_contrasts"]={}
    def genotype_pair_difference(frame,value,a,b):
        wide=frame.pivot_table(index="genotype",columns="treatment",values=value,aggfunc="first")
        if a not in wide or b not in wide:return {"n_shared_genotypes":0}
        pairs=wide[[a,b]].dropna()
        if len(pairs)<5:return {"n_shared_genotypes":len(pairs)}
        delta=(pairs[a]-pairs[b]).to_numpy(float)
        boot=delta[RNG.integers(0,len(delta),size=(N_BOOT,len(delta)))].mean(axis=1)
        signs=RNG.choice([-1.,1.],size=(N_PERM,len(delta)))
        null=(signs*delta).mean(axis=1)
        pval=(int((np.abs(null)>=abs(delta.mean())-1e-12).sum())+1)/(N_PERM+1)
        return {"n_shared_genotypes":len(pairs),"mean_a_minus_b":v(delta.mean()),
           "genotype_pair_bootstrap95":[v(z) for z in np.quantile(boot,[.025,.975])],
           "within_genotype_signflip_p":v(pval),
           "site_confounded_observational_contrast":True}
    totalpairs=valid[["genotype","treatment_name","TFS","inflorno"]].rename(columns={"treatment_name":"treatment"})
    for a,b in (("Net","Shade"),("Net","Open"),("Shade","Open")):
        contrast={}
        for field in ("TFS","inflorno"):
            contrast[field]=genotype_pair_difference(totalpairs,field,a,b)
        for field in ("first_fruit_fraction","later_fruit_fraction","late_minus_first_fraction"):
            contrast[field]=genotype_pair_difference(matched,field,a,b)
        report["cross_site_common_genotype_contrasts"][a+"_minus_"+b]=contrast
    report["competing_denominators"]={
       "plant_total_TFS":"plant-level total seeds in original authors source; not sum of fruit rates",
       "fruit_in_inflorescence":"fruits from female flowers in a single inflorescence, not paternal fitness",
       "late_result":"within-plant after first inflorescence, only among plants with first and later female flowers",
       "growth_or_genotype_selection":"survivorship and missing flowering outcomes can bias conditional contrast"}
    report["status"]="success_original_figshare_sequential_fruiting_analysis"
    (OUT/"sagittaria_trifolia_2026_original_sequential_compensation_v09.json").write_text(
       json.dumps(v(report),indent=2,ensure_ascii=False)+"\n")
    print("SAGITTARIA_ORIGINAL_REANALYSIS",json.dumps(v(report),ensure_ascii=False)[:30000],flush=True)
if __name__=="__main__":main()
