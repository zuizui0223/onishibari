#!/usr/bin/env python3
"""v0.21: matched-sex floral opportunity reanalysis of original Daphne jezoensis 2018 censuses.

Quantities: observed-window integral of fraction of initial-inflorescence flowers open.
NOT exact individual floral longevity, NOT dry-flower maintenance cost, NOT seed fitness.
No raw plant IDs, coordinates or source rows are emitted.
"""
from __future__ import annotations
import datetime,hashlib,io,itertools,json,urllib.request,zipfile
from pathlib import Path
import numpy as np,pandas as pd
OUT=Path("analysis_results_daphne_flower_days_v21")
ZIP_SHA="97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8"
CSV_SHA="140aeae0064ff9e13ff4199794abe7dc533d5ac4f3fa0e929fc815eb0d4df017"
RNG=np.random.default_rng(20261116)
BOOTS=18000

def pull(url,limit):
    req=urllib.request.Request(url,headers={"User-Agent":"Daphne-published-open-flower-days-study-v21","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=80) as res:b=res.read(limit+1)
    if len(b)>limit:raise RuntimeError("Oversize public source")
    return b

def read():
    j=json.loads(pull("https://zenodo.org/api/records/8253921",2_000_000))
    fs=j["files"];entries=fs.get("entries",{}) if isinstance(fs,dict) else {x["key"]:x for x in fs}
    a=entries["Data for Figures_Analyses .zip"]["links"]
    z=pull(a.get("content") or a["self"],1_000_000)
    assert hashlib.sha256(z).hexdigest()==ZIP_SHA
    with zipfile.ZipFile(io.BytesIO(z)) as f:
        names=[n for n in f.namelist() if Path(n).name=="FlowerLongevity.csv"]
        assert len(names)==1
        b=f.read(names[0])
    assert hashlib.sha256(b).hexdigest()==CSV_SHA
    d=pd.read_csv(io.BytesIO(b))
    assert len(d)==408 and set(d.Sex)=={"F","H"} and d.RametID.nunique()==14
    assert set(d.Treatment)=={"natural","bag","outcross","emasculate_bag","emasculate_outcross","self"}
    d["count"]=pd.to_numeric(d.FloweringNo,errors="coerce")
    d["total"]=pd.to_numeric(d.FlowerNo,errors="raise")
    assert d["count"].isna().sum()==9
    assert d["count"].dropna().ge(0).all()
    assert (d["count"].isna() | d["count"].le(d["total"])).all()
    assert d.total.gt(0).all()
    # Original date source is M/D with NO calendar year; don't interpret as actual year 0001.
    def day(s):
        v=datetime.datetime.strptime("2020/"+str(s),"%Y/%m/%d")
        return (v-datetime.datetime(2020,4,16)).days
    d["day"]=d.Date.map(day).astype(int)
    days=sorted(d.day.unique())
    assert len(days)==8 and days[0]==0 and days[-1]==24
    units=["RametID","Sex","Treatment","InfloID"]
    assert not d.duplicated(units+["day"]).any()
    n=d.groupby(units).size()
    assert len(n)==51 and n.eq(8).all()
    # The known full flower count must be constant per inflorescence/census.
    constant=d.groupby(units).total.nunique()
    assert constant.eq(1).all(),"FlowerNo must be a fixed denominator per original inflorescence"
    return d,days

def trap(y,x):
    return float(np.trapezoid(np.asarray(y,dtype=float),np.asarray(x,dtype=float)))

def summarize_units(d,days):
    rows=[]
    for (plant,sex,treat,inflo),g in d.groupby(["RametID","Sex","Treatment","InfloID"]):
        g=g.sort_values("day")
        assert list(g.day)==days
        den=float(g.total.iloc[0])
        fraction=g["count"].to_numpy(dtype=float)/den
        mask=np.isnan(fraction)
        lower=np.where(mask,0,fraction)
        upper=np.where(mask,1,fraction)
        if mask.any():
            q=pd.Series(fraction,index=days).interpolate(method="index",limit_direction="both").to_numpy()
        else:q=fraction
        # 'exposure' in observed window, not intrinsic longevity.
        rows.append({"private_id":plant,"sex":sex,"treatment":treat,
           "n_missing":int(mask.sum()),"flowers_total":int(den),
           "auc_low":trap(lower,days),"auc_high":trap(upper,days),
           "auc_interpolate":trap(q,days),
           "observed_days_fraction_positive":int(np.sum(np.nan_to_num(fraction)>0)),
           "last_positive_day":int(max(np.asarray(days)[(np.nan_to_num(fraction)>0)])) if np.any(np.nan_to_num(fraction)>0) else None,
           "first_positive_day":int(min(np.asarray(days)[(np.nan_to_num(fraction)>0)])) if np.any(np.nan_to_num(fraction)>0) else None})
    return pd.DataFrame(rows)

def describe(g):
    return {"n_ramets":int(g.private_id.nunique()),
      "n_inflo":len(g),
      "n_missing":int(g.n_missing.sum()),
      "n_inflo_with_missing":int((g.n_missing>0).sum()),
      "mean_auc_low":round(float(g.auc_low.mean()),6),
      "mean_auc_high":round(float(g.auc_high.mean()),6),
      "mean_auc_interpolated_sensitivity":round(float(g.auc_interpolate.mean()),6),
      "sd_auc_interpolated":round(float(g.auc_interpolate.std(ddof=1)),6) if len(g)>1 else None,
      "median_auc_interpolated":round(float(g.auc_interpolate.median()),6),
      "mean_last_positive_census_day":round(float(g.last_positive_day.mean()),5)}

def bootstrap_unpaired(a,b):
    x=a.to_numpy(dtype=float);y=b.to_numpy(dtype=float)
    ia=RNG.integers(0,len(x),(BOOTS,len(x)))
    ib=RNG.integers(0,len(y),(BOOTS,len(y)))
    boot=x[ia].mean(axis=1)-y[ib].mean(axis=1)
    return [round(float(v),6) for v in np.quantile(boot,[.025,.975])]

def contrast_within_plant(df,sex,treatment0,treatment1):
    a=df[(df.sex==sex)&(df.treatment==treatment0)]
    b=df[(df.sex==sex)&(df.treatment==treatment1)]
    assert not a.private_id.duplicated().any() and not b.private_id.duplicated().any()
    p=a.merge(b,on="private_id",suffixes=("_ref","_treat"),validate="1:1")
    assert len(p)==(5 if sex=="F" else 9)
    lower=(p.auc_low_treat-p.auc_high_ref).to_numpy()
    upper=(p.auc_high_treat-p.auc_low_ref).to_numpy()
    point=(p.auc_interpolate_treat-p.auc_interpolate_ref).to_numpy()
    return {"sex":sex,"contrast":f"{treatment1} minus {treatment0}",
      "matched_ramets":len(p),
      "mean_lower_bound":round(float(lower.mean()),6),
      "mean_upper_bound":round(float(upper.mean()),6),
      "interpolated_mean":round(float(point.mean()),6),
      "interpolated_median":round(float(np.median(point)),6),
      "n_treat_larger":int((point>0).sum()),
      "n_treat_smaller":int((point<0).sum()),
      "n_equal":int((point==0).sum()),
      "n_matched_with_any_missing":int(((p.n_missing_ref+p.n_missing_treat)>0).sum()),
      "interpretation":"Exploratory label comparison. No source-verified random allocation or hand-pollination timing; not a causal contrast."}

def main():
    OUT.mkdir(exist_ok=True)
    d,days=read()
    u=summarize_units(d,days)
    groups={(sex,treat):describe(g) for (sex,treat),g in u.groupby(["sex","treatment"])}
    f=u[(u.sex=="F")&(u.treatment=="bag")]
    h=u[(u.sex=="H")&(u.treatment=="emasculate_bag")]
    assert len(f)==5 and len(h)==9 and not (f.n_missing.sum()+h.n_missing.sum())
    res={
      "status":"source_original_2018_flower_days_v21",
      "paper_doi":"10.1002/ajb2.16239",
      "zenodo_doi":"10.5281/zenodo.8253921",
      "source_zip_sha256":ZIP_SHA,"source_table_sha256":CSV_SHA,
      "n_source_censuses":len(d),"n_inflorescences":len(u),"n_unique_ramets":int(u.private_id.nunique()),
      "n_missing_open_flower_counts":int(d["count"].isna().sum()),
      "source_date_format":"M/D; no calendar year; 0=April16 through 24=May10, all same 8 census days",
      "census_day_offsets":days,
      "units":"observed-window integral of recorded flowering flowers / original flowers, days",
      "denominator":"original total flowers per inflorescence, invariant over censuses",
      "groups":{sex+"|"+treat:v for (sex,treat),v in groups.items()},
      "primary_unpollinated_sex_comparison":{
        "F_n":len(f),"H_n":len(h),
        "F_minus_H_auc_days":round(float(f.auc_low.mean()-h.auc_low.mean()),6),
        "F_to_H_auc_ratio":round(float(f.auc_low.mean()/h.auc_low.mean()),6),
        "sex_group_boot_CI95":bootstrap_unpaired(f.auc_low,h.auc_low),
        "n_missing_primary":0,
        "not_causal_sex_effect":True,
        "H_emasculation_is_not_identical_manipulation_to_F_bagging":True},
      "same_plant_treatment_label_sensitivities":[
        contrast_within_plant(u,"F","bag","outcross"),
        contrast_within_plant(u,"F","bag","natural"),
        contrast_within_plant(u,"H","emasculate_bag","emasculate_outcross"),
        contrast_within_plant(u,"H","emasculate_bag","self"),
        contrast_within_plant(u,"H","emasculate_bag","natural")],
      "measurement_boundaries":{
        "not_exact_flower_individual_lifetime":True,
        "not_measured_resource_maintenance_mass_or_CO2":True,
        "not_fitness_or_compatible_pollination":True,
        "not_focal_Daphne_pseudomezereum":True,
        "literal_na_is_missing_not_zero":True}}
    p=OUT/"daphne_2018_flower_opportunity_auc_source_v21.json"
    p.write_text(json.dumps(res,ensure_ascii=False,indent=2)+"\n")
    print("DAPHNE_FLOWER_DAYS_ORIGINAL_V21",json.dumps(res,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
