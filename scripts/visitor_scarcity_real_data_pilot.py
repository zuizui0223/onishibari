#!/usr/bin/env python3
"""2025 Zurich garden quantitative visitation and reproduction data pilot.

Calculate visitor detections per RECORD (window duration unknown from matrix
alone); zero-visit row fraction; site × plant visitation vs within-plant
reproductive outcome for matched anonymized garden IDs. Garden coordinates
never copied/exported; 24 sites is a local ecological case study, not a
cross-angiosperm evolutionary model.
"""
import csv,hashlib,io,json,math,time,urllib.request,zipfile,datetime
from collections import Counter,defaultdict
from pathlib import Path
from scipy.stats import spearmanr

OUT=Path("analysis_results_visitor_link")
OUT.mkdir(exist_ok=True)
URL=("https://www.envidat.ch/dataset/d8e61f31-c2e5-4052-8784-b1624dfbfc40/"
 "resource/cef4a88c-0ef6-4327-a4a0-fc7d33b468e5/download/"
 "rejichacko_etal_2025_envidat.zip")
FILE="06_trait_data/species_temporal_flower_visitation_matrix.csv"
COMMON_PLANT={"daucus carota":"carrot", "onobrychis viciifolia":"sainfoin",
 "raphanus sativus":"radish", "symphytum officinale":"comfrey"}
SUCCESS={
 "daucus carota":("07_pollination_success/daucus_carota_seed_set.csv","seed_per_umbel"),
 "onobrychis viciifolia":("07_pollination_success/onobrychis_viciifolia_fruit_set.csv","fruit_per_flower"),
 "raphanus sativus":("07_pollination_success/raphanus_sativus_fruit_set.csv","fruit_per_flower"),
 "symphytum officinale":("07_pollination_success/symphytum_officinale_fruit_set.csv","seeded_flower_fraction")}

def getzip():
    for attempt in range(3):
        try:
            req=urllib.request.Request(URL,headers={"User-Agent":"onishibari-2025-visitation-audit/0.2"})
            with urllib.request.urlopen(req,timeout=100) as r:b=r.read(80_000_001)
            if len(b)>80_000_000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise ValueError("bad or oversized ZIP")
            return b
        except OSError:
            if attempt==2:raise
            time.sleep(2**attempt)

def rows(z,name):
    text=z.read(name).decode("utf-8-sig",errors="replace")
    return list(csv.DictReader(io.StringIO(text)))

def num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def excluded(row):
    return any(str(row.get(k,"")).strip().lower() in ("1","yes","true","excluded","exclude")
       for k in ("exclude","excluded"))

def site_garden(k):
    return str(k).strip().lower()

def plant(name):
    return str(name).strip().replace("_"," ").lower()

def main():
    summary={"source":"Reji Chacko, Moretti, Frey (2025), EnviDat DOI:10.16904/envidat.676",
      "study_scope":"24 Zurich home gardens, four standardized phytometer plant species",
      "sensitive_coordinate_export":False,
      "units":"visits per CSV-recorded observation window; actual durations require raw field protocol validation",
      "species_reproductive_outcomes_are_not_comparable_units":True,
      "no_direct_pollen_deposition":True,"status":"starting"}
    try:
        data=getzip()
        summary["source_zip_sha256"]=hashlib.sha256(data).hexdigest()
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            mat=rows(z,FILE)
            fields=list(mat[0])
            flower_taxa=[f for f in fields if f not in ("Id","phytometer_plant","capture_date","capture_window")]
            summary["matrix_rows"]=len(mat);summary["visitor_columns"]=len(flower_taxa)
            grid={}
            total_visits=0
            day_counts=Counter()
            key_counts=Counter()
            all_groups=defaultdict(list)
            missing_only=0
            for row in mat:
                site=site_garden(row.get("Id"));p=plant(row.get("phytometer_plant"))
                date=row.get("capture_date","");window=row.get("capture_window","")
                key_counts[(site,p,date,window)]+=1
                positives=[]
                observed=0
                for key in flower_taxa:
                    val=num(row.get(key))
                    if val is not None:
                        observed+=1
                        if val>0:positives.append((key,val))
                if observed==0:
                    missing_only+=1
                    continue
                vcount=sum(value for _,value in positives)
                total_visits+=vcount
                agg=all_groups[(site,p)]
                agg.append((vcount,len(positives),date,window))
                day_counts[date]+=1
            summary["n_dates"]=len(day_counts)
            parsed_dates=[]
            for name in day_counts:
                try:parsed_dates.append(datetime.datetime.strptime(name,"%d.%m.%Y").date())
                except ValueError:continue
            summary["date_distribution_summary"]={
                "earliest":min(parsed_dates).isoformat() if parsed_dates else "",
                "latest":max(parsed_dates).isoformat() if parsed_dates else "",
                "n_date_strings_parsed":len(parsed_dates),
                "max_observation_rows_one_day":max(day_counts.values()) if day_counts else 0}
            summary["garden_plant_groups"]=len(all_groups)
            summary["matrix_nonmissing_observation_rows"]=sum(len(x) for x in all_groups.values())
            summary["matrix_rows_all_taxa_missing"]=missing_only
            summary["total_visitor_count_from_matrix"]=total_visits
            summary["unique_identifiers_site_plant_date_window"]=len(key_counts)
            summary["duplicate_key_rows"]=sum(v-1 for v in key_counts.values() if v>1)
            summary["zero_records_n_and_fraction"]={
               "n":sum(x[0]==0 for v in all_groups.values() for x in v),
               "fraction":round(sum(x[0]==0 for v in all_groups.values() for x in v)/
                                max(1,sum(len(v) for v in all_groups.values())),5)}
            byplant=defaultdict(list)
            for (site,p),v in all_groups.items():
                byplant[p].append(dict(site=site,n_windows=len(v),
                     n_visits=sum(x[0] for x in v),
                     p_nonzero=sum(x[0]>0 for x in v)/len(v),
                     visits_per_record=sum(x[0] for x in v)/len(v)))
            summary["visitation_summary_by_plant"]={}
            for p, v in sorted(byplant.items()):
                n=sum(r["n_windows"] for r in v)
                N=sum(r["n_visits"] for r in v)
                summary["visitation_summary_by_plant"][p]={
                  "unique_gardens":len(v),"recorded_windows":n,
                  "total_visits":N,
                  "average_visits_per_record":round(N/max(n,1),5),
                  "sites_with_zero_total_visits":sum(x["n_visits"]==0 for x in v)}
            results={}
            for p,(filename,outcome) in SUCCESS.items():
                sr=rows(z,filename)
                numerator=defaultdict(float);denominator=defaultdict(float)
                discarded=0
                for row in sr:
                    if excluded(row):discarded+=1;continue
                    site=site_garden(row.get("Id"))
                    if outcome=="seed_per_umbel":
                        count=num(row.get("n_seeds"))
                        if count is None:continue
                        numerator[site]+=count;denominator[site]+=1
                    elif outcome=="fruit_per_flower":
                        yes=num(row.get("n_flowers_with_fruits"))
                        no=num(row.get("n_flowers_without_fruits"))
                        if yes is None or no is None or yes+no==0:continue
                        numerator[site]+=yes;denominator[site]+=yes+no
                    elif outcome=="seeded_flower_fraction":
                        yes=num(row.get("n_flowers_with_seeds"))
                        no=num(row.get("n_flowers_without_seeds"))
                        if yes is None or no is None or yes+no==0:continue
                        numerator[site]+=yes;denominator[site]+=yes+no
                joined=[]
                for x in byplant.get(COMMON_PLANT[p],[]):
                    site=x["site"]
                    if denominator[site]>0:
                        joined.append((x["visits_per_record"],numerator[site]/denominator[site]))
                corr=spearmanr([x[0] for x in joined],[x[1] for x in joined]) if len(joined)>=6 else None
                results[p]={
                  "reproductive_endpoint":outcome,
                  "original_reproductive_table_rows":len(sr),
                  "excluded_records_count":discarded,
                  "n_gardens_with_reproductive_data":sum(v>0 for v in denominator.values()),
                  "n_sites_with_both_visitor_windows_and_success":len(joined),
                  "visit_vs_reproductive_outcome_spearman_rho":
                     round(float(corr.statistic),5) if corr is not None and math.isfinite(float(corr.statistic)) else None,
                  "unadjusted_p":round(float(corr.pvalue),6) if corr is not None and math.isfinite(float(corr.pvalue)) else None,
                  "warning":"No effort duration, plant flower number, species-specific visitation efficiency, urban covariates, repeated dates or site hierarchy adjusted."
                }
            summary["site_level_visits_vs_reproduction"]=results
        summary["status"]="success"
    except Exception as ex:
        summary["status"]="failure";summary["error"]=type(ex).__name__+": "+str(ex)
    (OUT/"visitation_scarcity_pilot.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(summary,ensure_ascii=False,indent=2)[:25000])
    if summary["status"]!="success":raise RuntimeError(summary["error"])
if __name__=="__main__":main()
