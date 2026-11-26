#!/usr/bin/env python3
"""Date-matched flower display and visitor-count source audit: Zurich 2016.

Identifies 2016 garden × fieldwork date × phytometer species from the study's
documented *per-date* flower counts and observation effort in minutes.
Joins the EnviDat time-window visitation matrix to the same garden/day/plant.

The original paper already tested visitor functional groups vs reproduction
(Reji Chacko et al. 2026, DOI 10.1111/1365-2664.70384). We only test an
unanswered display-normalized rate sensitivity, not claim new pollination
causality or re-label it original data discovery.

No private garden location or Id is exported; all primary rows stay in memory.
"""
from __future__ import annotations
import csv
import datetime as dt
import hashlib
import io
import json
import math
import time
import urllib.request
import zipfile
from collections import Counter,defaultdict
from pathlib import Path

OUT=Path("analysis_results_zurich_display_date")
OUT.mkdir(exist_ok=True)

STUDY="https://www.envidat.ch/dataset/d8e61f31-c2e5-4052-8784-b1624dfbfc40/resource/cef4a88c-0ef6-4327-a4a0-fc7d33b468e5/download/rejichacko_etal_2025_envidat.zip"
SOURCE="https://raw.githubusercontent.com/mrejichacko/BetterBlooms/main/cleaned_data/sampling_effort_gardens_plants_aggregated_per_date.txt"
FILE="06_trait_data/species_temporal_flower_visitation_matrix.csv"
META=("Id","phytometer_plant","capture_date","capture_window")
PLANTS={
 "carrot":("n_flowers_carrot","sampling_effort_min_carrot"),
 "radish":("n_flowers_radish","sampling_effort_min_radish"),
 "sainfoin":("n_flowers_sainfoin","sampling_effort_min_sainfoin"),
 "comfrey":("n_flowers_comfrey","sampling_effort_min_comfrey")
}
def load(url,limit):
    for attempt in range(3):
        try:
            q=urllib.request.Request(url,headers={"User-Agent":"Daphne-flower-resource-date-match/0.1 (research)"})
            with urllib.request.urlopen(q,timeout=120) as r:data=r.read(limit+1)
            if len(data)>limit:raise ValueError("Oversized response")
            return data
        except OSError:
            if attempt==2:raise
            time.sleep(2**attempt)

def numeric(x):
    try:
        v=float(str(x).strip())
        if math.isfinite(v):return v
    except (ValueError,TypeError):pass
    return None

def cleaned_rows():
    a=load(SOURCE,1_000_000)
    x=list(csv.DictReader(io.StringIO(a.decode("utf-8-sig","replace")),delimiter=";"))
    b=load(STUDY,80_000_000)
    if not zipfile.is_zipfile(io.BytesIO(b)):raise ValueError("EnviDat archive not ZIP")
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        raw=z.read(FILE)
    f=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace"))))
    return x,f,{"fieldwork_aggregated_source":SOURCE,"fieldwork_sha256":hashlib.sha256(a).hexdigest(),
                "visitation_source":STUDY,"visitation_table_sha256":hashlib.sha256(raw).hexdigest(),
                "visitation_zip_sha256":hashlib.sha256(b).hexdigest()}

def prepare():
    field,visits,source=cleaned_rows()
    tags=[k for k in visits[0] if k not in META]
    aggregate=defaultdict(lambda:[0.0,0])
    unknown=Counter()
    totals=0.
    field_dates=set()
    for rec in visits:
        try:day=dt.datetime.strptime(rec["capture_date"],"%d.%m.%Y").date()
        except ValueError:raise ValueError("Unexpected window date "+rec["capture_date"])
        plant=str(rec["phytometer_plant"]).strip().lower()
        if plant not in PLANTS:
            unknown[plant]+=1
            continue
        site=str(rec["Id"]).strip()
        vals=[numeric(rec.get(c)) for c in tags]
        if all(v is None for v in vals):
            unknown["all_missing_taxa_in_window"]+=1
            continue
        if any(v is not None and v<0 for v in vals):
            raise ValueError("Negative counted visitor")
        n=sum(v for v in vals if v is not None)
        aggregate[(site,day,plant)][0]+=n
        aggregate[(site,day,plant)][1]+=1
        totals+=n
        field_dates.add(day)
    if unknown.get("carrot") or unknown.get("radish"):
        raise ValueError("Unexpected plant labels")
    by_source={}
    for r in field:
        key=(str(r["Id"]).strip(),dt.date.fromisoformat(r["capture_period"]))
        if key in by_source: raise ValueError("Duplicate garden/date source")
        by_source[key]=r
    rows=[];mismatches=Counter();source_keys=set()
    for (site,day),r in by_source.items():
        for plant,(fkey,ekey) in PLANTS.items():
            f=numeric(r.get(fkey))
            effort=numeric(r.get(ekey))
            nwin,n=aggregate.get((site,day,plant),[0.0,0])
            if f is None or effort is None:raise ValueError("Missing numeric flower effort")
            if f<0 or effort<0:raise ValueError("Negative flowers or effort")
            if effort>0 and f>0:
                if n==0:mismatches["flowering_effort_without_matching_visit_window"]+=1
                rows.append({"site":site,"day":day,"plant":plant,
                    "flower_units":f,"effort_minutes":effort,
                    "all_visitor_count":nwin,"n_observation_windows":n,
                    "day_of_year":day.timetuple().tm_yday})
                source_keys.add((site,day,plant))
            else:
                if nwin>0:mismatches["visits_without_positive_flower_effort"]+=nwin
                if n>0:mismatches["windows_without_positive_flower_effort"]+=n
    orphan=[(key,num[0]) for key,num in aggregate.items() if key not in source_keys]
    return rows,{
        "fieldwork_site_date_rows":len(field),
        "matrix_rows":len(visits),
        "matrix_taxon_columns":len(tags),
        "n_visitors_total_from_matrix":totals,
        "matrix_unique_site_day_plant":len(aggregate),
        "fieldwork_unique_site_day":len(by_source),
        "site_day_plant_positive_F_and_effort":len(rows),
        "unique_sites":len({r["site"] for r in rows}),
        "unique_field_dates":len(field_dates),
        "field_date_range":[min(field_dates).isoformat(),max(field_dates).isoformat()],
        "joined_visitors":sum(r["all_visitor_count"] for r in rows),
        "joined_observation_windows":sum(r["n_observation_windows"] for r in rows),
        "number_of_missing_matching_window_rows":sum(r["n_observation_windows"]==0 for r in rows),
        "orphan_matrix_keys_with_no_positive_flower_effort":len(orphan),
        "orphan_counts_by_problem":dict(mismatches),
        "unknown_taxa_or_plant_entries":dict(unknown),
        "by_plant":{
            p:{"n_site_day_rows":sum(r["plant"]==p for r in rows),
               "n_positive_visitor_rows":sum(r["plant"]==p and r["all_visitor_count"]>0 for r in rows),
               "sum_observed_visitors":sum(r["all_visitor_count"] for r in rows if r["plant"]==p),
               "flower_units_range":[min(r["flower_units"] for r in rows if r["plant"]==p),
                                     max(r["flower_units"] for r in rows if r["plant"]==p)]}
            for p in PLANTS}
    },source

def main():
    result={"status":"failed","data_grain":"source garden×actual 2016 observation day×plant",
      "scope":"Genuine date-matched flowers and observation effort; NOT measured pollen effectiveness or final fitness",
      "source_study":"Reji Chacko et al. (2026), DOI 10.1111/1365-2664.70384; original EnviDat DOI 10.16904/envidat.676",
      "privacy":"Garden IDs/locations not exported",
      "caveats":["Original phytometer groups vary in umbels, flowers and species composition",
        "Garden site/date observations do not track genetic individuals or true seed-level pollen deposition",
        "Some 2016 visitors include recorded illegitimate pollinators and visitors not effective at pollination",
        "Sampling effort minutes depend on volunteer coverage; require offset",
        "Visit matrix capture-window rows need not each represent one equal amount of minutes",
        "The authors previously tested pollinator guilds vs fruit and seed set; avoid novelty claim"]}
    try:
        rows,qa,provenance=prepare()
        result["source"]=provenance
        result["qa"]=qa
        if len(rows)<180:raise ValueError("Not enough date-matched records for within-site analysis")
        if qa["number_of_missing_matching_window_rows"]>15:
            raise ValueError("Date/plant alignment failure")
        result["status"]="success"
    except Exception as exc:
        result["error"]=type(exc).__name__+": "+str(exc)
    (OUT/"site_date_flower_effort_visitors_qa.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,indent=2,ensure_ascii=False)[:25000])
    if result["status"]!="success":raise RuntimeError(result.get("error"))
if __name__=="__main__":main()
