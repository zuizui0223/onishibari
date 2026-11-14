#!/usr/bin/env python3
"""Source-backed, privacy-preserving GBIF spatial-overlap audit for Daphne.

All individual occurrence coordinates are processed in memory and discarded.
Exports only aggregate cell counts and pairwise overlap, not point locations.
This assesses sampling coverage, NOT a rigorous species distribution model.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ENDPOINT="https://api.gbif.org/v1/occurrence/search"
SPECIES={
    "Daphne pseudomezereum":5524168,
    "Daphne jezoensis":5523246,
    "Daphne kiusiana":5523189
}
HEADERS={"User-Agent":"onishibari-GBIF-spatial-audit/0.1 (research)"}
PAGE=300
GRID_DEGREES=0.25

def get(query):
    url=ENDPOINT+"?"+urlencode(query)
    last=None
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers=HEADERS),timeout=45) as res:
                return json.loads(res.read().decode("utf-8"))
        except (HTTPError,URLError,TimeoutError) as exc:
            last=exc
            time.sleep(2**attempt)
    raise RuntimeError(f"GBIF endpoint failed: {last}")

def download_taxon(taxon_key):
    params={"taxonKey":taxon_key,"country":"JP","hasCoordinate":"true",
            "occurrenceStatus":"PRESENT","limit":PAGE}
    first=get(dict(params,offset=0))
    count=first["count"]
    results=list(first.get("results",[]))
    while len(results)<count:
        offset=len(results)
        page=get(dict(params,offset=offset))
        rows=page.get("results",[])
        if not rows:
            raise RuntimeError("Unexpected end of pagination")
        results.extend(rows)
    if len(results)!=count:
        raise ValueError(f"Pagination inconsistent: {len(results)} != {count}")
    return results, count

def valid(row,key):
    try:
        lat=float(row["decimalLatitude"])
        lon=float(row["decimalLongitude"])
    except (TypeError,ValueError,KeyError):
        return "missing_coordinate"
    if not (20<=lat<=47 and 120<=lon<=155):
        return "outside_broad_japan_bbox"
    if row.get("hasGeospatialIssues"):
        return "gbif_geospatial_issues"
    if row.get("speciesKey") is not None and int(row["speciesKey"])!=key:
        return "species_key_mismatch"
    uncertainty=row.get("coordinateUncertaintyInMeters")
    if uncertainty is not None:
        try:
            if float(uncertainty)>10000:
                return "uncertainty_above_10km"
        except (ValueError,TypeError):
            return "unparseable_uncertainty"
    return None

def cell(row):
    return (int(math.floor(float(row["decimalLatitude"])/GRID_DEGREES)),
            int(math.floor(float(row["decimalLongitude"])/GRID_DEGREES)))

def lat_summary(cells):
    if not cells:
        return None
    lat=sorted((x[0]+0.5)*GRID_DEGREES for x in cells)
    def at(q):
        i=(len(lat)-1)*q
        lo=int(math.floor(i))
        hi=int(math.ceil(i))
        return round(lat[lo]*(hi-i)+lat[hi]*(i-lo),2) if lo!=hi else round(lat[lo],2)
    return {"p10":at(.1),"median":at(.5),"p90":at(.9)}

def one_species(rows, taxon_key, count):
    reason=Counter()
    accepted=[]
    seen_keys=set()
    for row in rows:
        k=row.get("key")
        if k in seen_keys:
            reason["duplicate_gbif_key"]+=1
            continue
        seen_keys.add(k)
        issue=valid(row,taxon_key)
        if issue:
            reason[issue]+=1
            continue
        accepted.append(row)
    cells={cell(row) for row in accepted}
    modern={cell(row) for row in accepted if str(row.get("year",""))>="2000"}
    known_uncertainty=sum(r.get("coordinateUncertaintyInMeters") is not None for r in accepted)
    bs=Counter(r.get("basisOfRecord") or "UNKNOWN" for r in accepted)
    return {
        "requested_georeferenced_count": count,
        "retained_records":len(accepted),
        "excluded_by_reason":dict(reason),
        "unique_025deg_cells":len(cells),
        "modern_post2000_unique_cells":len(modern),
        "latitude_of_occupied_cell_centers":lat_summary(cells),
        "basis_of_record":dict(bs),
        "uncertainty_reported_in_retained":known_uncertainty,
    }, cells, modern

def pairs(keys,sets):
    output={}
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            both=len(sets[a]&sets[b])
            union=len(sets[a]|sets[b])
            output[a+" / "+b]={
                "shared_025deg_cells":both,
                "union_025deg_cells":union,
                "jaccard":round(both/union,4) if union else None,
                "overlap_fraction_a":round(both/len(sets[a]),4) if sets[a] else None,
                "overlap_fraction_b":round(both/len(sets[b]),4) if sets[b] else None,
            }
    return output

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="data/derived/gbif_daphne_spatial_aggregates.json")
    args=parser.parse_args()
    summary,allcells,moderncells={}, {}, {}
    for name,key in SPECIES.items():
        rows,count=download_taxon(key)
        result,cells,newcells=one_species(rows,key,count)
        summary[name]=result
        allcells[name]=cells
        moderncells[name]=newcells
        print("SPATIAL_SPECIES="+json.dumps({"name":name,**result},sort_keys=True),flush=True)
    report={
        "generated_utc":datetime.now(timezone.utc).isoformat(),
        "source":ENDPOINT,
        "query_filters":{"country":"JP","hasCoordinate":"true",
                         "occurrenceStatus":"PRESENT","taxonKey":SPECIES},
        "analysis_filters":{
            "bbox":"20–47 N, 120–155 E (broad conservative Japan envelope)",
            "hasGeospatialIssues":False,
            "max_known_coord_uncertainty_m":10000,
            "coordinate_uncertainty_missing":"retained but flagged",
            "grid_degrees":GRID_DEGREES,
            "modern_year_min":2000,
            "note":"No coordinates or spatial cell identifiers exported"
        },
        "species":summary,
        "pairwise_all_dates":pairs(list(SPECIES),allcells),
        "pairwise_year_2000_plus":pairs(list(SPECIES),moderncells),
        "limitations":[
            "A 0.25-degree cell represents presence of at least one reported location, not abundance.",
            "GBIF records combine different survey effort, years, coordinates and taxonomy.",
            "Cultivated, misidentified and repeated occurrences can remain after quality filters.",
            "0.25-degree spatial overlap is scale-dependent and cannot establish true biological sympatry.",
            "Latitude is not temperature; do not infer causality in bloom timing or reproduction.",
            "No exact rare plant coordinates, site names, GBIF record IDs or location grids exported."
        ]
    }
    dest=Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    for k,v in report["pairwise_all_dates"].items():
        print("SPATIAL_PAIR="+json.dumps({"pair":k,**v},sort_keys=True),flush=True)
    print("SPATIAL_OUTPUT="+str(dest),flush=True)

if __name__=="__main__":
    main()
