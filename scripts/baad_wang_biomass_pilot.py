#!/usr/bin/env python3
"""BAAD directly measured woody dry mass × Wang sexual systems (exploratory).

BAAD v1.0.0 source individual plants, study/sampling selection not balanced.
Population/species allocation, floral investment, pollen delivery and NDVI are
NOT inferred. A pooled height-mass relationship is only used to show how
ontogeny shifts apparent sex-system mass differences; no cause inference.
"""
from collections import Counter,defaultdict
import csv,hashlib,io,json,math,statistics,time,urllib.request,zipfile
from pathlib import Path
import numpy as np

OUT=Path("analysis_results_baad_wang")
OUT.mkdir(exist_ok=True)
BAAD="https://github.com/dfalster/baad/releases/download/v1.0.0/baad_data.zip"
WANG="https://zenodo.org/api/records/4895147"
DATA="baad_data/baad_data.csv"
DICTIONARY="baad_data/baad_dictionary.csv"

def grab(url,maxsize=20_000_000):
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Daphne-BAAD-allometry-ecology/0.2"})
            with urllib.request.urlopen(req,timeout=120) as r:b=r.read(maxsize+1)
            if len(b)>maxsize:raise ValueError("Download too large")
            return b
        except OSError:
            if i==2:raise
            time.sleep(2**i)

def csvbytes(data):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig",errors="replace"))))

def clean(name):
    return " ".join(str(name or "").replace("_"," ").split()).lower()

def pos(x):
    try:
        v=float(x)
        if math.isfinite(v) and v>0:return v
    except (TypeError,ValueError):pass
    return None

def median(vals):
    return round(statistics.median(vals),5) if vals else None

def main():
    result={"source":"Falster et al 2015 BAAD version 1.0.0, Wang et al 2021 deposited sexual system labels",
      "status":"failed","mass_measure":"individual m.to dry biomass, NOT floral allocation or site biomass density",
      "biases":["individuals from different ontogenetic ages/studies, species are NOT sampled to typical size",
         "woody plants only, no herbaceous extrapolation",
         "species labels in Wang represent sex system and not functional performance",
         "family, study, environment and height confounding",
         "height-adjusted residuals are exploratory rather than causal sex-system effects",
         "BAAD individual biomass doesn't measure reproductive allocation or flower biomass"]}
    try:
        raw=grab(BAAD)
        if not zipfile.is_zipfile(io.BytesIO(raw)):raise ValueError("BAAD not ZIP")
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            data=csvbytes(z.read(DATA));dct=csvbytes(z.read(DICTIONARY))
        meta=json.loads(grab(WANG,3_000_000))
        x=next(f for f in meta["files"] if f.get("key")=="SplistWithSexSyst.csv")
        wblob=grab(x["links"]["self"],6_000_000)
        wang=csvbytes(wblob)
        classmap={}
        for w in wang:
            name=clean(w.get("Genus_species"))
            sex=clean(w.get("SexSyst"))
            if name and sex in ("dioecy","monoecy","hermaphroditism"):
                if name in classmap and classmap[name]!=sex:raise ValueError("Wang conflict")
                classmap[name]=sex
        result["BAAD_records"]=len(data)
        result["Wang_records"]=len(wang)
        result["mass_field_definitions"]=[r for r in dct if r.get("variable") in ("m.to","m.lf","m.st","m.rt","h.t")]
        result["source_sha256"]={"BAADzip":hashlib.sha256(raw).hexdigest(),
            "WangSpecies":hashlib.sha256(wblob).hexdigest()}
        seen=Counter();records=[]
        for r in data:
            name=clean(r.get("speciesMatched") or r.get("species"))
            if not name:continue
            seen[name]+=1
            mass=pos(r.get("m.to"))
            height=pos(r.get("h.t"))
            if mass is not None:
                records.append({"species":name,"mass":mass,"height":height,
                    "family":clean(r.get("family")),"study":r.get("studyName"),"sex":classmap.get(name)})
        result["BAAD_species_unique"]=len(seen)
        result["BAAD_species_matching_Wang"]=len(set(seen)&set(classmap))
        result["individuals_with_valid_total_mass"]=len(records)
        matched=[r for r in records if r["sex"]]
        result["individuals_with_mass_and_sex_class"]=len(matched)
        result["species_with_mass_and_sex_class"]=len({r["species"] for r in matched})
        result["individuals_with_mass_height_and_sex_class"]=sum(r["height"] is not None for r in matched)
        result["sex_system_counts_by_individual"]=dict(Counter(r["sex"] for r in matched))
        result["sex_system_counts_by_species"]=dict(Counter(
            next(r["sex"] for r in matched if r["species"]==n) for n in {r["species"] for r in matched}))
        withheight=[r for r in matched if r["height"] is not None]
        if len(withheight)<50:raise ValueError("Mass and height coverage insufficient")
        X=np.array([[1,math.log(r["height"])] for r in withheight])
        y=np.array([math.log(r["mass"]) for r in withheight])
        beta=np.linalg.lstsq(X,y,rcond=None)[0]
        result["pooled_log_mass_height_slope"]=round(float(beta[1]),5)
        result["pooled_intercept"]=round(float(beta[0]),5)
        S=defaultdict(list)
        for r in withheight:
            res=math.log(r["mass"])-float(beta[0]+beta[1]*math.log(r["height"]))
            S[r["species"]].append({"sex":r["sex"],"family":r["family"],
                "study":r["study"],"mass_g":r["mass"],"height_m":r["height"],
                "residual":res})
        ssummary=[]
        for name,arr in sorted(S.items()):
            ssummary.append({"species":name,"family":arr[0]["family"],"sex":arr[0]["sex"],
                "sampled_individuals":len(arr),
                "species_sample_median_mass_kg":median([x["mass_g"] for x in arr]),
                "species_sample_median_height_m":median([x["height_m"] for x in arr]),
                "pooled_height_adjusted_log_mass_residual":median([x["residual"] for x in arr]),
                "n_studies":len(set(x["study"] for x in arr))})
        result["species_with_valid_total_mass_and_height"]=len(ssummary)
        result["comparisons"]={}
        for sex in ["hermaphroditism","dioecy","monoecy"]:
            gr=[s for s in ssummary if s["sex"]==sex]
            result["comparisons"][sex]={
                "n_species":len(gr),
                "n_individuals":sum(s["sampled_individuals"] for s in gr),
                "median_species_sample_mass_kg":median([s["species_sample_median_mass_kg"] for s in gr]),
                "median_species_sample_height_m":median([s["species_sample_median_height_m"] for s in gr]),
                "median_height_adjusted_log_mass_residual":median([s["pooled_height_adjusted_log_mass_residual"] for s in gr]),
                "n_distinct_families":len(set(s["family"] for s in gr))}
        families=defaultdict(set)
        for s in ssummary:families[s["family"]].add(s["sex"])
        result["families_with_more_than_one_sex_system"]=sum(len(x)>1 for x in families.values())
        by_family={}
        for s in ssummary:
            family=s["family"]
            by_family.setdefault(family,{})
            by_family[family].setdefault(s["sex"],[])
            by_family[family][s["sex"]].append(s)
        result["multi_sex_system_families"]={
            family:{sex:{"n_species":len(spp),
              "median_height_adjusted_log_mass_residual":
                median([q["pooled_height_adjusted_log_mass_residual"] for q in spp])}
              for sex,spp in value.items()}
            for family,value in sorted(by_family.items()) if len(value)>1
        }
        if len(ssummary)>=20:
            fields=list(ssummary[0])
            with (OUT/"height_adjusted_species_baad_mass.csv").open("w",newline="") as o:
                wr=csv.DictWriter(o,fieldnames=fields);wr.writeheader();wr.writerows(ssummary)
        result["status"]="success"
    except Exception as e:
        result["error"]=type(e).__name__+": "+str(e)
    (OUT/"BAAD_vs_Wang_direct_biomass_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2)[:14000])
    if result["status"]!="success":raise RuntimeError(result["error"])
if __name__=="__main__":main()
