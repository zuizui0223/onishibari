#!/usr/bin/env python3
"""Inspect paired published Mercurialis sex-fitness scenarios without exporting IDs.
Source: Chen & Pannell 2025 (New Phytol), pinned Zenodo v15644955.
This is a schema/semantics probe; no causal effects are estimated here.
"""
import csv, hashlib, io, json, urllib.request, zipfile, statistics, math, collections
U="https://zenodo.org/api/records/15644955/files/Geito_NewPhyto_v2.zip/content"
SHA="a18b7d296110be1e9817d6e1a3994365a4e7ae0aab50435e7e6883e1f3de4042"
b=urllib.request.urlopen(urllib.request.Request(U,headers={"User-Agent":"Source-meaning-v0.52"}),timeout=90).read()
assert hashlib.sha256(b).hexdigest()==SHA
z=zipfile.ZipFile(io.BytesIO(b))
def read(f):return list(csv.DictReader(io.StringIO(z.read("Geito_NewPhyto_v2/"+f).decode("utf-8-sig"))))
field=read("3G.csv");fitness=read("3G_ID01.csv")
def num(v):
 try:
  x=float(v)
  return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None
parents={(r["Population"],r["ID"]):r for r in field}
pairs=collections.defaultdict(dict)
for r in fitness:
 population = r["Population"].replace("3G_P","Population ")
 pairs[(population,r["ID"])][r["Inbreeding"]]=r
assert len(pairs)==181 and all(set(v)=={"ID=0","ID=1"} for v in pairs.values())
if set(pairs)!=set(parents):
 print("JOIN_KEY_DIAGNOSTIC",json.dumps({
  "field_population_categories":dict(collections.Counter(k[0] for k in parents)),
  "fitness_population_categories":dict(collections.Counter(k[0] for k in pairs)),
  "shared_exact_keys":len(set(pairs)&set(parents)),
  "field_ID_numeric_format":sum(k[1].isdigit() for k in parents),
  "fitness_ID_numeric_format":sum(k[1].isdigit() for k in pairs),
  "same_ID_ignore_population":len({k[1] for k in pairs}&{k[1] for k in parents}),
  "field_ID_unique":len({k[1] for k in parents}),
  "fitness_ID_unique":len({k[1] for k in pairs}),
 }))
 raise SystemExit(0) # diagnostic-only; never misjoin rows
cols=["Male_fitness","Female_fitness","Total_fitness"]
summary={}
for col in cols:
 a=[num(p["ID=0"][col]) for p in pairs.values()]
 b=[num(p["ID=1"][col]) for p in pairs.values()]
 q=[(x,y) for x,y in zip(a,b) if x is not None and y is not None]
 d=[x-y for x,y in q]
 summary[col]={"paired_nonmissing":len(q),"nonneg0":sum(x>=0 for x,y in q),"nonneg1":sum(y>=0 for x,y in q),"scenario0_mean":sum(x for x,y in q)/len(q),"scenario1_mean":sum(y for x,y in q)/len(q),"delta0_minus_delta1_min":min(d),"delta0_minus_delta1_max":max(d),"delta0_minus_delta1_mean":sum(d)/len(d),"delta0_minus_delta1_zero":sum(abs(x)<1e-10 for x in d)}
recorded={"n":len(fitness),"unique_parents":len(pairs),"modelled_field_parents":len(parents),"sex_fitness":summary}
for col in cols:
 print("FITNESS_PAIR_PROBE",json.dumps({col:summary[col]},ensure_ascii=False))
print("SCENARIO_SUMMARY",json.dumps({k:v for k,v in recorded.items() if k!="sex_fitness"}))
# Test whether published selfing reduction is exactly the selfing fraction in
# total harvested seed, without assuming any causal transfer.
ratios=[]
for key,v in pairs.items():
 s=num(parents[key]["Selfing_rate"])
 f0=num(v["ID=0"]["Female_fitness"])
 f1=num(v["ID=1"]["Female_fitness"])
 if s is not None and f0 is not None and f0>0 and f1 is not None:
  ratios.append(((f0-f1)/f0, s))
if ratios:
 dif=[x-y for x,y in ratios]
 print("MATERNAL_SELFING_DISCOUNT_AUDIT",json.dumps({"pairs":len(dif),"max_abs_difference":max(abs(x) for x in dif),"median_abs_difference":statistics.median(map(abs,dif)),"mean_signed_difference":statistics.mean(dif)}))
# summary of male-fitness values relative to male flower count, not a model
for mode in ["ID=0","ID=1"]:
 for col in ["Male_fitness","Female_fitness","Total_fitness"]:
  x=[(num(parents[k]["Male_flower_total"]),num(v[mode][col])) for k,v in pairs.items()]
  x=[(a,b) for a,b in x if a is not None and b is not None]
  print("OUTCOME_SUPPORT",json.dumps({"mode":mode,"field":col,"n":len(x),"zeros":sum(y==0 for a,y in x),"max_y":max(y for a,y in x),"mean_y":statistics.mean(y for a,y in x),"nonintegral_y":sum(abs(y-round(y))>1e-8 for a,y in x)}))
