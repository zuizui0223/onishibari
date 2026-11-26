#!/usr/bin/env python3
"""Public archive column and original R-analysis glossary audit, v0.52.
Never print original maternal/offspring IDs, exact coordinates, or genotype rows.
Only file names, column names, available n, and nonidentifying R expressions.
This script is a SCHEMA DISCOVERY stage, not a fit to fitness data.
"""
import csv, hashlib, io, json, re, urllib.request, zipfile
URL="https://zenodo.org/api/records/15644955/files/Geito_NewPhyto_v2.zip/content"
SHA="a18b7d296110be1e9817d6e1a3994365a4e7ae0aab50435e7e6883e1f3de4042"
PREFIX="Geito_NewPhyto_v2/"
raw=urllib.request.urlopen(urllib.request.Request(URL,headers={"User-Agent":"Primary-archive-audit-v0.52"}),timeout=90).read()
assert hashlib.sha256(raw).hexdigest()==SHA
z=zipfile.ZipFile(io.BytesIO(raw))
for file in ["3G.csv","3G_ID01.csv","Genotype_3P.csv"]:
    r=list(csv.DictReader(io.StringIO(z.read(PREFIX+file).decode("utf-8-sig"))))
    cols=list(r[0]) if r else []
    print("SOURCE_SCHEMA",json.dumps({"file":file,"n_rows":len(r),"n_columns":len(cols),"columns":cols},ensure_ascii=False))
    if file=="3G.csv":
        for key in cols:
            if re.search(r"^ID$|^Population$|^Sex$|father|mother|patern|outcross|self|fitness|gender|male|female|flower|biomass|seed|offspring",key,re.I):
                vals=[t[key] for t in r if t.get(key) not in ("","NA","NaN")]
                numeric=[]
                for v in vals:
                    try:numeric.append(float(v))
                    except (TypeError,ValueError):pass
                d={"column":key,"nonmissing":len(vals),"numeric_n":len(numeric)}
                if numeric:d.update(min=min(numeric),max=max(numeric),mean=sum(numeric)/len(numeric))
                print("COLUMN_SUMMARY",json.dumps(d,ensure_ascii=False))
for file in ["Selfing.R","Fitness.R","Trade-off.R","README.md.txt"]:
    text=z.read(PREFIX+file).decode("utf-8-sig",errors="replace")
    matches=[]
    for line in text.splitlines():
        line=line.strip()
        if re.search(r"Selfed|Genotyped|Male_|Female_|Outcross|Sired|Father|Fathering|Paternity|offspring|fitness|pollen|RS_|δ|delta",line,re.I) and len(line)<250 and not line.startswith("#"):
            matches.append(line)
    print("R_EXPRESSION_AUDIT",json.dumps({"file":file,"matched_expressions":matches[:35],"total_matches":len(matches)},ensure_ascii=False))
