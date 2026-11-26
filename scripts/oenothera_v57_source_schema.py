#!/usr/bin/env python3
"""v57 audit of *independent* Oenothera harringtonii parent + offspring SSR
archive (Rhodes/Fant/Skogen 2017 Mol Ecol; Dryad p24q3, preserved Zenodo
4942351). Cross-species comparator, NOT a Daphne parentage refit.

Read .xlsx with Python stdlib OOXML only; source file never redistributed,
all sample IDs, SSR alleles and individual metadata stay out of output.
"""
import collections,csv,hashlib,io,json,re,urllib.request,zipfile
from xml.etree import ElementTree as ET
BASE="https://zenodo.org/api/records/4942351/files/"
FILES={
  "genotypes.xlsx":"487a6b4e01d30876ea17bb9c31e54bc3",
  "multiplePaternity.csv":"600f6f370ffa8ad205d0ccb6bc92ab65",
  "pollinatorBouts.csv":"d57f7bcbbaf486b01d99bf13e343d3a4",
  "pollenDispersal.csv":"11d75d055360b0d1d4a908019581faa3",
  "Rhodes et al. 2017.R":"a17ce0fb33e9f7f60e94b295b7d539ed"
}
MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
N={"m":MAIN}
def download(name,md5):
    url=BASE+urllib.parse.quote(name)+"/content"
    req=urllib.request.Request(url,headers={"User-Agent":"independent-source-validation-v57/1.0"})
    b=urllib.request.urlopen(req,timeout=90).read()
    if hashlib.md5(b).hexdigest()!=md5:raise RuntimeError(f"Wrong source digest for {name}")
    return b
def col(ix):
    s="";i=ix+1
    while i:s=chr(65+(i-1)%26)+s;i=(i-1)//26
    return s
def spreadsheet(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=z.namelist()
        shared=[]
        if "xl/sharedStrings.xml" in names:
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(x.itertext()) for x in root.findall("m:si",N)]
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=wb.find("m:sheets",N)
        labels=[s.attrib.get("name","") for s in sheets]
        out=[]
        paths=sorted([x for x in names if re.match(r"xl/worksheets/sheet\d+\.xml$",x)],key=lambda p:int(re.search(r"sheet(\d+)",p).group(1)))
        for n,path in enumerate(paths):
            root=ET.fromstring(z.read(path))
            lines=root.findall("m:sheetData/m:row",N)
            rows=[]
            for row in lines:
                vals={}
                for c in row.findall("m:c",N):
                    addr=c.get("r","");mat=re.match("([A-Z]+)",addr)
                    if not mat:continue
                    ci=0
                    for ch in mat.group(1):ci=ci*26+ord(ch)-64
                    ci-=1
                    v=c.find("m:v",N)
                    if c.get("t")=="inlineStr":
                        inline=c.find("m:is",N)
                        value="" if inline is None else "".join(inline.itertext())
                    elif v is None:value=""
                    elif c.get("t")=="s":
                        value=shared[int(v.text)]
                    else:value=v.text or ""
                    vals[ci]=value
                if vals:rows.append(vals)
            maxc=max((max(r) for r in rows if r),default=-1)
            first=rows[0] if rows else {}
            title=[str(first.get(i,"")) for i in range(min(maxc+1,28))]
            # Avoid outputting sample rows or genotype alleles, only first row
            # if human-readable header and only letters/non-IDs.
            title=[re.sub(r"[\d]{5,}","[digits]",x)[:50] for x in title]
            summary={
                "sheet":labels[n] if n<len(labels) else "sheet",
                "xml_rows_nonempty":len(rows),
                "max_column_index":maxc+1,
                "first_row_headers_redacted":title,
                "header_first_cell":title[0] if title else "",
                "rows_numeric_cells":sum(sum(1 for x in r.values() if re.fullmatch(r"\d+(?:\.\d+)?",str(x))) for r in rows),
                "first_10_row_nonempty_cell_counts":[len(r) for r in rows[:10]],
                "rows_available_with_two_plus_cells":sum(len(r)>=2 for r in rows)
            }
            out.append(summary)
    return out
def csv_summary(name,raw):
    text=raw.decode("utf-8-sig")
    text=text.replace("\r\n","\n").replace("\r","\n")
    rr=list(csv.DictReader(io.StringIO(text,newline="")))
    cols=list(rr[0]) if rr else []
    out={"file":name,"rows":len(rr),"columns":cols}
    for k in cols:
        if re.search(r"treatment|visitor|pollinator|time",k,re.I):
            frequency=collections.Counter(str(r.get(k,"")) for r in rr)
            if len(frequency)<30:out[k+"_count"]={a:int(b) for a,b in frequency.items()}
        if re.search(r"correlatedPaternity|isolation20|flowersVisited",k,re.I):
            v=[]
            for r in rr:
                try:v.append(float(r.get(k,"")))
                except (ValueError,TypeError):pass
            if v:out[k+"_summary"]={"count":len(v),"mean":sum(v)/len(v),"min":min(v),"max":max(v)}
    return out
def main():
    out={"source_DOI":"10.5061/dryad.p24q3",
         "source_article":"10.1111/mec.14115",
         "source_zenodo_pinned":"4942351",
         "files":{},"spreadsheet":None,"csv":{}}
    for name,md5 in FILES.items():
        data=download(name,md5)
        out["files"][name]={"md5":md5,"sha256":hashlib.sha256(data).hexdigest(),"bytes":len(data)}
        if name.endswith(".xlsx"):out["spreadsheet"]=spreadsheet(data)
        elif name.endswith(".csv"):out["csv"][name]=csv_summary(name,data)
        elif name.endswith(".R"):
            text=data.decode("utf-8-sig",errors="replace")
            out["original_analysis_features"]={
               x:len(re.findall(x,text,re.I))
               for x in ["MLTR","CERVUS","genotypes","offspring","glmer","lm\\(","correlatedPaternity","isolation20","treatment"]}
    return out
if __name__=="__main__":
    x=main()
    print("OENOTHERA_V57_SOURCE_SCHEMA",json.dumps(x,ensure_ascii=False,allow_nan=False),flush=True)
