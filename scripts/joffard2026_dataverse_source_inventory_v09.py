#!/usr/bin/env python3
"""Audit public 2026 Joffard et al. Silene dioica paternal + pollen carryover data.

Read-only official Dataverse API without token/cookies, downloads only two
unrestricted tabular files. Never export mother, flower, bee or father IDs.
Source: https://doi.org/10.57745/CF1VNW
Metadata is the authority for file IDs and published version.
"""
import csv,hashlib,io,json,urllib.request,urllib.parse,urllib.error,zipfile
from pathlib import Path
from collections import Counter

OUT=Path("analysis_results_joffard2026_public_source_v09")
OUT.mkdir(exist_ok=True)
HOST="https://entrepot.recherche.data.gouv.fr"
DOI="doi:10.57745/CF1VNW"
FILES={"paternity_assignment_results.tab","pollen_counts_on_stigmas.tab"}
CAP=600000

def read(url,accept):
    req=urllib.request.Request(url,headers={"User-Agent":"Silene-pollen-transfer-original-data-audit/0.1","Accept":accept})
    try:
        with urllib.request.urlopen(req,timeout=45) as r:
            blob=r.read(CAP+1)
            if len(blob)>CAP:raise RuntimeError("response over safety cap")
            return blob,{"http_status":r.status,"content_type":r.headers.get("content-type"),"url":r.geturl()}
    except urllib.error.HTTPError as err:
        return None,{"http_status":err.code,"error":"public_http_error"}
    except Exception as err:
        return None,{"error_type":type(err).__name__}
def inventory(data):
    txt=data.decode("utf-8-sig","replace").replace("\r\n","\n").replace("\r","\n")
    typ=csv.excel_tab if "\t" in txt.partition("\n")[0] else csv.excel
    rr=csv.DictReader(io.StringIO(txt),dialect=typ)
    rows=list(rr)
    cols=rr.fieldnames or []
    out={"n_records":len(rows),"columns":cols,"n_missing":{},"distinct_counts":{},"small_categorical_frequency":{},
         "numeric_summaries":{}}
    for k in cols:
        v=[r.get(k,"") for r in rows]
        miss=sum(x is None or not str(x).strip() or x in ("NA","NaN","nan") for x in v)
        out["n_missing"][k]=miss
        vals=[str(x).strip() for x in v if x is not None and str(x).strip() not in ("","NA","NaN","nan")]
        out["distinct_counts"][k]=len(set(vals))
        if len(set(vals))<=20:out["small_categorical_frequency"][k]=dict(Counter(vals))
        numeric=[]
        for x in vals:
            try:numeric.append(float(x))
            except ValueError:pass
        if numeric and len(numeric)>=.9*len(vals):
            out["numeric_summaries"][k]={"min":min(numeric),"max":max(numeric),
                                       "sum":sum(numeric),"n_numeric":len(numeric),
                                       "n_zero":sum(x==0 for x in numeric)}
    return out

def main():
    meta_url=HOST+"/api/datasets/:persistentId/?"+urllib.parse.urlencode({"persistentId":DOI})
    data,stat=read(meta_url,"application/json")
    report={"status":"failed","dataset":DOI,"metadata_probe":stat,"files":[],
      "data_consent":"only openly available deposited files; no authorization or access request",
      "source_rows_never_republished":True}
    if not data:
        report["status"]="metadata_unavailable"
        (OUT/"joffard2026_metadata_inventory.json").write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps(report,indent=2))
        return
    doc=json.loads(data)
    if doc.get("status")!="OK":raise ValueError("unexpected Dataverse response")
    vers=doc["data"]["latestVersion"]
    items=vers.get("files",[])
    report["datasetVersion"]={"major":vers.get("versionNumber"),"minor":vers.get("versionMinorNumber"),
                              "releaseTime":vers.get("releaseTime"),"files":len(items)}
    for item in items:
        f=item["dataFile"]
        label=item.get("label","")
        if label not in FILES:continue
        fid=f["id"]
        ent={"name":label,"file_id":fid,"source_bytes_declared":f.get("filesize"),
             "original_md5":f.get("md5"),"metadata_variables":f.get("variables"),
             "restricted":f.get("restricted")}
        url=HOST+"/api/access/datafile/"+str(fid)+"?"+urllib.parse.urlencode({"format":"original"})
        raw,access=read(url,"text/plain, text/tab-separated-values, application/octet-stream")
        ent["original_access"]=access
        if raw is not None:
            ent["bytes"]=len(raw)
            ent["sha256"]=hashlib.sha256(raw).hexdigest()
            if raw[:2]==b"PK" and zipfile.is_zipfile(io.BytesIO(raw)):
                ent["error"]="unexpected archive instead of table"
            else:ent["table"]=inventory(raw)
        report["files"].append(ent)
        print("PUBLIC_ORIGINAL_FILE",json.dumps(ent,ensure_ascii=False)[:13000],flush=True)
    if len(report["files"])!=2:raise ValueError("expected two Joffard study source files; got "+str(len(report["files"])))
    report["status"]="all_original_tables_loaded" if all("table" in f for f in report["files"]) else "partial_or_inaccessible_source"
    (OUT/"joffard2026_metadata_inventory.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("DATAVERSE_ORIGINAL_AUDIT",report["status"],flush=True)

if __name__=="__main__":main()
