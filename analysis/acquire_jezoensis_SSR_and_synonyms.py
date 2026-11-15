#!/usr/bin/env python3
"""Acquire *Daphne jezoensis* published SSR fragments and audit synonym records.

Not a species-tree analysis: AB897817-AB897826 are 2014 microsatellite
fragments, NOT nuclear ITS; transferability to *D. pseudomezereum* untested.
"""
import argparse,csv,hashlib,io,json,re,time
from datetime import datetime,timezone
from pathlib import Path
import urllib.parse,urllib.request
from Bio import SeqIO

SOURCE=Path(__file__).resolve().parents[1]/"data/derived/jezoensis_2014_SSR_primer_candidates.csv"
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NAMES=(
 "Daphne jezoensis",
 "Daphne kamtschatica var. jezoensis",
 "Daphne kamtschatica subsp. jezoensis",
 "Daphne pseudomezereum var. jezoensis",
 "Daphne pseudomezereum subsp. jezoensis",
 "Daphne rebunensis",
)
EXPECTED=[f"AB8978{x}" for x in range(17,27)]

def accepted_taxon(name):
    return name.strip().casefold() in {x.casefold() for x in NAMES}

def looks_like_its(title):
    return bool(re.search(r"internal\s+transcribed\s+spacer|\bITS[12]?\b",title,re.I))

def manifest():
    with SOURCE.open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    if sorted(x["accession"] for x in rows)!=sorted(EXPECTED):
        raise ValueError("Unexpected accession inventory")
    if len({x["locus"] for x in rows})!=10:
        raise ValueError("Duplicate locus")
    return rows

def fetch(url):
    err=None
    for attempt in range(3):
        try:
            request=urllib.request.Request(url,headers={
                "User-Agent":"Daphne-jezoensis-SSR-audit/2026"})
            with urllib.request.urlopen(request,timeout=45) as stream:
                return stream.read(12000000).decode("utf-8"),stream.geturl()
        except Exception as e:
            err=e
            time.sleep(2**attempt)
    raise RuntimeError(f"{type(err).__name__}: {str(err)[:180]}")

def url_for(name,parameters):
    return BASE+name+"?"+urllib.parse.urlencode(parameters)

def genbank(ids):
    u=url_for("efetch.fcgi",dict(db="nuccore",id=",".join(ids),
                                   rettype="gb",retmode="text"))
    data,u=fetch(u)
    return list(SeqIO.parse(io.StringIO(data),"genbank")),u

def selftest():
    assert len(manifest())==10
    assert accepted_taxon("Daphne kamtschatica var. jezoensis")
    assert not accepted_taxon("Daphne pseudomezereum")
    assert not accepted_taxon("Daphne pseudomezereum var. koreana")
    assert not looks_like_its("Daphne jezoensis microsatellite Dp171 genomic sequence")
    assert looks_like_its("internal transcribed spacer ITS1")
    print("JEZOENSIS_SSR_AUDIT_SELFTEST=SUCCESS: 10 genomic SSR markers, not ITS")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=Path("outputs/jezoensis_nuclear_ssr_20261115"))
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args();selftest()
    if args.self_test:return
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Snapshot output directory must be empty")
    args.out.mkdir(parents=True,exist_ok=True)
    entries,origin=genbank(EXPECTED)
    by={r.id.split(".")[0]:r for r in entries}
    verified=[];failed=[]
    for original in manifest():
        acc=original["accession"]
        record=by.get(acc)
        if record is None:
            failed.append(dict(accession=acc,error="missing_accession"));continue
        organism=record.annotations.get("organism","")
        seq=str(record.seq).upper()
        if (not accepted_taxon(organism) or looks_like_its(record.description)
            or not 70<=len(seq)<=5000 or set(seq)-set("ACGTRYSWKMBDHVN")):
            failed.append(dict(accession=acc,error="taxon_or_marker_mismatch",
                               organism=organism,description=record.description))
            continue
        verified.append(dict(accession=record.id,source_paper=original["source_paper"],
             organism=organism,locus=original["locus"],length_bp=len(seq),
             sha256=hashlib.sha256(seq.encode()).hexdigest(),
             record_description=record.description,source_url=origin,sequence=seq))
    with (args.out/"verified_jezoensis_ssr.fasta").open("w") as f:
        for r in verified:
            f.write(">"+r["accession"]+"|"+r["locus"]+"|SSR_NOT_ITS\n")
            s=r["sequence"]
            f.write("\n".join(s[i:i+80] for i in range(0,len(s),80))+"\n")
    with (args.out/"verified_ssr_records.csv").open("w",newline="",encoding="utf-8") as f:
        keys=("accession","organism","locus","length_bp","sha256",
              "record_description","source_paper","source_url")
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for x in verified:w.writerow({k:x[k] for k in keys})
    synonym_search=[]
    for name in NAMES:
        query=f'"{name}"[Organism]'
        u=url_for("esearch.fcgi",dict(db="nuccore",term=query,
                                         retmode="json",retmax=0))
        q=dict(taxon_query=name,query=query,url=u)
        try:
            content,_=fetch(u)
            n=int(json.loads(content)["esearchresult"]["count"])
            q.update(status="success",reported_record_count=n)
        except Exception as e:q.update(status="error",error=str(e))
        synonym_search.append(q)
        time.sleep(0.4)
    result=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        primary_source="https://doi.org/10.3732/apps.1400001",
        original_published_ssr=10,verified_ssr=len(verified),
        ssr_errors=failed,synonym_search=synonym_search,
        note="SSR fragments are not ITS; cannot add to 24-tip ITS phylogeny",
        target_cross_amplification="NOT_TESTED",
        sex_system_of_original_target="GYNODIOECY",
        original_species_experiment="47 individuals, two D. jezoensis populations")
    (args.out/"summary.json").write_text(
         json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("JEZOENSIS_SSR_AUDIT="+json.dumps(dict(
        verified_ssr=len(verified),failed=len(failed),
        query_counts=[(q["taxon_query"],q.get("reported_record_count"),
                       q["status"]) for q in synonym_search]),ensure_ascii=False),flush=True)
    if not verified:raise SystemExit("No independently verified SSR accession")

if __name__=="__main__":main()
