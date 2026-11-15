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
    matching_ids=set()
    for name in NAMES:
        query=f'"{name}"[Organism]'
        u=url_for("esearch.fcgi",dict(db="nuccore",term=query,
                                         retmode="json",retmax=200))
        q=dict(taxon_query=name,query=query,url=u)
        try:
            content,_=fetch(u)
            result=json.loads(content)["esearchresult"]
            n=int(result["count"]); ids=result.get("idlist",[])
            if n>len(ids):raise ValueError("NCBI result truncated; retry with pagination")
            matching_ids.update(ids)
            q.update(status="success",reported_record_count=n,returned_ids=len(ids))
        except Exception as e:q.update(status="error",error=str(e))
        synonym_search.append(q)
        time.sleep(0.4)
    discovered=[]
    classification_errors=[]
    plastid_expected={"LC680075":"matK","LC693863":"rbcL",
                      "OK419833":"trnH-psbA"}
    plastid_verified={}
    # An ESearch synonym match can carry the wrong modern organism name.
    # Fetch records and check original registered taxon before classifying.
    for uid_batch_start in range(0,len(matching_ids),50):
        batch=sorted(matching_ids)[uid_batch_start:uid_batch_start+50]
        try:
            records,source=genbank(batch)
            for record in records:
                organism=record.annotations.get("organism","")
                description=record.description
                locus_kind=("ITS" if looks_like_its(description) else
                            "published_2014_SSR" if record.id.split(".")[0] in EXPECTED
                            else "unresolved_other")
                acc=record.id.split(".")[0]
                if acc in plastid_expected and accepted_taxon(organism):
                    marker=plastid_expected[acc]
                    desc_norm=description.lower()
                    if marker.lower() not in desc_norm and not (
                            marker=="trnH-psbA" and "trnh-psba" in desc_norm):
                        raise ValueError("Expected plastid region not found: "+record.id)
                    seq=str(record.seq).upper()
                    if set(seq)-set("ACGTRYSWKMBDHVN"):
                        raise ValueError("Invalid plastid DNA alphabet")
                    plastid_verified[acc]=(record.id,marker,seq)
                discovered.append(dict(accession=record.id,organism=organism,
                    description=description,locus_kind=locus_kind,
                    length_bp=len(record.seq),matches_expected_taxon=accepted_taxon(organism),
                    source_url=source))
        except Exception as e:
            classification_errors.append(str(e))
    with (args.out/"NCBI_jezoensis_record_types.csv").open("w",encoding="utf-8",newline="") as fh:
        keys=("accession","organism","description","locus_kind","length_bp",
              "matches_expected_taxon","source_url")
        w=csv.DictWriter(fh,fieldnames=keys);w.writeheader();w.writerows(discovered)
    with (args.out/"jezoensis_3_plastid_fragments.fasta").open("w") as fh:
        for acc in plastid_expected:
            if acc not in plastid_verified:continue
            ident,marker,seq=plastid_verified[acc]
            fh.write(">"+ident+"|"+marker+"|Daphne_jezoensis|nakata074\n")
            fh.write("\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n")
    result=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        primary_source="https://doi.org/10.3732/apps.1400001",
        original_published_ssr=10,verified_ssr=len(verified),
        ssr_errors=failed,synonym_search=synonym_search,
        ncbi_record_types=discovered,
        ncbi_organism_matched_count=sum(x["matches_expected_taxon"] for x in discovered),
        ncbi_ITS_matched_count=sum(x["matches_expected_taxon"] and x["locus_kind"]=="ITS" for x in discovered),
        classification_errors=classification_errors,
        plastid_verified={k:dict(accession=v[0],marker=v[1],length_bp=len(v[2]),
            sha256=hashlib.sha256(v[2].encode()).hexdigest())
            for k,v in plastid_verified.items()},
        note="SSR fragments are not ITS; cannot add to 24-tip ITS phylogeny",
        target_cross_amplification="NOT_TESTED",
        sex_system_of_original_target="GYNODIOECY",
        original_species_experiment="47 individuals, two D. jezoensis populations")
    (args.out/"summary.json").write_text(
         json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("JEZOENSIS_SSR_AUDIT="+json.dumps(dict(
        verified_ssr=len(verified),failed=len(failed),
        locus_types=[(x["accession"],x["locus_kind"]) for x in discovered],
        queried_ITS_count=result["ncbi_ITS_matched_count"],
        plastid_markers_verified=len(plastid_verified),
        query_counts=[(q["taxon_query"],q.get("reported_record_count"),
                       q["status"]) for q in synonym_search]),ensure_ascii=False),flush=True)
    if not verified:raise SystemExit("No independently verified SSR accession")

if __name__=="__main__":main()
