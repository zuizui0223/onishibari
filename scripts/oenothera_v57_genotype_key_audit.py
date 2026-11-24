#!/usr/bin/env python3
"""v57 audit of Oenothera 323 parents, 521 offspring, seven SSR genotype pairs.
Parse original XLSX with stdlib XML only. Print strictly aggregates; never print
raw parent/offspring IDs, allele values, or maternal family coordinates.
"""
import collections,hashlib,io,json,re,zipfile,urllib.request,urllib.parse
from xml.etree import ElementTree as ET
from oenothera_v57_source_schema import BASE,MAIN,N,download,FILES
def records(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            tree=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(x.itertext()) for x in tree.findall("m:si",N)]
        out={}
        for label,number in (("parents",1),("offspring",2)):
            root=ET.fromstring(z.read(f"xl/worksheets/sheet{number}.xml"))
            rows=[]
            for row in root.findall("m:sheetData/m:row",N):
                values=[""]*15
                for c in row.findall("m:c",N):
                    addr=c.get("r","");part=re.match(r"([A-Z]+)",addr)
                    if not part:continue
                    idx=0
                    for ch in part.group(1):idx=26*idx+ord(ch)-64
                    idx-=1
                    val=c.find("m:v",N)
                    if c.get("t")=="inlineStr":
                        inline=c.find("m:is",N)
                        v="" if inline is None else "".join(inline.itertext())
                    elif val is None:v=""
                    elif c.get("t")=="s":v=shared[int(val.text)]
                    else:v=str(val.text)
                    if idx<len(values):values[idx]=v
                if values[0] or any(values[1:]):rows.append(values)
            assert len(rows)==(324 if label=="parents" else 522)
            out[label]=rows[1:]
    return out
def looksnumeric(v):
    try:
        x=float(v)
        return x.is_integer()
    except (TypeError,ValueError):return False
def parentage_keys(p,o):
    mothers=set(r[0] for r in p)
    assert len(mothers)==len(p) and all(mothers)
    rawparents=[str(s) for s in mothers]
    offs=[str(r[0]) for r in o]
    exact_parent_prefix=[]
    treatment_suffix=[]
    diagnostics=collections.Counter()
    malformed=[]
    for oid in offs:
        candidates=[p for p in rawparents if oid.startswith(p)]
        if candidates:
            parent=max(candidates,key=len)
            suffix=oid[len(parent):]
            exact_parent_prefix.append((parent,suffix))
            found=re.search(r"^[^a-zA-Z]*([CDN])(?:[^a-zA-Z]*)(\d+)[^a-zA-Z]*$",suffix,re.I)
            if found:
                treatment_suffix.append(found.group(1).upper())
            else:
                treatment_suffix.append("unknown")
                diagnostics["suffix_redacted_"+re.sub("[0-9]","#",re.sub("[a-z]","a",re.sub("[A-Z]","A",suffix)))]+=1
        else:
            malformed.append(oid)
    def mask(z):
        return "".join("#" if ch.isdigit() else "A" if ch.isalpha() else ch for ch in str(z))
    parent_patterns=collections.Counter(mask(v) for v in rawparents)
    offspring_patterns=collections.Counter(mask(v) for v in offs)
    mother_labels=[v.split("_")[0] for v in offs]
    # Text normalization rules are checked only as cohort-level match counts;
    # no original genet or offspring identifiers are emitted.
    from oenothera_v57_source_schema import download,FILES
    import csv
    dry=download("multiplePaternity.csv",FILES["multiplePaternity.csv"]).decode("utf-8-sig")
    dry=dry.replace(chr(13)+chr(10),chr(10)).replace(chr(13),chr(10))
    meta=list(csv.DictReader(io.StringIO(dry,newline="")))
    plantIDs=[str(t["plantID"]) for t in meta]
    def numeric(x):return "".join(ch for ch in str(x) if ch.isdigit())
    knownp=set(rawparents)
    plant_parent_straight=sum(x in knownp for x in plantIDs)
    plant_off_mother_straight=sum(x in set(mother_labels) for x in plantIDs)
    match_off_to_parent_digits=sum(any(numeric(y)==numeric(x) for y in rawparents) for x in mother_labels)
    match_off_to_parent_first4=sum(any(y[:4]==x for y in rawparents) for x in mother_labels)
    match_off_to_parent_last4=sum(any(y[-4:]==x for y in rawparents) for x in mother_labels)
    parent_ids=collections.Counter(map(len,rawparents))
    output={
       "original_parent_rows":len(p),"unique_parent_ids":len(mothers),
       "parent_id_shape_counts":dict(parent_patterns),
       "offspring_id_shape_counts":dict(offspring_patterns),
       "offspring_maternal_prefix_unique":len(set(mother_labels)),
       "offspring_maternal_prefix_length_counts":dict(collections.Counter(len(x) for x in mother_labels)),
       "dryad_plant_id_shape_counts":dict(collections.Counter(mask(x) for x in plantIDs)),
       "dryad_plantID_exact_in_parent":plant_parent_straight,
       "dryad_plantID_exact_in_offspring_mother_prefix":plant_off_mother_straight,
       "offspring_maternal_match_by_digit_only":match_off_to_parent_digits,
       "offspring_maternal_match_to_parent_first4":match_off_to_parent_first4,
       "offspring_maternal_match_to_parent_last4":match_off_to_parent_last4,
       "original_offspring_rows":len(o),
       "parent_id_length_counts":dict(parent_ids),
       "offspring_exact_parent_prefix_matches":len(exact_parent_prefix),
       "offspring_unmatched_count":len(malformed),
       "offspring_unmatched_mask_patterns":dict(collections.Counter(re.sub("[0-9]","#",re.sub("[a-z]","a",re.sub("[A-Z]","A",x))) for x in malformed).most_common(10)),
       "treatment_from_matched_suffix_counts":dict(collections.Counter(treatment_suffix)),
       "matched_suffix_shape_frequencies":dict(diagnostics),
       "parent_genotype_nonempty_allele_cells":sum(sum(bool(x) for x in row[1:15]) for row in p),
       "offspring_genotype_nonempty_allele_cells":sum(sum(bool(x) for x in row[1:15]) for row in o),
       "parent_genotype_nonnumeric_nonempty_cells":sum(sum(bool(x) and not looksnumeric(x) for x in row[1:15]) for row in p),
       "offspring_genotype_nonnumeric_nonempty_cells":sum(sum(bool(x) and not looksnumeric(x) for x in row[1:15]) for row in o),
       "offspring_allele_missing_distribution":dict(collections.Counter(sum(not bool(x) for x in row[1:15]) for row in o)),
       "parent_allele_missing_distribution":dict(collections.Counter(sum(not bool(x) for x in row[1:15]) for row in p)),
       "offspring_unique_ids":len(set(offs))
    }
    return output
if __name__=="__main__":
    data=download("genotypes.xlsx",FILES["genotypes.xlsx"])
    rows=records(data)
    print("OENOTHERA_V57_ALLELE_KEY_AUDIT",json.dumps(parentage_keys(rows["parents"],rows["offspring"]),ensure_ascii=False,allow_nan=False),flush=True)
