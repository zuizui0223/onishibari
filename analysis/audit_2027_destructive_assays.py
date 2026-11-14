#!/usr/bin/env python3
"""Audit separately sampled stigma, pollen-tube, and embryo assay cohorts.

This does NOT infer self-incompatibility, paternity, or transitions between
different flowers. Every input/output must live under a private directory.
No actual Onishibari samples are included; --self-test uses synthetic fixtures.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import date
from pathlib import Path

from randomize_2027_pilot import MORPHS

KEY = ("population_id", "genet_id", "branch_id", "flower_id")
HEADERS = KEY + (
    "sex_morph", "treatment", "assay_stage", "assay_status",
    "assay_date", "hours_since_pollination", "pollen_qc",
    "pollen_donor_genet_id", "pollen_grains_count",
    "pollen_grains_germinated", "pollen_tubes_at_ovary",
    "ovules_examined", "fertilized_ovules", "intact_embryos",
    "aborted_embryos", "specimen_evidence_id", "notes",
)
STAGES = ("stigma", "pollen_tube", "embryo")
SUMMARY = ("sex_morph", "treatment", "assay_stage", "entries",
           "sampled_valid", "genets_sampled_valid", "lost_pending",
           "failed_pollen_or_missing_evidence", "sum_grains", "sum_germinated",
           "sum_tubes_ovary", "sum_ovules", "sum_fertilized_ovules",
           "sum_intact_embryos", "sum_aborted_embryos", "interpretation")
METRICS = (
    "pollen_grains_count", "pollen_grains_germinated", "pollen_tubes_at_ovary",
    "ovules_examined", "fertilized_ovules", "intact_embryos", "aborted_embryos"
)
ALLOWED_TREATMENTS = set.union(*(set(xs) for xs in MORPHS.values()))
ALLOWED_STATUSES = ("sampled", "pending", "lost")

def parse(rows):
    seen=set()
    results=[]
    for row in rows:
        key=tuple(row.get(k,"").strip() for k in KEY)
        if any(not x for x in key) or key in seen:
            raise ValueError("Missing or duplicate dissection cohort flower ID")
        seen.add(key)
        stage=row.get("assay_stage")
        status=row.get("assay_status")
        morph=row.get("sex_morph")
        treatment=row.get("treatment")
        if stage not in STAGES or status not in ALLOWED_STATUSES:
            raise ValueError("Invalid assay_stage/assay_status")
        if morph not in MORPHS or treatment not in MORPHS[morph]:
            raise ValueError("Unverified morph/treatment; cannot score")
        if status=="sampled":
            if not row.get("assay_date") or not row.get("specimen_evidence_id"):
                raise ValueError("Sampled tissue needs date and specimen evidence ID")
            date.fromisoformat(row["assay_date"])
            if not row.get("hours_since_pollination"):
                raise ValueError("Specify hours since pollination")
            if float(row["hours_since_pollination"])<0:
                raise ValueError("Negative assay time")
        else:
            if any(row.get(x,"").strip() not in ("","NA") for x in METRICS):
                raise ValueError("Pending/lost samples cannot have scored counts")
        if treatment in ("bag_hand_self","bag_hand_cross","open_plus_cross") and status=="sampled":
            donor=row.get("pollen_donor_genet_id","").strip()
            if not donor:
                raise ValueError("Hand pollination requires identified donor")
            if treatment=="bag_hand_self" and donor!=row["genet_id"]:
                raise ValueError("Self treatment donor differs from maternal genet")
            if treatment!="bag_hand_self" and donor==row["genet_id"]:
                raise ValueError("Cross treatment donor equals maternal genet")
        quality=row.get("pollen_qc")
        if quality not in ("passed","failed","not_checked","not_applicable"):
            raise ValueError("Invalid pollen QC value")
        values={}
        if status=="sampled":
            for metric in METRICS:
                raw=row.get(metric,"").strip()
                if raw in ("","NA"):
                    values[metric]=None
                else:
                    try:
                        count=int(raw)
                    except ValueError:
                        raise ValueError("Non-integer assay count: "+metric)
                    if count<0:
                        raise ValueError("Negative assay count")
                    values[metric]=count
            def require(what):
                if any(values[n] is None for n in what):
                    raise ValueError("Sampled stage missing required counts "+stage)
            if stage=="stigma":
                require(("pollen_grains_count","pollen_grains_germinated"))
                if values["pollen_grains_germinated"]>values["pollen_grains_count"]:
                    raise ValueError("Germinated grains exceed total")
            elif stage=="pollen_tube":
                require(("pollen_grains_count","pollen_tubes_at_ovary"))
                if values["pollen_tubes_at_ovary"]>values["pollen_grains_count"]:
                    raise ValueError("Tubes exceed observed grains")
            else:
                require(("ovules_examined","fertilized_ovules",
                         "intact_embryos","aborted_embryos"))
                if not (values["intact_embryos"]+values["aborted_embryos"]
                        <=values["fertilized_ovules"]<=values["ovules_examined"]):
                    raise ValueError("Impossible embryo counts")
            allowed={
                "stigma":{"pollen_grains_count","pollen_grains_germinated"},
                "pollen_tube":{"pollen_grains_count","pollen_tubes_at_ovary"},
                "embryo":{"ovules_examined","fertilized_ovules",
                          "intact_embryos","aborted_embryos"},
            }[stage]
            if any(values[n] is not None for n in METRICS if n not in allowed):
                raise ValueError("Metrics from different dissection stages mixed")
        valid=(status=="sampled" and bool(row.get("specimen_evidence_id")) and (
            treatment not in ("bag_hand_cross","bag_hand_self","open_plus_cross")
            or quality=="passed"))
        results.append({**row, "_valid":valid, "_values":values})
    return results

def summarise(data):
    groups=defaultdict(list)
    for r in data:
        groups[(r["sex_morph"],r["treatment"],r["assay_stage"])].append(r)
    out=[]
    for key,records in sorted(groups.items()):
        good=[r for r in records if r["_valid"]]
        def total(metric):
            return sum(r["_values"].get(metric) or 0 for r in good)
        out.append({
            "sex_morph":key[0],"treatment":key[1],"assay_stage":key[2],
            "entries":len(records),"sampled_valid":len(good),
            "genets_sampled_valid":len({(r["population_id"],r["genet_id"]) for r in good}),
            "lost_pending":sum(r["assay_status"] in ("pending","lost") for r in records),
            "failed_pollen_or_missing_evidence":sum(
                r["assay_status"]=="sampled" and not r["_valid"] for r in records),
            "sum_grains":total("pollen_grains_count"),
            "sum_germinated":total("pollen_grains_germinated"),
            "sum_tubes_ovary":total("pollen_tubes_at_ovary"),
            "sum_ovules":total("ovules_examined"),
            "sum_fertilized_ovules":total("fertilized_ovules"),
            "sum_intact_embryos":total("intact_embryos"),
            "sum_aborted_embryos":total("aborted_embryos"),
            "interpretation":"DIFFERENT_FLOWERS_NOT_STAGE_TRANSITION_CAUSALITY",
        })
    return out

def self_test():
    common=dict(population_id="SYNTHETIC",genet_id="H",branch_id="B",
                sex_morph="male_biased_morph",treatment="bag_hand_cross",
                assay_status="sampled",assay_date="2027-02-22",
                hours_since_pollination="24",
                pollen_qc="passed",pollen_donor_genet_id="DONOR_X",
                specimen_evidence_id="SYNTHETIC-VOUCHER",notes="SYNTHETIC")
    cases=[
        dict(common,flower_id="F01",assay_stage="stigma",
             pollen_grains_count="24",pollen_grains_germinated="12"),
        dict(common,flower_id="F02",assay_stage="pollen_tube",
             pollen_grains_count="24",pollen_tubes_at_ovary="6"),
        dict(common,flower_id="F03",assay_stage="embryo",
             ovules_examined="2",fertilized_ovules="1",
             intact_embryos="0",aborted_embryos="1"),
        dict(common,flower_id="F04",assay_stage="embryo",
             assay_status="lost",assay_date="",
             specimen_evidence_id="",hours_since_pollination=""),
    ]
    out=summarise(parse(cases))
    assert len(out)==3
    assert sum(r["sampled_valid"] for r in out)==3
    assert sum(r["lost_pending"] for r in out)==1
    corrupt=[dict(x) for x in cases]
    corrupt[0]["pollen_grains_germinated"]="25"
    try:
        parse(corrupt)
    except ValueError:
        pass
    else:
        raise AssertionError("Impossible pollen count was accepted")
    corrupt=[dict(x) for x in cases]
    corrupt[1]["flower_id"]="F01"
    try:
        parse(corrupt)
    except ValueError:
        pass
    else:
        raise AssertionError("Same flower used for two destructive stages")
    corrupt=[dict(x) for x in cases]
    corrupt[2]["intact_embryos"]="2"
    try:
        parse(corrupt)
    except ValueError:
        pass
    else:
        raise AssertionError("Impossible embryo progression was accepted")
    print("STAGE_ASSAY_QC_TEST=SUCCESS; SYNTHETIC destructive flowers only",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--assays",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        self_test();return
    if not args.assays or not args.output:
        p.error("Both --assays and --output required")
    if not all("private" in str(x).split("/") for x in (args.assays,args.output)):
        p.error("Raw assay records and per-genet outputs must stay under private/")
    with args.assays.open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or not set(HEADERS).issubset(reader.fieldnames):
            raise ValueError("Missing required fields in dissection sample CSV")
        rows=list(reader)
    if not rows:
        raise ValueError("No real assay rows; do not substitute synthetic results")
    summary=summarise(parse(rows))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=SUMMARY)
        writer.writeheader()
        writer.writerows(summary)
    print("STAGE_ASSAY_AUDIT=AGGREGATED; no fertility mechanism inferred",flush=True)
    print("Rows:",len(rows),"Grouped stages:",len(summary),flush=True)

if __name__=="__main__":
    main()
