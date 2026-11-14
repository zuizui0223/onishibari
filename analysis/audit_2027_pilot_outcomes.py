#!/usr/bin/env python3
"""Audit execution fidelity and stage-specific endpoints of the 2027 Onishibari pilot.

No network and no artificial observations in production.
Inputs are PRIVATE field records. Outputs contain genet and flower IDs and
must likewise stay in PRIVATE storage; public GitHub stores only source code.
Use --self-test to exercise synthetic fixtures which are NEVER field evidence.

The script does NOT claim a fertility mechanism; its role is to prevent false
positives and false negatives from design drift, missingness and pseudoreplication.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import date
from pathlib import Path
import tempfile

from randomize_2027_pilot import MORPHS, allocate, COLUMNS

ID = ("population_id", "genet_id", "branch_id", "flower_id")
EXEC_HEADER = ID + (
    "actual_treatment", "execution_status", "execution_date",
    "bag_integrity", "pollen_viability_qc", "cross_donor_genet_id",
    "self_donor_genet_id", "operator_notes",
)
FRUIT_HEADER = ID + (
    "assessment_date", "followup_status", "fruit_mature", "mature_seed_count",
    "seed_germinated_count", "observer_notes",
)
VALID_EXECUTION = {"performed", "missed", "deviation"}
VALID_BAG = {"intact", "breached", "unknown", "not_applicable"}
VALID_POLLEN = {"passed", "failed", "not_checked", "not_applicable"}
VALID_FOLLOWUP = {"completed", "lost", "pending", "destructive_sampling"}
HAND_CROSS = {"bag_hand_cross", "open_plus_cross"}
HAND_SELF = {"bag_hand_self"}
BAG = {"bag_no_hand", "bag_sham", "bag_hand_self", "bag_hand_cross"}
REPORT_HEADER = ("sex_morph","treatment","planned_flowers","genets_assigned",
                 "performed_as_assigned","protocol_deviations",
                 "missing_execution","lost_or_pending",
                 "completed_itt_flowers","completed_pp_flowers",
                 "pp_mature_fruit","pp_seed_bearing_flowers",
                 "pp_total_mature_seeds","pp_germinated_seeds",
                 "pp_genets_with_followup","pp_genets_with_seed",
                 "gate_status")

def load(path, required, allow_empty=False):
    with Path(path).open(encoding="utf-8",newline="") as stream:
        rdr=csv.DictReader(stream)
        if rdr.fieldnames is None or not set(required).issubset(rdr.fieldnames):
            raise ValueError("Invalid header in "+str(path))
        rows=[{x:(str(v).strip() if v is not None else "")
               for x,v in raw.items() if x is not None} for raw in rdr]
    if not rows and not allow_empty:
        raise ValueError("No rows in "+str(path)+"; no observations will be fabricated")
    return rows

def key(row):
    k=tuple(row.get(x,"") for x in ID)
    if any(not value for value in k):
        raise ValueError("Empty/missing population/genet/branch/flower ID")
    return k

def index(rows, source):
    out={}
    for r in rows:
        k=key(r)
        if k in out:
            raise ValueError(f"Duplicate observation ID in {source}: {k}")
        out[k]=r
    return out

def is_int(x):
    try:
        n=int(x)
    except (TypeError,ValueError):
        raise ValueError("Malformed seed count: "+str(x))
    if n<0:
        raise ValueError("Negative seed count")
    return n

def valid_date(s,label):
    if not s:
        raise ValueError("Missing "+label+" date for completed action")
    return date.fromisoformat(s)

def execution_quality(row, planned):
    if row is None:
        return ["execution_unrecorded"], False
    if row.get("execution_status") not in VALID_EXECUTION:
        raise ValueError("Unexpected execution_status: "+str(row.get("execution_status")))
    if row.get("bag_integrity") not in VALID_BAG or row.get("pollen_viability_qc") not in VALID_POLLEN:
        raise ValueError("Missing or unknown bag/pollen quality category")
    status=row["execution_status"]
    if status=="missed":
        if row.get("actual_treatment"):
            raise ValueError("Missed action cannot have an actual treatment")
        return ["assigned_but_not_executed"],False
    valid_date(row.get("execution_date"),"execution")
    notes=[]
    if status=="deviation" or row.get("actual_treatment")!=planned:
        notes.append("different_treatment_or_deviation")
    if planned in BAG and row["bag_integrity"]!="intact":
        notes.append("bag_not_confirmed_intact")
    if planned not in BAG and row["bag_integrity"] not in ("not_applicable","intact"):
        notes.append("open_control_bag_status_unexpected")
    if planned in HAND_CROSS|HAND_SELF and row["pollen_viability_qc"]!="passed":
        notes.append("pollen_viability_not_passed")
    if planned in HAND_CROSS:
        donor=row.get("cross_donor_genet_id","")
        if not donor or donor==row["genet_id"]:
            notes.append("cross_donor_missing_or_same_genet")
    if planned in HAND_SELF:
        donor=row.get("self_donor_genet_id","")
        if not donor or donor!=row["genet_id"]:
            notes.append("self_pollen_donor_unverified")
    if planned not in HAND_CROSS and row.get("cross_donor_genet_id"):
        notes.append("unexpected_cross_donor")
    if planned not in HAND_SELF and row.get("self_donor_genet_id"):
        notes.append("unexpected_self_donor")
    if planned not in HAND_CROSS|HAND_SELF and row["pollen_viability_qc"]=="passed":
        notes.append("pollen_qc_without_hand_pollination")
    return notes, len(notes)==0

def fruit_quality(row):
    if row is None:
        return ["followup_unrecorded"],False,None
    if row.get("followup_status") not in VALID_FOLLOWUP:
        raise ValueError("Unknown followup_status: "+str(row.get("followup_status")))
    if row["followup_status"] != "completed":
        if any(row.get(f) not in ("", "NA") for f in
               ("fruit_mature","mature_seed_count","seed_germinated_count")):
            raise ValueError("Incomplete followup has numerical fruit/seed outcome")
        return ["outcome_pending_or_lost"],False,None
    valid_date(row.get("assessment_date"),"outcome assessment")
    if row.get("fruit_mature") not in ("0","1"):
        raise ValueError("Completed fruit followup must be explicitly 0 or 1")
    fruit=int(row["fruit_mature"])
    seeds=is_int(row.get("mature_seed_count"))
    if not fruit and seeds:
        raise ValueError("Mature seeds are not compatible with fruit_mature=0")
    if row.get("seed_germinated_count") in ("", "NA"):
        germinated=None
    else:
        germinated=is_int(row["seed_germinated_count"])
        if germinated>seeds:
            raise ValueError("Germination exceeds mature seeds")
    return [],True,{"fruit":fruit,"seeds":seeds,"germinated":germinated}

def run_audit(assignments,execution,followup):
    aa=index(assignments,"assignments")
    ee=index(execution,"execution")
    ff=index(followup,"followup")
    for inp,name in ((ee,"execution"),(ff,"followup")):
        extra=set(inp)-set(aa)
        if extra:
            raise ValueError("Unknown "+name+" IDs not present in assignment list")
    out_rows=[]
    blocks=defaultdict(list)
    for k,plan in sorted(aa.items()):
        if plan.get("not_actual_pollination")!="TRUE":
            raise ValueError("Input must be a PRE-ASSIGNMENT from randomizer")
        if plan.get("allocation_status")!="ASSIGNED_NOT_APPLIED":
            if k in ee or k in ff:
                raise ValueError("Unassigned/nonverified flower appears in execution/followup")
            continue
        morph=plan.get("sex_morph")
        treatment=plan.get("treatment")
        if morph not in MORPHS or treatment not in MORPHS[morph]:
            raise ValueError("Invalid morph/treatment combination")
        er=ee.get(k)
        fr=ff.get(k)
        flags,pp1=execution_quality(er,treatment)
        followflags,completed,values=fruit_quality(fr)
        flags+=followflags
        performed_as_assigned=bool(er and er.get("execution_status")=="performed"
            and er.get("actual_treatment")==treatment)
        if fr and er and er.get("execution_date") and fr.get("assessment_date") and completed:
            if date.fromisoformat(fr["assessment_date"])<date.fromisoformat(er["execution_date"]):
                raise ValueError("Fruit assessment predates experimental pollination")
        per_protocol=bool(pp1 and completed)
        r={
            "population_id":k[0],"genet_id":k[1],"branch_id":k[2],"flower_id":k[3],
            "sex_morph":morph,"assigned_treatment":treatment,
            "execution_status":er.get("execution_status") if er else "unrecorded",
            "performed_as_assigned":int(performed_as_assigned),
            "per_protocol_valid":int(per_protocol),
            "intention_to_treat_observed":int(completed),
            "mature_fruit":values["fruit"] if values is not None else "",
            "mature_seeds":values["seeds"] if values is not None else "",
            "germinated_seeds":values["germinated"] if values is not None and values["germinated"] is not None else "",
            "quality_flags":";".join(flags) if flags else "none",
        }
        out_rows.append(r)
        blocks[(morph,treatment)].append(r)
    report=[]
    for morph,arms in MORPHS.items():
        for arm in arms:
            data=blocks[(morph,arm)]
            pp=[r for r in data if r["per_protocol_valid"]]
            observed=[r for r in data if r["intention_to_treat_observed"]]
            seed_genets={ (r["population_id"],r["genet_id"]) for r in pp
                         if r["mature_seeds"]>0}
            pp_genets={ (r["population_id"],r["genet_id"]) for r in pp}
            gate="EXPLORATORY_ONLY"
            if morph=="female_morph" and arm=="bag_hand_cross":
                gate="POSITIVE_CONTROL_SEED_OBSERVED" if seed_genets else "POSITIVE_CONTROL_NOT_OBSERVED"
            if morph=="male_biased_morph" and arm in ("bag_hand_self","bag_hand_cross"):
                gate="SEED_PRODUCED_NOT_PARENTHOOD_VERIFIED" if seed_genets else "NO_SEED_OBSERVED_NOT_PROOF_OF_STERILITY"
            report.append(dict(
                sex_morph=morph,treatment=arm,planned_flowers=len(data),
                genets_assigned=len({(r["population_id"],r["genet_id"]) for r in data}),
                performed_as_assigned=sum(r["performed_as_assigned"] for r in data),
                protocol_deviations=sum(r["execution_status"]=="deviation" or
                    "different_treatment" in r["quality_flags"] for r in data),
                missing_execution=sum(r["execution_status"]=="unrecorded" for r in data),
                lost_or_pending=sum("outcome_pending_or_lost" in r["quality_flags"] or
                    "followup_unrecorded" in r["quality_flags"] for r in data),
                completed_itt_flowers=len(observed),completed_pp_flowers=len(pp),
                pp_mature_fruit=sum(r["mature_fruit"] for r in pp),
                pp_seed_bearing_flowers=sum(r["mature_seeds"]>0 for r in pp),
                pp_total_mature_seeds=sum(r["mature_seeds"] for r in pp),
                pp_germinated_seeds=sum(r["germinated_seeds"] for r in pp
                                       if r["germinated_seeds"]!=""),
                pp_genets_with_followup=len(pp_genets),
                pp_genets_with_seed=len(seed_genets),
                gate_status=gate,
            ))
    ctrl=next(r for r in report if r["sex_morph"]=="female_morph"
              and r["treatment"]=="bag_hand_cross")
    if ctrl["pp_seed_bearing_flowers"]==0:
        for row in report:
            if row["sex_morph"]=="male_biased_morph" and row["treatment"]=="bag_hand_cross":
                row["gate_status"]="H_CROSS_NOT_INTERPRETABLE_WITHOUT_F_CROSS_POSITIVE"
    return out_rows,report

def write_csv(rows,fieldnames,path):
    with Path(path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def self_test():
    buds=[]
    for morph,genet in (("male_biased_morph","H"),
                       ("female_morph","F")):
        for i in range(10):
            buds.append(dict(population_id="SYNTHETIC",genet_id=genet,
                             branch_id="b"+str(i%2),flower_id=f"f{i:02d}",
                             sex_morph=morph,voucher_checked="yes",
                             flower_stage="predehiscence"))
    aa=allocate(buds,20271114)
    exe,fr=[],[]
    for r in aa:
        arm=r["treatment"]
        er={field:r[field] for field in ID}
        er.update(actual_treatment=arm,execution_status="performed",
                  execution_date="2027-02-20",
                  bag_integrity="intact" if arm in BAG else "not_applicable",
                  pollen_viability_qc="passed" if arm in HAND_CROSS|HAND_SELF else "not_applicable",
                  cross_donor_genet_id="X" if arm in HAND_CROSS else "",
                  self_donor_genet_id=r["genet_id"] if arm in HAND_SELF else "",
                  operator_notes="")
        f={field:r[field] for field in ID}
        fruit=int((r["genet_id"]=="F" and arm=="bag_hand_cross")
                or (r["genet_id"]=="H" and arm=="bag_hand_self"))
        f.update(assessment_date="2027-06-20",followup_status="completed",
                 fruit_mature=str(fruit),mature_seed_count=str(fruit),
                 seed_germinated_count="",observer_notes="SYNTHETIC")
        exe.append(er);fr.append(f)
    x,summary=run_audit(aa,exe,fr)
    assert len(x)==20
    cross_F=next(s for s in summary if s["sex_morph"]=="female_morph" and
                 s["treatment"]=="bag_hand_cross")
    assert cross_F["pp_seed_bearing_flowers"]>=1
    assert cross_F["gate_status"]=="POSITIVE_CONTROL_SEED_OBSERVED"
    self_H=next(s for s in summary if s["sex_morph"]=="male_biased_morph" and
                 s["treatment"]=="bag_hand_self")
    assert self_H["pp_seed_bearing_flowers"]>=1
    # A pollen donor matching the maternal genet invalidates per-protocol outcross.
    donor_test=[r.copy() for r in exe]
    hit=next(r for r in donor_test if r["actual_treatment"]=="bag_hand_cross")
    hit["cross_donor_genet_id"]=hit["genet_id"]
    _,s2=run_audit(aa,donor_test,fr)
    assert sum(row["completed_pp_flowers"] for row in s2)==19
    # Missing outcome means unknown, not 0 fruit.
    _,s3=run_audit(aa,exe,fr[:-1])
    assert sum(row["completed_itt_flowers"] for row in s3)==19
    # Duplicate execution must fail.
    try:
        run_audit(aa,exe+[exe[0]],fr)
    except ValueError:
        pass
    else:
        raise AssertionError("Duplicate flower was accepted")
    # All 20 actual fruits are SYNTHETIC test fixtures only.
    print("PILOT_OUTCOME_QC_TEST=SUCCESS; 20 SYNTHETIC flowers only; "
          "no data from real Onishibari plants",flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--assignment",type=Path)
    parser.add_argument("--execution",type=Path)
    parser.add_argument("--followup",type=Path)
    parser.add_argument("--out-qc",type=Path)
    parser.add_argument("--out-summary",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        self_test();return
    for field in ("assignment","execution","followup","out_qc","out_summary"):
        if getattr(args,field) is None:
            parser.error("All --assignment --execution --followup --out-qc --out-summary required")
    if not all("private" in str(getattr(args,x)).split("/") for x in
               ("assignment","execution","followup","out_qc","out_summary")):
        parser.error("All inputs/outputs must live inside a directory named private")
    aa=load(args.assignment, ID+("sex_morph","treatment","allocation_status","not_actual_pollination"))
    ex=load(args.execution,EXEC_HEADER,allow_empty=True)
    fu=load(args.followup,FRUIT_HEADER,allow_empty=True)
    checked,summary=run_audit(aa,ex,fu)
    args.out_qc.parent.mkdir(parents=True,exist_ok=True)
    args.out_summary.parent.mkdir(parents=True,exist_ok=True)
    write_csv(checked,tuple(checked[0]) if checked else
              ("population_id","genet_id","flower_id","quality_flags"),args.out_qc)
    write_csv(summary,REPORT_HEADER,args.out_summary)
    print("PILOT_QC=AGGREGATED; status=NOT_CLINICAL_OR_CAUSAL_PROOF; "
          "data is user-provided field observations",flush=True)
    print("QC rows",len(checked),"Study arms",len(summary),flush=True)
    print("Outputs private:",args.out_qc,args.out_summary,flush=True)

if __name__=="__main__":
    main()
