#!/usr/bin/env python3
"""v0.47 optional 2027 Onishibari self-pollen-source execution + offspring QC.

NO field data are created or held here. All real input and output files must
remain PRIVATE. --self-test runs synthetic scenarios, including intentional
contamination, missingness and source-fidelity failures. Original v1 pollination
audit remains unchanged and incompatible with the optional v45 treatment arm.

Audits intended vs performed donor identity, dedicated same-genet source,
mature seed, germination and exact parental-route/genotype-supported *germinated*
offspring. Does not infer evolutionary fitness or sterility from zeros.
"""
from __future__ import annotations
import argparse,csv,json,tempfile
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path
from randomize_2027_self_pollen_site_v45 import (
    H_ARMS,F_ARMS,ARMS,VERSION,OUTPUT_COLUMNS,assign,validate
)
from randomize_2027_pilot import COLUMNS

ID=("population_id","genet_id","branch_id","flower_id")
EXEC_HEADER=ID+(
 "actual_treatment","execution_status","execution_date",
 "bag_integrity","pollen_viability_qc","pollen_source_genet_id",
 "pollen_source_branch_id","pollen_source_flower_id",
 "pollen_source_genet_identity_confirmed","pollen_quantity_count",
 "operator_notes",
)
FOLLOWUP_HEADER=ID+(
 "assessment_date","followup_status","fruit_mature",
 "mature_seed_count","germinated_seed_count",
 "genetically_confirmed_expected_parentage_germinated_count",
 "offspring_parentage_qc","observer_notes",
)
VALID_EXEC={"performed","missed","deviation"}
VALID_BAG={"intact","breached","unknown","not_applicable"}
VALID_POLLEN={"passed","failed","not_checked","not_applicable"}
VALID_DONOR={"yes","no","not_checked","not_applicable"}
VALID_FOLLOW={"completed","pending","lost"}
HAND_SELF={"bag_self_same_flower","bag_self_other_flower"}
HAND_CROSS={"bag_hand_cross","open_plus_cross"}
HAND=HAND_SELF|HAND_CROSS
BAG={"bag_no_hand","bag_sham","bag_self_same_flower","bag_self_other_flower","bag_hand_cross"}
DONOR="pollen_donor_same_genet"

def read_csv(path,required,allow_empty=False):
    with Path(path).open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames is None or not set(required).issubset(set(reader.fieldnames)):
            raise ValueError(f"Missing required columns: {path}")
        rows=[{k:(v.strip() if v is not None else "") for k,v in r.items()}
              for r in reader]
    if not rows and not allow_empty:raise ValueError(f"Empty file: {path}")
    return rows

def key(r):
    tup=tuple(r.get(k,"") for k in ID)
    if any(not x for x in tup):raise ValueError("Incomplete four-field floral ID")
    return tup

def index(rows,name):
    result={}
    for r in rows:
        k=key(r)
        if k in result:raise ValueError(f"Duplicate {name} flower id")
        result[k]=r
    return result

def nonnegative(value,name):
    if value is None or value=="" or value=="NA":
        raise ValueError(f"Missing required numeric outcome: {name}")
    try: n=int(value)
    except (TypeError,ValueError):raise ValueError(f"Invalid integer count {name}")
    if n<0:raise ValueError(f"Negative {name}")
    return n

def dateval(value,name):
    if not value:raise ValueError("Missing required date: "+name)
    return date.fromisoformat(value)

def schedule_guard(schedule):
    group=defaultdict(list)
    for s in schedule:
        if s.get("protocol_version")!=VERSION or s.get("not_actual_pollination")!="TRUE":
            raise ValueError("This QC only accepts the original v45 unperformed schedule")
        group[(s["population_id"],s["genet_id"])].append(s)
    for genet,ss in group.items():
        ready=[r for r in ss if r["flower_stage"]=="predehiscence"
               and r["voucher_checked"]=="yes"]
        chosen=[r for r in ss if r["allocation_status"]=="ASSIGNED_NOT_APPLIED"]
        if not chosen:continue
        morph=chosen[0]["sex_morph"]
        if morph not in ARMS or any(r["sex_morph"]!=morph for r in ss):
            raise ValueError("Conflicting genet morph labels")
        if len(chosen)!=len(ARMS[morph]) or set(r["planned_treatment"] for r in chosen)!=set(ARMS[morph]):
            raise ValueError("Arm set is not complete for genet")
        if len(chosen)>10 or len(chosen)>len(ready)//2:raise ValueError("Unsafe genotypic bud allocation")
        if morph=="male_biased_morph":
            src=next(r for r in chosen if r["planned_treatment"]==DONOR)
            recv=next(r for r in chosen if r["planned_treatment"]=="bag_self_other_flower")
            same=next(r for r in chosen if r["planned_treatment"]=="bag_self_same_flower")
            if src["mature_seed_followup_expected"]!="NO" or src["role"]!="pollen_source_only":
                raise ValueError("Donor cannot have a maternal endpoint")
            if (recv["same_genet_donor_branch_id"],recv["same_genet_donor_flower_id"])!=(src["branch_id"],src["flower_id"]):
                raise ValueError("Same-genet-other-flower source mismatched to dedicated donor")
            if (same["same_genet_donor_branch_id"],same["same_genet_donor_flower_id"])!=(same["branch_id"],same["flower_id"]):
                raise ValueError("Within-flower pollen source not self")
        for r in chosen:
            if r["planned_treatment"]!=DONOR and r["mature_seed_followup_expected"]!="YES":
                raise ValueError("Assigned recipient has no endpoint")
    return group

def exec_flags(plan,event,exec_index):
    if event is None:return ["execution_missing"],False,None
    planned=plan["planned_treatment"]
    status=event.get("execution_status")
    if status not in VALID_EXEC:raise ValueError("Unrecognized execution status")
    if event.get("bag_integrity") not in VALID_BAG or event.get("pollen_viability_qc") not in VALID_POLLEN:
        raise ValueError("Unrecognized bag or pollen quality code")
    if event.get("pollen_source_genet_identity_confirmed") not in VALID_DONOR:
        raise ValueError("Unrecognized source genet verification")
    if status=="missed":
        if event.get("actual_treatment") not in ("","NA"):
            raise ValueError("Missed pollination cannot record actual treatment")
        return ["unperformed"],False,None
    d=dateval(event.get("execution_date"),"execution_date")
    flags=[]
    if status=="deviation" or event.get("actual_treatment")!=planned:
        flags.append("different_arm_or_deviation")
    if planned in BAG and event["bag_integrity"]!="intact":
        flags.append("bag_integrity_unverified")
    if planned in HAND and event["pollen_viability_qc"]!="passed":
        flags.append("pollen_viability_unverified")
    if planned in HAND:
        if event.get("pollen_source_genet_identity_confirmed")!="yes":
            flags.append("donor_genet_identity_unverified")
        try:
            if nonnegative(event.get("pollen_quantity_count"),"pollen_quantity_count")<=0:
                flags.append("nonpositive_pollen_quantity")
        except ValueError:flags.append("missing_or_invalid_pollen_quantity")
        srcgen=event.get("pollen_source_genet_id","")
        srcbranch=event.get("pollen_source_branch_id","")
        srcflower=event.get("pollen_source_flower_id","")
        if planned in HAND_SELF and srcgen!=plan["genet_id"]:
            flags.append("self_not_same_genet")
        if planned in HAND_CROSS and (not srcgen or srcgen==plan["genet_id"]):
            flags.append("outcross_not_different_genet")
        if not srcbranch or not srcflower:
            flags.append("pollen_source_flower_identity_missing")
        if planned=="bag_self_same_flower" and (srcbranch,srcflower)!=(plan["branch_id"],plan["flower_id"]):
            flags.append("within_flower_source_not_recipient")
        if planned=="bag_self_other_flower":
            expect=(plan["same_genet_donor_branch_id"],plan["same_genet_donor_flower_id"])
            if (srcbranch,srcflower)!=expect:
                flags.append("different_flower_source_not_preadjusted_dedicated_donor")
            else:
                donor_key=(plan["population_id"],plan["genet_id"],srcbranch,srcflower)
                donor_event=exec_index.get(donor_key)
                if donor_event is None or donor_event.get("execution_status")!="performed" or donor_event.get("actual_treatment")!=DONOR:
                    flags.append("dedicated_pollen_donor_execution_unverified")
                elif dateval(donor_event.get("execution_date"),"dedicated donor execution")>d:
                    flags.append("dedicated_donor_used_after_recipient_pollination")
    else:
        if event.get("pollen_source_genet_id") not in ("","NA") and planned!=DONOR:
            flags.append("pollen_source_populated_on_no_hand_arm")
        if planned==DONOR:
            # This is NOT a maternal recipient. Source must have passed pollen QC.
            if event.get("pollen_viability_qc")!="passed":
                flags.append("dedicated_donor_pollen_viability_unverified")
    return flags,len(flags)==0,d

def outcome(row,planned):
    if row is None:return ["followup_missing"],False,None
    state=row.get("followup_status")
    if state not in VALID_FOLLOW:raise ValueError("Invalid followup status")
    if state!="completed":
        if any(row.get(k) not in ("","NA") for k in (
            "fruit_mature","mature_seed_count","germinated_seed_count",
            "genetically_confirmed_expected_parentage_germinated_count")):
            raise ValueError("Incomplete followup cannot be coded as observed seed zero")
        return ["outcome_pending_or_lost"],False,None
    assessment=dateval(row.get("assessment_date"),"assessment_date")
    if row.get("fruit_mature") not in ("0","1"):
        raise ValueError("Mature fruit must explicitly be 0/1 when complete")
    fruit=int(row["fruit_mature"])
    mature=nonnegative(row.get("mature_seed_count"),"mature_seed_count")
    if fruit==0 and mature>0:raise ValueError("Mature seeds with fruit=0")
    germ_raw=row.get("germinated_seed_count")
    germ=None if germ_raw in ("","NA") else nonnegative(germ_raw,"germinated_seed_count")
    if germ is not None and germ>mature:raise ValueError("Germination > mature seeds")
    verified_raw=row.get("genetically_confirmed_expected_parentage_germinated_count")
    verified=None if verified_raw in ("","NA") else nonnegative(verified_raw,"verified parental viable count")
    if verified is not None:
        if germ is None:raise ValueError("Genotyped viable count requires measured germination")
        if verified>germ:raise ValueError("Genotyped germinated seed count > germinated seeds")
        if row.get("offspring_parentage_qc")!="passed":
            raise ValueError("Genetically verified seeds require passed offspring parentage QC")
    if row.get("offspring_parentage_qc") not in ("passed","failed","pending","not_done"):
        raise ValueError("Invalid offspring parentage QC")
    return [],True,{"fruit":fruit,"mature":mature,"germinated":germ,
                     "strict":verified if verified is not None else None,"date":assessment}

def audit(schedule,events,followups):
    aa=index(schedule,"allocation")
    ee=index(events,"execution")
    ff=index(followups,"followup")
    schedule_guard(schedule)
    assigned={k:v for k,v in aa.items() if v["allocation_status"]=="ASSIGNED_NOT_APPLIED"}
    if not set(ee)<=set(assigned) or not set(ff)<=set(assigned):
        raise ValueError("Execution or outcome contains unassigned source flower")
    rows=[]
    by_arm=defaultdict(list)
    for k,s in sorted(assigned.items()):
        plan=s["planned_treatment"]
        e=ee.get(k)
        r=ff.get(k)
        flags,ok,edate=exec_flags(s,e,ee)
        if plan==DONOR and r is not None:raise ValueError("Pollen-only donor wrongly has maternal fruit followup")
        if plan==DONOR:
            complete=False
            seed=None
        else:
            more,complete,seed=outcome(r,plan)
            flags+=more
            if seed and edate is not None and seed["date"]<edate:
                raise ValueError("Fruit/seed followup predates source treatment")
            if seed and (not ok):
                flags.append("seed_result_cannot_establish_assigned_pollen_source")
        good=bool(plan!=DONOR and ok and complete and seed is not None)
        strict=seed["strict"] if good and seed is not None else None
        row={"sex_morph":s["sex_morph"],"planned_arm":plan,
             "genet_private":s["genet_id"],"passed_execution":int(ok),
             "completed":int(complete),"per_protocol_completed":int(good),
             "fruit":seed["fruit"] if seed else None,
             "mature":seed["mature"] if seed else None,
             "germinated":seed["germinated"] if seed else None,
             "strict_genetically_confirmed_germinated":strict,
             "flags":flags}
        rows.append(row)
        by_arm[(s["sex_morph"],plan)].append(row)
    report=[]
    for key,group in sorted(by_arm.items()):
        strict=[r for r in group if r["strict_genetically_confirmed_germinated"] is not None]
        perprotocol=[r for r in group if r["per_protocol_completed"]]
        report.append({"sex_morph":key[0],"arm":key[1],"assigned_genets":len(group),
          "source_only_donor":key[1]==DONOR,
          "executed_as_planned":sum(r["passed_execution"] for r in group),
          "completed_followup":sum(r["completed"] for r in group),
          "per_protocol_complete_flowers":len(perprotocol),
          "reported_mature_seeds_per_protocol":sum(r["mature"] for r in perprotocol),
          "with_genetically_verified_viable_outcome":len(strict),
          "positive_genets_with_genetically_confirmed_expected_parentage_germinated_seed":
            len({r["genet_private"] for r in strict if r["strict_genetically_confirmed_germinated"]>0}),
          "confirmed_expected_parentage_germinated_seed_total":
            sum(r["strict_genetically_confirmed_germinated"] for r in strict),
          "missing_or_incomplete_count":len(group)-len(perprotocol),
          "notes":"GENET UNIT; raw flower rows and missingness must not be treated as independent successes"})
    return {"status":"AUDIT_ONLY_NO_CAUSAL_INTERPRETATION",
      "source_genets":len({(x["population_id"],x["genet_id"]) for x in assigned.values()}),
      "assigned_role_count":len(assigned),"arm_groups":report,
      "n_flagged_assigned_roles":sum(bool(x["flags"]) for x in rows),
      "n_unknown_outcome_roles":sum(x["completed"]==0 and x["planned_arm"]!=DONOR for x in rows),
      "no_observational_zero_filled":True,
      "stratify_by_genet_and_pollen_source_not_flower_replicate":True}

def save_summary(summary,path):
    out=Path(path)
    if "data/private/" not in str(out).replace("\\","/"):
        raise ValueError("Audit output must be in data/private/, never public repo")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")

def selftest():
    # All observations here are SYNTHETIC fixtures only; not field evidence.
    buds=[]
    for sex,g,n in (("male_biased_morph","HFAKE",14),("female_morph","FFAKE",10)):
        for i in range(n):
            buds.append(dict(population_id="SYNTHETIC",genet_id=g,
             branch_id=f"b{i%2}",flower_id=f"f{i:02d}",
             sex_morph=sex,voucher_checked="yes",flower_stage="predehiscence"))
    plans=assign(buds,20271117)
    assert validate(buds,plans)["assigned_total_plant_bud_roles"]==12
    active=[s for s in plans if s["allocation_status"]=="ASSIGNED_NOT_APPLIED"]
    donor=next(s for s in active if s["planned_treatment"]==DONOR)
    events=[]
    outcomes=[]
    for p in active:
        arm=p["planned_treatment"]
        er={k:p[k] for k in ID}
        er.update(actual_treatment=arm,execution_status="performed",
          execution_date="2027-02-23",bag_integrity="intact" if arm in BAG else "not_applicable",
          pollen_viability_qc="passed" if arm in HAND or arm==DONOR else "not_applicable",
          pollen_source_genet_id=(p["genet_id"] if arm in HAND_SELF else
             "OTHER_GENET" if arm in HAND_CROSS else ""),
          pollen_source_branch_id=(p["same_genet_donor_branch_id"] if arm in HAND_SELF else
             "crossbud" if arm in HAND_CROSS else ""),
          pollen_source_flower_id=(p["same_genet_donor_flower_id"] if arm in HAND_SELF else
             "otherflower" if arm in HAND_CROSS else ""),
          pollen_source_genet_identity_confirmed="yes" if arm in HAND else "not_applicable",
          pollen_quantity_count="20" if arm in HAND else "",operator_notes="SYNTHETIC")
        events.append(er)
        if arm==DONOR:continue
        rr={k:p[k] for k in ID}
        vi=int(arm=="bag_hand_cross" and p["genet_id"]=="HFAKE")
        rr.update(assessment_date="2027-06-10",followup_status="completed",
          fruit_mature=str(vi),mature_seed_count=str(vi),
          germinated_seed_count=str(vi),
          genetically_confirmed_expected_parentage_germinated_count=str(vi),
          offspring_parentage_qc="passed",observer_notes="SYNTHETIC")
        outcomes.append(rr)
    result=audit(plans,events,outcomes)
    assert result["assigned_role_count"]==12
    assert result["n_unknown_outcome_roles"]==0
    assert result["n_flagged_assigned_roles"]==0
    hcross=next(x for x in result["arm_groups"] if x["sex_morph"]=="male_biased_morph" and x["arm"]=="bag_hand_cross")
    assert hcross["positive_genets_with_genetically_confirmed_expected_parentage_germinated_seed"]==1
    donorgrp=next(x for x in result["arm_groups"] if x["arm"]==DONOR)
    assert donorgrp["source_only_donor"] and donorgrp["completed_followup"]==0
    # Source identity tampering must drop PP outcome, never count success.
    wrong=[dict(x) for x in events]
    geit=next(x for x in wrong if x["actual_treatment"]=="bag_self_other_flower")
    geit["pollen_source_flower_id"]="INCORRECT"
    bad=audit(plans,wrong,outcomes)
    row=next(x for x in bad["arm_groups"] if x["arm"]=="bag_self_other_flower")
    assert row["per_protocol_complete_flowers"]==0
    # Absent followups are not zero observed seeds.
    missing=audit(plans,events,outcomes[:-1])
    assert missing["n_unknown_outcome_roles"]==1
    # Incorrect positive seed evidence cannot exceed observed germination.
    impossible=[dict(x) for x in outcomes]
    impossible[0]["genetically_confirmed_expected_parentage_germinated_count"]="2"
    try:audit(plans,events,impossible)
    except ValueError:pass
    else:raise AssertionError("Impossible viable offspring count was accepted")
    # Dedicated source flower cannot be reused as a recipient outcome.
    stolen=[dict(x) for x in outcomes]+[{**{k:donor[k] for k in ID},
        "assessment_date":"2027-06-10","followup_status":"completed",
        "fruit_mature":"0","mature_seed_count":"0","germinated_seed_count":"0",
        "genetically_confirmed_expected_parentage_germinated_count":"0",
        "offspring_parentage_qc":"passed","observer_notes":"SYNTHETIC"}]
    try:audit(plans,events,stolen)
    except ValueError:pass
    else:raise AssertionError("Dedicated donor used as maternal fruit endpoint")
    return {"status":"PASS_V47_SYNTHETIC_SOURCE_AND_OFFSPRING_QC",
       "fake_total_buds":len(buds),"fake_genets":2,
       "fake_assigned_roles":12,"fake_maternal_recipients":11,
       "test_source_mismatch_rejected_from_per_protocol":True,
       "test_missing_followup_not_zero_filled":True,
       "test_more_verified_seeds_than_germinated_hard_fails":True,
       "test_dedicated_pollen_donor_cannot_be_maternal_receiver":True,
       "not_observed_2027_focal_field_data":True}

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--schedule",type=Path)
    a.add_argument("--execution",type=Path)
    a.add_argument("--followup",type=Path)
    a.add_argument("--summary-private",type=Path)
    a.add_argument("--self-test",action="store_true")
    p=a.parse_args()
    if p.self_test:
        print("ONISHIBARI_SOURCE_2027_QC_V47",json.dumps(selftest()),flush=True)
        return
    if not all((p.schedule,p.execution,p.followup,p.summary_private)):
        a.error("Supply private v45 schedule, execution, followup and summary paths or --self-test")
    summary=audit(read_csv(p.schedule,OUTPUT_COLUMNS),
                  read_csv(p.execution,EXEC_HEADER,allow_empty=True),
                  read_csv(p.followup,FOLLOWUP_HEADER,allow_empty=True))
    save_summary(summary,p.summary_private)
    print("ONISHIBARI_SOURCE_2027_QC_V47",json.dumps({
        "status":"audit_success_summary_private",
        "n_assigned_roles":summary["assigned_role_count"],
        "n_missing_recipient_followup":summary["n_unknown_outcome_roles"],
        "n_flagged_role":summary["n_flagged_assigned_roles"],
        "never_treat_missing_as_zero":True},ensure_ascii=False))

if __name__=="__main__":main()
