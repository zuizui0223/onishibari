#!/usr/bin/env python3
"""v0.49: PRIVATE seed-to-reproductive-adult longitudinal stage audit.

Future Onishibari D. pseudomezereum v45 optional source-pollen pilot.
Joins planned recipient flower, actual pollen source execution (v47),
completed mother-flower fruit/germination/parentage counts (v47), and
repeated fate records of IDENTIFIABLE offspring (new cohort ledger).

Does NOT create observed data, infer absolute lifetime fitness from a flower,
or turn missing seedlings into deaths. A named seed and the mother plant are
both nested; only genets are independent field sampling units. The strict
endpoint is cross/self lineage validated + germinated/recruited/reproduced;
flowering adult is not independently confirmed mature offspring.

All REAL inputs and generated summaries must stay under data/private/.
Only self-test with synthetic fake plants is allowed in public CI.
"""
from __future__ import annotations
import argparse,csv,json,tempfile
from collections import defaultdict
from datetime import date
from pathlib import Path
from audit_2027_self_pollen_site_v47 import (
  ID,EXEC_HEADER,FOLLOWUP_HEADER,read_csv,index,schedule_guard,
  exec_flags,outcome,HAND_SELF,HAND_CROSS,DONOR
)
from randomize_2027_self_pollen_site_v45 import assign,validate,OUTPUT_COLUMNS,VERSION

LEDGER_HEADER=ID+(
 "seed_id","observation_date","stage","fate_status",
 "mating_origin","mating_origin_genotype_qc",
 "reproductive_evidence","growth_context","observer_notes",
)
STAGES={"mature_seed":0,"germinated":1,"established_seedling":2,
        "juvenile":3,"reproductive_adult":4}
FATES={"alive","dead","lost","pending"}
ORIGINS={"self","outcross","unknown"}
ORIGIN_QC={"passed","failed","pending","not_done"}
GROWTH={"natural_field","common_garden","growth_chamber","transplanted","unknown"}
REPRODUCTIVE={"none","flowering_observed","offspring_parentage_confirmed"}

def seedkey(r):
    if not r.get("seed_id"):raise ValueError("Missing seed identity")
    return (*tuple(r[k] for k in ID),r["seed_id"])

def arm_compatible(arm,origin):
    if arm in HAND_SELF:return origin=="self"
    if arm in HAND_CROSS:return origin=="outcross"
    if arm in ("bag_no_hand","bag_sham","open"):
        return origin in ("self","outcross")
    return origin in ORIGINS

def time_and_stage(rows):
    byseed=defaultdict(list)
    per_flower=defaultdict(set)
    for r in rows:
        sk=seedkey(r)
        if r["stage"] not in STAGES:raise ValueError("Unknown offspring development stage")
        if r["fate_status"] not in FATES:raise ValueError("Unknown offspring fate")
        if r["mating_origin"] not in ORIGINS or r["mating_origin_genotype_qc"] not in ORIGIN_QC:
            raise ValueError("Unknown genetic origin/QC")
        if r["growth_context"] not in GROWTH or r["reproductive_evidence"] not in REPRODUCTIVE:
            raise ValueError("Unknown growth context or adult reproduction evidence")
        dt=date.fromisoformat(r["observation_date"])
        if r["fate_status"]!="alive" and r["reproductive_evidence"]!="none":
            raise ValueError("Only observed alive seedling/adult can provide reproductive evidence")
        if r["reproductive_evidence"]!="none" and r["stage"]!="reproductive_adult":
            raise ValueError("Reproductive claim requires observed reproductive adult stage")
        if r["stage"]=="reproductive_adult" and r["fate_status"]=="alive" and r["reproductive_evidence"]=="none":
            raise ValueError("Reproductive adult cannot be a size-only classification")
        byseed[sk].append((dt,r))
        per_flower[sk[:4]].add(sk)
    for sk,history in byseed.items():
        history.sort(key=lambda x:x[0])
        if len({x[0] for x in history})!=len(history):
            raise ValueError("Duplicate source seed/census date")
        last_stage=-1;was_dead=False
        origin_set=set();contexts=set()
        for dt,r in history:
            f=r["fate_status"]
            stage=STAGES[r["stage"]]
            if was_dead and f=="alive":raise ValueError("Offspring revived after confirmed death")
            if f=="alive":
                if stage<last_stage:raise ValueError("Development stage regressed in a living offspring")
                last_stage=stage
            if f=="dead":was_dead=True
            if r["mating_origin_genotype_qc"]=="passed":
                if r["mating_origin"]=="unknown":
                    raise ValueError("Passed genetic origin QC cannot have unknown mating route")
                origin_set.add(r["mating_origin"])
            contexts.add(r["growth_context"])
        if len(origin_set)>1:raise ValueError("Contradictory verified self versus outcross provenance for same seed")
    return byseed,per_flower

def audit(schedule,executions,followups,ledger):
    assignments=index(schedule,"v45 assigned")
    execs=index(executions,"v47 executed")
    follows=index(followups,"v47 completed fruit")
    schedule_guard(schedule)
    assigned={k:v for k,v in assignments.items()
              if v["allocation_status"]=="ASSIGNED_NOT_APPLIED"}
    for table in (execs,follows):
        if not set(table)<=set(assigned):
            raise ValueError("Execution or mother fruit outcomes refer to unallocated flower")
    byseed,byflower=time_and_stage(ledger)
    for flower in byflower:
        if flower not in assigned:raise ValueError("Seed source flower was unallocated")
        if assigned[flower]["planned_treatment"]==DONOR:
            raise ValueError("Dedicated pollen-only flower falsely produced maternal offspring records")
    group=defaultdict(list)
    failed_parentage=[]
    all_seed_count=0
    n_ever_mature=n_ever_germ=n_ever_established=n_ever_adult=n_ever_confirmed_offspring=0
    n_lost_or_pending_seed_histories=0
    n_strict_germ,n_strict_established,n_strict_adult,n_strict_repro=0,0,0,0
    censored_lost=0
    for flower,seedset in byflower.items():
        plan=assigned[flower]
        arm=plan["planned_treatment"]
        ef=execs.get(flower)
        exec_issues,execution_ok,executed_on=exec_flags(plan,ef,execs)
        fr=follows.get(flower)
        follow_issues,completed,counts=outcome(fr,arm)
        # Preserve source/aggregate QC defects as noncausal, nonzero flags.
        failed_parentage.extend("execution_"+str(flag) for flag in exec_issues)
        failed_parentage.extend("followup_"+str(flag) for flag in follow_issues)
        if completed and executed_on is not None and counts["date"]<executed_on:
            raise ValueError("Mother seed/germination followup predates actual pollination")
        # An individual offspring ledger is NOT assumed complete. Its known
        # germinated and mature seed counts cannot exceed the mother's
        # verified completed aggregate. All unobserved seeds stay unknown.
        seeds=list(seedset)
        earliest_germ=0;verified_strict_germ=0;all_mature=0
        for sk in seeds:
            h=byseed[sk]
            if executed_on is not None and any(dt<executed_on and r["fate_status"]=="alive" for dt,r in h):
                raise ValueError("Observed seed life stage predates actual pollination")
            alive=[r for _,r in h if r["fate_status"]=="alive"]
            if not alive:continue
            all_mature+=1
            seen_germ=any(STAGES[r["stage"]]>=1 for r in alive)
            if seen_germ:earliest_germ+=1
            verified_route=next((r["mating_origin"] for _,r in reversed(h)
              if r["mating_origin_genotype_qc"]=="passed"),None)
            if (seen_germ and verified_route is not None and
                arm_compatible(arm,verified_route)):
                verified_strict_germ+=1
        if completed:
            # Only when all seeds and their exact germination events are
            # accounted for, completed aggregate cannot precede germination.
            direct_germ_dates=[]
            for sk in seeds:
                stamps=[dt for dt,r in byseed[sk]
                        if r["stage"]=="germinated" and r["fate_status"]=="alive"]
                if stamps:direct_germ_dates.append(min(stamps))
            if (len(seeds)==counts["mature"] and earliest_germ>0 and
                counts["germinated"] is not None and
                earliest_germ==counts["germinated"] and
                len(direct_germ_dates)==earliest_germ and
                max(direct_germ_dates)>counts["date"]):
                raise ValueError("Mother aggregate counts seed germination before explicit seed event")
            if all_mature>counts["mature"]:
                raise ValueError("Individual seed identities exceed mother's recorded mature seed count")
            if counts["germinated"] is not None and earliest_germ>counts["germinated"]:
                raise ValueError("Offspring histories exceed mother followup germination count")
            if counts["strict"] is not None and verified_strict_germ>counts["strict"]:
                raise ValueError("Matching verified germinated seeds exceed mother's strict molecular count")
        elif earliest_germ or all_mature:
            # Missing/late aggregate followup is a QC gap; do not silently
            # classify it as successful per protocol until corrected.
            failed_parentage.append("seed_history_present_without_completed_mother_aggregate")
        for sk in seeds:
            h=byseed[sk]
            n=1;all_seed_count+=1
            alive=[r for _,r in h if r["fate_status"]=="alive"]
            stage_high=max([STAGES[r["stage"]] for r in alive],default=-1)
            if stage_high>=0:n_ever_mature+=1
            if stage_high>=1:n_ever_germ+=1
            if stage_high>=2:n_ever_established+=1
            if stage_high>=4:n_ever_adult+=1
            latest=h[-1][1]
            if latest["fate_status"] in ("lost","pending"):
                n_lost_or_pending_seed_histories+=1
                if latest["fate_status"]=="lost":censored_lost+=1
            verified_route=next((r["mating_origin"] for _,r in reversed(h)
              if r["mating_origin_genotype_qc"]=="passed"),None)
            if verified_route is not None and not arm_compatible(arm,verified_route):
                failed_parentage.append("mating_origin_conflicts_with_experimental_arm")
            strict=bool(execution_ok and completed and counts is not None and
              counts["strict"] is not None and verified_route is not None and
              arm_compatible(arm,verified_route))
            if strict and stage_high>=1:n_strict_germ+=1
            if strict and stage_high>=2:n_strict_established+=1
            if strict and stage_high>=4:n_strict_adult+=1
            if strict and stage_high>=4 and any(
                  r["reproductive_evidence"]=="offspring_parentage_confirmed" and
                  r["fate_status"]=="alive" for _,r in h):
                n_strict_repro+=1
            if stage_high>=4 and any(
                  r["reproductive_evidence"]=="offspring_parentage_confirmed" and
                  r["fate_status"]=="alive" for _,r in h):
                n_ever_confirmed_offspring+=1
            if strict and latest["fate_status"] in ("lost","pending"):
                # Never mark censored outcome as dead or as a completed
                # lifetime survival fraction.
                pass
            group[(plan["sex_morph"],arm)].append({
              "stage_high":stage_high,"strict":strict,
              "censored":latest["fate_status"] in ("lost","pending")})
    # The stage counts are reported globally without source ids, and never
    # divided by the total assignments as if followup were complete.
    summary={
      "status":"V49_LONGITUDINAL_SEED_STAGE_QC_NOT_NEW_FIELD_RESULTS",
      "source_unique_identified_seed_records":all_seed_count,
      "ever_mature_identified":n_ever_mature,
      "ever_germinated_identified":n_ever_germ,
      "ever_established_seedlings_identified":n_ever_established,
      "ever_reproductive_adults_identified":n_ever_adult,
      "ever_confirmed_own_offspring_reproduced":n_ever_confirmed_offspring,
      "source_route_checked_and_per_protocol_strict_counts":{
          "germinated":n_strict_germ,
          "established_seedling":n_strict_established,
          "reproductive_adult":n_strict_adult,
          "offspring_parentage_confirmed_reproduction":n_strict_repro},
      "missing_or_censored_seed_histories":n_lost_or_pending_seed_histories,
      "lost_or_unlocated_not_assumed_dead":censored_lost,
      "known_source_or_aggregate_problems_count":len(failed_parentage),
      "group_summaries":[{
         "morph":k[0],"planned_treatment":k[1],
         "n_identified_seeds":len(v),
         "n_strict_source_matched":sum(r["strict"] for r in v),
         "n_already_reproductive_adults_strict":sum(r["strict"] and r["stage_high"]>=4 for r in v),
         "n_censored_seed_histories":sum(r["censored"] for r in v)}
         for k,v in sorted(group.items())],
      "not_a_recruitment_survival_rate_without_complete_censoring_and_denominators":True,
      "not_a_lifetime_fitness_measurement":True,
      "no_source_private_ids_or_exact_coordinates_in_summary":True
    }
    return summary

def save_private(out,path):
    path=Path(path)
    if "data/private/" not in str(path).replace("\\","/"):
        raise ValueError("Longitudinal raw or aggregate outputs must be saved under data/private/")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")

def synthetic_test():
    bud=[]
    for morph,g,n in (("male_biased_morph","H_SYNTH",14),("female_morph","F_SYNTH",10)):
        for i in range(n):
            bud.append({"population_id":"SYNTHETIC","genet_id":g,
               "branch_id":f"b{i%2}","flower_id":f"f{i:02d}",
               "sex_morph":morph,"voucher_checked":"yes",
               "flower_stage":"predehiscence"})
    plans=assign(bud,20271117)
    validate(bud,plans)
    active=[a for a in plans if a["allocation_status"]=="ASSIGNED_NOT_APPLIED"]
    donor=next(a for a in active if a["planned_treatment"]==DONOR)
    execs=[];outcomes=[]
    for a in active:
        arm=a["planned_treatment"]
        e={k:a[k] for k in ID}
        gen=a["genet_id"]
        e.update(actual_treatment=arm,execution_status="performed",
          execution_date="2027-03-01",
          bag_integrity="intact",
          pollen_viability_qc="passed" if arm in HAND_SELF|HAND_CROSS|{DONOR} else "not_applicable",
          pollen_source_genet_id=(gen if arm in HAND_SELF else
             "OTHER" if arm in HAND_CROSS else ""),
          pollen_source_branch_id=(
             a["same_genet_donor_branch_id"] if arm in HAND_SELF else
             "OTHER_BRANCH" if arm in HAND_CROSS else ""),
          pollen_source_flower_id=(
             a["same_genet_donor_flower_id"] if arm in HAND_SELF else
             "OTHER_FLOWER" if arm in HAND_CROSS else ""),
          pollen_source_genet_identity_confirmed=(
             "yes" if arm in HAND_SELF|HAND_CROSS else "not_applicable"),
          pollen_quantity_count="20" if arm in HAND_SELF|HAND_CROSS else "",
          operator_notes="SYNTHETIC")
        execs.append(e)
        if arm==DONOR:continue
        f={k:a[k] for k in ID}
        f.update(assessment_date="2028-04-15",followup_status="completed",
          fruit_mature="1" if arm=="bag_hand_cross" else "0",
          mature_seed_count="1" if arm=="bag_hand_cross" else "0",
          germinated_seed_count="1" if arm=="bag_hand_cross" else "0",
          genetically_confirmed_expected_parentage_germinated_count="1" if arm=="bag_hand_cross" else "0",
          offspring_parentage_qc="passed",
          observer_notes="SYNTHETIC")
        outcomes.append(f)
    h=next(a for a in active if a["planned_treatment"]=="bag_hand_cross" and a["genet_id"]=="H_SYNTH")
    k={n:h[n] for n in ID}
    def row(date_,stage,fate,rc="none",origin="outcross",qc="passed"):
        return {**k,"seed_id":"SEED_FAKE_1","observation_date":date_,
           "stage":stage,"fate_status":fate,"mating_origin":origin,
           "mating_origin_genotype_qc":qc,"reproductive_evidence":rc,
           "growth_context":"natural_field","observer_notes":"SYNTHETIC"}
    ledger=[
      row("2027-07-15","mature_seed","alive"),
      row("2028-03-30","germinated","alive"),
      row("2028-07-01","established_seedling","alive"),
      row("2029-05-01","juvenile","alive"),
      row("2032-04-20","reproductive_adult","alive",rc="flowering_observed"),
      row("2033-04-20","reproductive_adult","alive",rc="offspring_parentage_confirmed")
    ]
    # A completed maternal count cannot claim the exact known child's
    # germination occurred before the documented first germination date.
    too_early=[dict(x) for x in outcomes]
    for x in too_early:
        if x["genet_id"]=="H_SYNTH" and x["fruit_mature"]=="1":
            x["assessment_date"]="2027-07-15"
    try:audit(plans,execs,too_early,ledger)
    except ValueError:pass
    else:raise AssertionError("Future germination incorrectly reported as an early observed fact")
    a=audit(plans,execs,outcomes,ledger)
    assert a["source_route_checked_and_per_protocol_strict_counts"]["germinated"]==1
    assert a["source_route_checked_and_per_protocol_strict_counts"]["reproductive_adult"]==1
    assert a["source_route_checked_and_per_protocol_strict_counts"]["offspring_parentage_confirmed_reproduction"]==1
    assert a["known_source_or_aggregate_problems_count"]==0
    # Distinct pollen origin in genetically proven H hand-cross must never be
    # labelled strictly confirmed if source genotype indicates self.
    bad=[dict(x) for x in ledger]
    for x in bad:x["mating_origin"]="self"
    b=audit(plans,execs,outcomes,bad)
    assert b["source_route_checked_and_per_protocol_strict_counts"]["germinated"]==0
    assert b["known_source_or_aggregate_problems_count"]>=1
    # A disappeared recruit must be censored, not automatically coded dead.
    missing=ledger[:3]+[row("2029-05-01","established_seedling","lost")]
    m=audit(plans,execs,outcomes,missing)
    assert m["missing_or_censored_seed_histories"]==1
    assert m["ever_reproductive_adults_identified"]==0
    # Confirmed death followed by an alive juvenile must throw.
    dead=ledger[:3]+[row("2028-09-01","established_seedling","dead")]+ledger[3:]
    try:audit(plans,execs,outcomes,dead)
    except ValueError:pass
    else:raise AssertionError("Resurrection after dead accepted")
    # Donor-only reproductive record must throw.
    donorrow={**row("2027-07-15","mature_seed","alive")}
    for name in ID:donorrow[name]=donor[name]
    try:audit(plans,execs,outcomes,[donorrow])
    except ValueError:pass
    else:raise AssertionError("Donor credited as seed parent")
    # Seed count beyond original mother fruit followup must throw.
    twin=[dict(x) for x in ledger]
    for x in twin:x["seed_id"]="SEED_FAKE_2"
    try:audit(plans,execs,outcomes,ledger+twin)
    except ValueError:pass
    else:raise AssertionError("More seeds than source fruit record accepted")
    return {"status":"PASS_V49_SYNTHETIC_ONISHIBARI_SEED_TO_ADULT_QC",
      "fake_mother_genets":2,"fake_assigned_roles":12,
      "fake_offspring_seed_count":1,"fake_time_points":len(ledger),
      "verified_outcross_reproductive_adult_synthetic":1,
      "origin_mismatch_not_strict":True,
      "missing_offspring_not_dead":True,
      "dead_offspring_does_not_revive":True,
      "dedicated_pollen_source_never_maternal":True,
      "no_more_seed_identities_than_mother_fruit_count":True,
      "source_maternal_followup_must_postdate_explicit_seed_germination":True,
      "no_actual_2027_experiment_observed":True}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--schedule",type=Path)
    p.add_argument("--execution",type=Path)
    p.add_argument("--mother-followup",type=Path)
    p.add_argument("--seed-ledger",type=Path)
    p.add_argument("--output-private",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        print("ONISHIBARI_OFFSPRING_LIFE_STAGE_V49",json.dumps(synthetic_test()),flush=True)
        return
    if not all((a.schedule,a.execution,a.mother_followup,a.seed_ledger,a.output_private)):
        p.error("Use --self-test or supply 4 PRIVATE field input files and private output path")
    result=audit(read_csv(a.schedule,OUTPUT_COLUMNS),
       read_csv(a.execution,EXEC_HEADER,allow_empty=True),
       read_csv(a.mother_followup,FOLLOWUP_HEADER,allow_empty=True),
       read_csv(a.seed_ledger,LEDGER_HEADER,allow_empty=True))
    save_private(result,a.output_private)
    print("ONISHIBARI_OFFSPRING_LIFE_STAGE_V49",json.dumps({
       "status":"audited_longitudinal_private_sources",
       "n_seed_ids":result["source_unique_identified_seed_records"],
       "n_censored_not_dead":result["missing_or_censored_seed_histories"],
       "not_lifetime_survival_probability":True,
       "no_source_identifiers_logged":True},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
