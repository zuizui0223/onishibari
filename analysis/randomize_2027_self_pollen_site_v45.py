#!/usr/bin/env python3
"""Optional 2027 Daphne pseudomezereum H-morph self-pollen site pilot v0.45.

Preserves existing v1 randomizer and 5-arm protocol. Does not create any
observed field data. Input: real, verified private flower inventory; output
remains PRIVATE. In H, distinguish within-flower self pollen vs geitonogamous
self pollen from a DIFFERENT same-genet flower, using a DEDICATED donor bud.
The donor bud counts toward the conservation manipulation ceiling and is
NEVER an unmanipulated natural-open female-outcome replicate.

Prior 2008 FFPRI D. jezoensis observation (not Onishibari):
same-flower forced self: no fruits; same-plant different-flower self: fruits.
First-party source: FFPRI Annual Report 2008 No.49, p.40,
https://www.ffpri.go.jp/fsm/research/pubs/nenpo/documents/annual_rep_049.pdf
Pollen source differs by flower age, quantity, and intrinsic contamination;
this is only a mechanistic hypothesis for Onishibari, not demonstrated SI.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import argparse,csv,hashlib,random,tempfile
from pathlib import Path
from randomize_2027_pilot import read_inventory,COLUMNS

VERSION="onishibari_2027_same_vs_other_flower_self_v45"
H_ARMS=(
  "open",                       # natural observation only; do NOT donor-sample
  "bag_no_hand",                # spontaneous selfing + apomixis/contamination
  "bag_sham",                   # brushing / manipulation without pollen
  "bag_self_same_flower",       # brush ripe pollen from the recipient's anthers
  "bag_self_other_flower",      # different pollen-source flower, same genet
  "bag_hand_cross",             # distinct genetically validated donor
  "pollen_donor_same_genet"     # pollen source only, not maternal outcome
)
F_ARMS=(
 "open","bag_no_hand","bag_sham","bag_hand_cross","open_plus_cross"
)
ARMS={"male_biased_morph":H_ARMS,"female_morph":F_ARMS}
OUTPUT_COLUMNS=COLUMNS+(
 "planned_treatment","allocation_status","role",
 "mature_seed_followup_expected",
 "same_genet_donor_branch_id","same_genet_donor_flower_id",
 "cross_donor_genet_id", "donor_genet_verification_status",
 "assigned_protocol_seed","protocol_version","not_actual_pollination"
)

def stable(seed,key):
    return int.from_bytes(hashlib.sha256(
      f"{seed}|{key[0]}|{key[1]}|{VERSION}".encode("utf8")).digest()[:8],"big")

def assign(rows,seed=20271117):
    if not isinstance(seed,int):raise TypeError("seed must be int")
    groups=defaultdict(list)
    for r in rows:groups[(r["population_id"],r["genet_id"])].append(r)
    out=[]
    for key,allrows in sorted(groups.items()):
        morph=allrows[0]["sex_morph"]
        arms=ARMS.get(morph)
        ready=[r for r in allrows if r["voucher_checked"]=="yes" and
                r["flower_stage"]=="predehiscence" and arms]
        # This ceiling counts every assigned H pollen donor, H receiver and F
        # receiver as a manipulation, even if only tagged / sampled.
        maxallowed=min(10,len(ready)//2)
        seed_local=stable(seed,key)
        rng=random.Random(seed_local)
        plans={}
        if arms and maxallowed>=len(arms):
            chosen=rng.sample(sorted(ready,key=lambda x:(x["branch_id"],x["flower_id"])),len(arms))
            ordered=list(arms);rng.shuffle(ordered)
            plans={(r["branch_id"],r["flower_id"]):t
                   for r,t in zip(chosen,ordered)}
            assert len(plans)==len(arms)
            assert set(plans.values())==set(arms)
            if morph=="male_biased_morph":
                donor=next(r for r in chosen if plans[(r["branch_id"],r["flower_id"])]==
                    "pollen_donor_same_genet")
                donor_id=(donor["branch_id"],donor["flower_id"])
                assert donor_id!=next((r["branch_id"],r["flower_id"]) for r in chosen
                  if plans[(r["branch_id"],r["flower_id"])]=="bag_self_other_flower")
            else:donor_id=None
        else:donor_id=None
        for r in allrows:
            id=(r["branch_id"],r["flower_id"])
            arm=plans.get(id)
            if arm:
                status="ASSIGNED_NOT_APPLIED"
                role=("pollen_source_only" if arm=="pollen_donor_same_genet"
                      else "female_maternal_endpoint")
            elif r["voucher_checked"]!="yes" or morph not in ARMS:
                status="SKIPPED_UNVERIFIED_MORPH_OR_VOUCHER";role="none"
            elif r["flower_stage"]!="predehiscence":
                status="SKIPPED_NOT_PREDEHISCENCE";role="none"
            elif arms and len(ready)//2<len(arms):
                status="SKIPPED_INSUFFICIENT_VERIFIED_BUDS_FOR_CAP";role="none"
            else:status="SKIPPED_RESERVED_FOR_NATURAL_REPRODUCTION";role="none"
            if arm=="bag_self_other_flower":
                assert donor_id is not None
                donor_branch,donor_flower=donor_id
            elif arm=="bag_self_same_flower":
                donor_branch,donor_flower=id
            else:donor_branch,donor_flower="",""
            out.append({**{name:r[name] for name in COLUMNS},
                "planned_treatment":arm or "",
                "allocation_status":status,
                "role":role,
                "mature_seed_followup_expected":(
                  "YES" if arm and arm!="pollen_donor_same_genet" else "NO"),
                "same_genet_donor_branch_id":donor_branch,
                "same_genet_donor_flower_id":donor_flower,
                "cross_donor_genet_id":"",  # never invent genotype source
                "donor_genet_verification_status":(
                  "NOT_YET_GENOTYPED" if arm and arm.startswith("bag_self")
                  else "NOT_APPLICABLE"),
                "assigned_protocol_seed":str(seed),
                "protocol_version":VERSION,
                "not_actual_pollination":"TRUE"
            })
        count=sum(r["allocation_status"]=="ASSIGNED_NOT_APPLIED" for r in out
                  if (r["population_id"],r["genet_id"])==key)
        assert count<=10 and count<=len(ready)//2
    out.sort(key=lambda x:tuple(x[col] for col in COLUMNS[:4]))
    assert len(out)==len(rows)
    return out

def validate(rows,plans):
    index={(r["population_id"],r["genet_id"],r["branch_id"],r["flower_id"]):r
        for r in plans}
    assert len(index)==len(rows)
    bygen=defaultdict(list)
    for x in plans:
        assert x["not_actual_pollination"]=="TRUE"
        bygen[(x["population_id"],x["genet_id"])].append(x)
    for key,group in bygen.items():
        orig=[r for r in rows if (r["population_id"],r["genet_id"])==key]
        nready=sum(r["voucher_checked"]=="yes" and
                   r["flower_stage"]=="predehiscence" for r in orig)
        act=[r for r in group if r["allocation_status"]=="ASSIGNED_NOT_APPLIED"]
        assert len(act)<=10 and len(act)<=nready//2
        if not act:continue
        morph=act[0]["sex_morph"]
        assert set(r["planned_treatment"] for r in act)==set(ARMS[morph])
        assert len(act)==len(ARMS[morph])
        if morph=="male_biased_morph":
            donor=next(r for r in act if r["planned_treatment"]=="pollen_donor_same_genet")
            other=next(r for r in act if r["planned_treatment"]=="bag_self_other_flower")
            same=next(r for r in act if r["planned_treatment"]=="bag_self_same_flower")
            assert donor["role"]=="pollen_source_only"
            assert donor["mature_seed_followup_expected"]=="NO"
            assert other["same_genet_donor_branch_id"]==donor["branch_id"]
            assert other["same_genet_donor_flower_id"]==donor["flower_id"]
            assert (other["branch_id"],other["flower_id"]) != (donor["branch_id"],donor["flower_id"])
            assert same["same_genet_donor_branch_id"]==same["branch_id"]
            assert same["same_genet_donor_flower_id"]==same["flower_id"]
            assert all(r["mature_seed_followup_expected"]=="YES"
                       for r in act if r["planned_treatment"]!="pollen_donor_same_genet")
        assert all(r["cross_donor_genet_id"]=="" for r in act)
    return {
       "status":"PASS_PROSPECTIVE_DAPHNE_PISTIL_POLLEN_SOURCE_ALLOCATION_V45",
       "bud_records":len(rows),
       "genets":len(bygen),
       "assigned_genets":sum(any(r["allocation_status"]=="ASSIGNED_NOT_APPLIED"
                              for r in a) for a in bygen.values()),
       "assigned_male_morph_genets":sum(any(r["planned_treatment"]=="pollen_donor_same_genet"
                                      for r in a) for a in bygen.values()),
       "assigned_female_morph_genets":sum(any(r["planned_treatment"]=="open_plus_cross"
                                      for r in a) for a in bygen.values()),
       "assigned_total_plant_bud_roles":sum(r["allocation_status"]=="ASSIGNED_NOT_APPLIED" for r in plans),
       "no_pollen_donor_confounded_with_natural_open_receiver":True,
       "no_cross_pollen_source_or_field_outcomes_imputed":True,
       "strict_conservation_budget_counts_dedicated_same_genet_donor":True,
    }

def save(plans,out):
    out=Path(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",newline="",encoding="utf8") as h:
        w=csv.DictWriter(h,fieldnames=OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(plans)

def synthetic_test():
    with tempfile.TemporaryDirectory() as tmp:
        file=Path(tmp)/"PRIVATE_SYNTHETIC_ONLY.csv"
        with file.open("w",newline="",encoding="utf8") as o:
            w=csv.DictWriter(o,fieldnames=COLUMNS);w.writeheader()
            for sex,genet,n in (
                  ("male_biased_morph","H_OK_14",14),
                  ("female_morph","F_OK_10",10),
                  ("male_biased_morph","H_INSUFFICIENT_12",12),
                  ("female_morph","F_INSUFFICIENT_9",9),
                  ("undetermined","UNKNOWN_14",14)):
                for k in range(n):
                    w.writerow({"population_id":"SYNTHETIC",
                        "genet_id":genet,"branch_id":f"b{k%3}",
                        "flower_id":f"flower{k:02d}",
                        "sex_morph":sex,"voucher_checked":"yes",
                        "flower_stage":"predehiscence"})
        raw=read_inventory(file)
        p=assign(raw,20271117)
        assert p==assign(list(reversed(raw)),20271117)
        s=validate(raw,p)
        assert s["assigned_male_morph_genets"]==1
        assert s["assigned_female_morph_genets"]==1
        assert s["assigned_total_plant_bud_roles"]==12
        h=[r for r in p if r["genet_id"]=="H_OK_14" and r["planned_treatment"]]
        assert len(h)==7
        assert sum(r["role"]=="pollen_source_only" for r in h)==1
        assert sum(r["mature_seed_followup_expected"]=="YES" for r in h)==6
        assert not any(r["planned_treatment"] for r in p
                       if r["genet_id"] in ("H_INSUFFICIENT_12","F_INSUFFICIENT_9","UNKNOWN_14"))
        print("DAPHNE_POLLEN_SITE_PILOT_V45",json_summary(s),flush=True)

def json_summary(x):
    import json
    return json.dumps({**x,"fixture_only_not_observed":True},ensure_ascii=False)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--seed",type=int,default=20271117)
    p.add_argument("--self-test",action="store_true")
    arg=p.parse_args()
    if arg.self_test:
        synthetic_test()
        return
    if not arg.inventory or not arg.output:
        p.error("Supply verified private inventory/output or run --self-test")
    plans=assign(read_inventory(arg.inventory),arg.seed)
    summ=validate(read_inventory(arg.inventory),plans)
    save(plans,arg.output)
    print("DAPHNE_POLLEN_SITE_PILOT_V45",json_summary(summ),flush=True)
    print("SCHEDULE_ONLY_NOT_ACTUAL_POLLINATION; PRIVATE_COORDS_AND_IDS_STAY_OUT_OF_GITHUB")

if __name__=="__main__":main()
