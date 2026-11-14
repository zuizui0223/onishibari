#!/usr/bin/env python3
"""Pre-register a small randomized Onishibari pollination pilot from REAL bud inventory.

Pure standard library. NO invented field data, NO inferred sex, NO records of
actual pollination. Never put rare plant coordinates in public GitHub.

Usage:
  python3 analysis/randomize_2027_pilot.py \
      --inventory data/private/2027_bud_inventory.csv \
      --output data/private/2027_randomized_assignments.csv --seed 20271114
  python3 analysis/randomize_2027_pilot.py --self-test
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
from pathlib import Path
import random
import tempfile

MORPHS = {
    "male_biased_morph": (
        "open", "bag_no_hand", "bag_sham", "bag_hand_self", "bag_hand_cross"
    ),
    "female_morph": (
        "open", "bag_no_hand", "bag_sham", "bag_hand_cross", "open_plus_cross"
    ),
}
COLUMNS = (
    "population_id", "genet_id", "branch_id", "flower_id", "sex_morph",
    "voucher_checked", "flower_stage"
)
SCHEDULE_COLUMNS = COLUMNS + (
    "treatment", "allocation_status", "outcome_cohort", "cross_donor_genet_id",
    "self_pollen_source_genet_id", "protocol_seed", "protocol_version",
    "not_actual_pollination"
)
BLOCKING_WARNING = (
    "Five planned arms per confirmed plant morph require >=5 predehiscence "
    "buds on a uniquely identified genet. All output is PRE-ASSIGNMENT."
)
VERSION = "onishibari_2027_core_v1"

def read_inventory(path):
    with Path(path).open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None or not set(COLUMNS).issubset(reader.fieldnames):
            raise ValueError("Missing required columns: " + ", ".join(COLUMNS))
        rows = [{name: str(row.get(name, "")).strip() for name in COLUMNS}
                for row in reader]
    if not rows:
        raise ValueError("No inventory rows; no fake field entries are generated")
    seen = set()
    morphs = defaultdict(set)
    for r in rows:
        for name in ("population_id", "genet_id", "branch_id", "flower_id"):
            if not r[name]:
                raise ValueError(f"Empty required identifier: {name}")
        key = tuple(r[x] for x in ("population_id","genet_id","branch_id","flower_id"))
        if key in seen:
            raise ValueError("Duplicate plant/branch/flower identifier: "+repr(key))
        seen.add(key)
        if r["sex_morph"] not in (*MORPHS, "undetermined"):
            raise ValueError("Invalid sex morphology; photos do not prove functional sex")
        if r["voucher_checked"] not in ("yes","no"):
            raise ValueError("voucher_checked must be yes/no")
        if r["flower_stage"] not in ("predehiscence","open","not_ready"):
            raise ValueError("flower_stage must be predehiscence/open/not_ready")
        morphs[(r["population_id"],r["genet_id"])].add(r["sex_morph"])
    if any(len(x)>1 for x in morphs.values()):
        raise ValueError("Same genet has conflicting morphological sex labels")
    return rows

def stable_seed(seed, population_id, genet_id):
    s=f"{seed}|{population_id}|{genet_id}|{VERSION}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(s).digest()[:8], "big")

def allocate(rows, seed):
    if not isinstance(seed, int):
        raise ValueError("Protocol seed must be integer")
    groups=defaultdict(list)
    for r in rows:
        groups[(r["population_id"],r["genet_id"])].append(r)
    result=[]
    for k in sorted(groups):
        group=groups[k]
        morph=group[0]["sex_morph"]
        ready=[r for r in group if r["flower_stage"]=="predehiscence"
               and r["voucher_checked"]=="yes" and morph in MORPHS]
        rng=random.Random(stable_seed(seed,*k))
        arms=MORPHS.get(morph)
        if arms and len(ready)>=len(arms):
            base, spare=divmod(len(ready),len(arms))
            targets={name:base for name in arms}
            extras=list(arms)
            rng.shuffle(extras)
            for name in extras[:spare]:
                targets[name]+=1
            counts=Counter()
            per_branch=defaultdict(Counter)
            plans={}
            by_branch=defaultdict(list)
            for r in ready:
                by_branch[r["branch_id"]].append(r)
            for branch in sorted(by_branch):
                buds=sorted(by_branch[branch],key=lambda r:r["flower_id"])
                rng.shuffle(buds)
                for r in buds:
                    available=[arm for arm in arms if counts[arm]<targets[arm]]
                    if not available:
                        raise AssertionError("Insufficient treatment slots")
                    low=min(per_branch[branch][arm] for arm in available)
                    candidates=[arm for arm in available if per_branch[branch][arm]==low]
                    rng.shuffle(candidates)
                    # Favor underrepresented treatments across the whole genet.
                    min_total=min(counts[arm] for arm in candidates)
                    candidates=[arm for arm in candidates if counts[arm]==min_total]
                    chosen=rng.choice(candidates)
                    plans[(r["branch_id"],r["flower_id"])]=chosen
                    counts[chosen]+=1
                    per_branch[branch][chosen]+=1
            if sum(counts.values())!=len(ready) or set(counts)!=set(arms):
                raise AssertionError("Unbalanced allocation")
            if max(counts.values())-min(counts.values())>1:
                raise AssertionError("Allocation is not balanced within genet")
        else:
            plans={}
        for r in group:
            e={name:r[name] for name in COLUMNS}
            key2=(r["branch_id"],r["flower_id"])
            planned=plans.get(key2)
            if planned is not None:
                status="ASSIGNED_NOT_APPLIED"
            elif r["voucher_checked"]=="no" or morph=="undetermined":
                status="SKIPPED_UNVERIFIED_MORPH"
            elif r["flower_stage"]!="predehiscence":
                status="SKIPPED_NOT_PREDEHISCENCE"
            else:
                status="SKIPPED_LT5_READY_BUDS_PER_GENET"
            e.update({
                "treatment":planned or "",
                "allocation_status":status,
                "outcome_cohort":"mature_seed_followup" if planned else "",
                "cross_donor_genet_id":"",
                "self_pollen_source_genet_id":"",
                "protocol_seed":str(seed),
                "protocol_version":VERSION,
                "not_actual_pollination":"TRUE",
            })
            result.append(e)
    result.sort(key=lambda x:tuple(x[c] for c in
                ("population_id","genet_id","branch_id","flower_id")))
    if len(result)!=len(rows):
        raise AssertionError("An inventory flower was lost")
    return result

def save_rows(rows, path):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=SCHEDULE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

def test():
    with tempfile.TemporaryDirectory() as temp:
        csvfile=Path(temp)/"fixture_synthetic_never_observed.csv"
        with csvfile.open("w",encoding="utf-8",newline="") as out:
            writer=csv.DictWriter(out,fieldnames=COLUMNS)
            writer.writeheader()
            for morph, n, g in [
                ("male_biased_morph",12,"H_FIXTURE"),
                ("female_morph",11,"F_FIXTURE"),
                ("undetermined",7,"U_FIXTURE"),
                ("male_biased_morph",4,"SPARSE_FIXTURE"),
            ]:
                for i in range(n):
                    writer.writerow(dict(population_id="SYNTHETIC",genet_id=g,
                                       branch_id=f"b{i%3}",flower_id=f"flower{i:02d}",
                                       sex_morph=morph,voucher_checked="yes",
                                       flower_stage="predehiscence"))
        rows=read_inventory(csvfile)
        a=allocate(rows,20271114)
        b=allocate(list(reversed(rows)),20271114)
        assert a==b
        assert len(a)==34
        for morph,g in (("male_biased_morph","H_FIXTURE"),
                        ("female_morph","F_FIXTURE")):
            alloc=[r for r in a if r["genet_id"]==g]
            cnt=Counter(r["treatment"] for r in alloc)
            assert set(cnt)==set(MORPHS[morph])
            assert max(cnt.values())-min(cnt.values())<=1
            assert not any(r["cross_donor_genet_id"] for r in alloc)
        assert all(x["allocation_status"].startswith("SKIPPED") for x in a
                   if x["genet_id"] in ("U_FIXTURE","SPARSE_FIXTURE"))
        assert all(x["not_actual_pollination"]=="TRUE" for x in a)
        # Duplicate units cannot be counted as independent flowers.
        with csvfile.open("a",encoding="utf-8",newline="") as out:
            writer=csv.DictWriter(out,fieldnames=COLUMNS)
            writer.writerow(rows[0])
        try:
            read_inventory(csvfile)
        except ValueError:
            pass
        else:
            raise AssertionError("Duplicate inventory row was silently accepted")
    print("PILOT_RANDOMIZATION_TEST=SUCCESS; SYNTHETIC fixture only",
          flush=True)

def main():
    p=argparse.ArgumentParser(description=BLOCKING_WARNING)
    p.add_argument("--inventory", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--seed",type=int,default=20271114)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        test()
        return
    if args.inventory is None or args.output is None:
        p.error("Provide a real bud inventory and an output file (or --self-test)")
    rows=allocate(read_inventory(args.inventory),args.seed)
    save_rows(rows,args.output)
    counts=Counter(x["allocation_status"] for x in rows)
    print("PILOT_RANDOMIZATION=ASSIGNMENT_ONLY; NOT OBSERVATIONS",flush=True)
    print("Output:",args.output, "Inventory:",len(rows),flush=True)
    print("Statuses:",dict(counts),flush=True)
    print("WARNING: donor IDs, pollen QC and treatment outcomes must be filled in "
          "in the field; do NOT interpret planned treatments as executed.",flush=True)

if __name__=="__main__":
    main()
