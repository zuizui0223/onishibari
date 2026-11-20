#!/usr/bin/env python3
"""Audit reported chi-square(1) likelihood-ratio p values in Barbot et al 2022.

Source https://pmc.ncbi.nlm.nih.gov/articles/PMC9709571/
The published text appears to contain p-value decimal transcription errors.
This check derives *only* the p corresponding to reported chi-squared and
df=1. It cannot verify original LRT inputs or fix the original manuscript.
"""
import csv
import json
import math
from pathlib import Path

OUT=Path("analysis_results_silene_stat_audit")
OUT.mkdir(exist_ok=True)

TESTS=[
    ("artificial_display_to_any_visit",1.82,0.17,"approx",False),
    ("artificial_display_to_pollinator_number",6.04,0.014,"approx",False),
    ("artificial_display_to_total_visits",6.91,0.01,"less_than",False),
    ("real_flower_count_to_any_visit",3.25,0.71,"approx",True),
    ("flower_size_PC1_to_total_visits",6.44,0.11,"approx",True),
    ("artificial_display_to_female_mating_success",4.97,0.26,"approx",True),
    ("artificial_display_to_male_reproductive_success",0.87,0.35,"approx",False),
    ("artificial_display_to_female_reproductive_success",1.23,0.27,"approx",False),
]


def sf_chisquare1(x):
    assert x>=0
    return math.erfc(math.sqrt(x/2))


def check():
    assert math.isclose(sf_chisquare1(0),1.0)
    assert abs(sf_chisquare1(3.84145882069)-.05)<2e-9
    assert .013<sf_chisquare1(6.04)<.015
    assert .024<sf_chisquare1(4.97)<.028


def main():
    check()
    rows=[]
    for label,stat,p_reported,relation,known_mismatch in TESTS:
        p_calc=sf_chisquare1(stat)
        reported_consistent=(p_calc < p_reported if relation=="less_than" else abs(p_calc-p_reported)<=max(.006,.15*p_calc))
        rows.append({
           "contrast":label,
           "reported_chi2":stat,
           "reported_df":1,
           "paper_p_text":("<"+str(p_reported)) if relation=="less_than" else str(p_reported),
           "p_recalculated_chi2_1":round(p_calc,7),
           "numerically_consistent":reported_consistent,
           "likely_decimal_transcription_issue":not reported_consistent,
           "checks_underlying_raw_model":False,
           "source":"Barbot_et_al_2022_doi_10.1098/rspb.2022.1987"
        })
        if known_mismatch:assert not reported_consistent,(label,p_calc,p_reported)
        else:assert reported_consistent,(label,p_calc,p_reported)
    assert sum(r["likely_decimal_transcription_issue"] for r in rows)==3
    with (OUT/"silene_chisquare_pvalue_consistency.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]))
        wr.writeheader();wr.writerows(rows)
    res={"status":"success_quoted_statistical_consistency_only",
         "original_study":"Barbot et al. 2022 Proc R Soc B 10.1098/rspb.2022.1987",
         "original_experiment":"43 males and 43 females per group; 2x visual display with dummy flowers; 20 hand-pollinated positive-control females; 2,030 progeny genotyped",
         "found_discordant_quoted_chisquare_and_p_count":3,
         "do_not_reinterpret_underlying_tests_without_raw_stat_models":True,
         "note":"Three reported decimal p values fail the df1 chi-square identity, likely article typographical issues, rather than evidence for a different result.",
         "rows":rows}
    (OUT/"silene_reported_pvalue_audit.json").write_text(json.dumps(res,indent=2)+"\n")
    print(json.dumps(res,indent=2),flush=True)


if __name__=="__main__":
    main()
