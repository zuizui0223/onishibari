#!/usr/bin/env python3
"""Reproduce Shibata 2020 thesis Table 3-2 paternity-assignment totals.

No new father assignment, pollen movement estimate or Onishibari genotype.
"""
import csv
from pathlib import Path

DATA=Path(__file__).resolve().parents[1]/"data/derived/shibata2020_paternity_assignment.csv"

def main():
    with DATA.open(encoding="utf-8", newline="") as stream:
        rows=list(csv.DictReader(stream))
    if len(rows)!=4:
        raise ValueError("Expected four population-year entries")
    keys={(r["site"],r["year"]) for r in rows}
    if keys!={("CTS","2015"),("CTS","2016"),("NOP","2015"),("NOP","2016")}:
        raise ValueError("Unexpected paternity design")
    for r in rows:
        if r["taxon"]!="Daphne jezoensis" or r["source_url"]!="https://doi.org/10.14943/doctoral.k14152":
            raise ValueError("Invalid source or taxon")
        if not 0 <= int(r["assigned_father"]) <= int(r["progeny_analyzed"]):
            raise ValueError("Impossible count")
        if not 0 < float(r["mean_parent_distance_m"]) < 100:
            raise ValueError("Distance outside source range")
    analyzed=sum(int(r["progeny_analyzed"]) for r in rows)
    assigned=sum(int(r["assigned_father"]) for r in rows)
    if (analyzed,assigned)!=(814,482):
        raise ValueError("Table 3-2 not reproduced")
    print("SHIBATA2020_PATERNITY_CHECK=SUCCESS"
          f"; analyzed={analyzed}; assigned={assigned}; proportion={assigned/analyzed:.4f}")
    print("NOTE=The remaining 332 progeny were unassigned at the study threshold; "
          "father-mother distances are conditional on assignment, not a complete "
          "pollen dispersal kernel.")

if __name__=="__main__":
    main()
