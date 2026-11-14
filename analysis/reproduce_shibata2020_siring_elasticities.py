#!/usr/bin/env python3
"""Check predicted male-fitness elasticities in Shibata (2020) Chapter 3.

NO Onishibari male-fitness observation. A doubling multiplier is computed from
a fitted log-relationship; this is not a causal estimate of flower effects.
"""
import csv
from math import exp, log
from pathlib import Path

DATA=Path(__file__).resolve().parents[1]/"data/derived/shibata2020_sire_success_coefficients.csv"
URL="https://doi.org/10.14943/doctoral.k14152"

def main():
    with DATA.open(encoding="utf-8",newline="") as file:
        rows=list(csv.DictReader(file))
    if len(rows)!=2 or {r["plot"] for r in rows}!={"CTS","NOP"}:
        raise ValueError("Expected two population-specific regression formulas")
    by={r["plot"]:r for r in rows}
    for r in rows:
        if r["taxon"]!="Daphne jezoensis" or r["source_url"]!=URL:
            raise ValueError("Wrong species or source")
    if float(by["CTS"]["log_own_flower_coefficient"])!=.16:
        raise ValueError("CTS flower slope mismatch")
    if float(by["CTS"]["log_local_female_flowers_coefficient"])!=.39:
        raise ValueError("CTS local female density slope mismatch")
    if float(by["NOP"]["log_own_flower_coefficient"])!=.18:
        raise ValueError("NOP flower slope mismatch")
    if by["NOP"]["log_local_female_flowers_coefficient"].strip():
        raise ValueError("No local female term was retained in the NOP regression")
    cts_flow=2**float(by["CTS"]["log_own_flower_coefficient"])
    nop_flow=2**float(by["NOP"]["log_own_flower_coefficient"])
    cts_local=2**float(by["CTS"]["log_local_female_flowers_coefficient"])
    assert 1.11<cts_flow<1.13 and 1.12<nop_flow<1.15
    assert 1.3<cts_local<1.33
    print(f"SHIBATA2020_SIRE_MODEL_CHECK=SUCCESS; "
          f"own_flowers_doubled_CTS_multiplier={cts_flow:.3f}; "
          f"own_flowers_doubled_NOP_multiplier={nop_flow:.3f}; "
          f"local_female_flowers_doubled_CTS_multiplier={cts_local:.3f}")
    print("CAUTION=Predictor associations, conditional on the original regression and "
          "selected parentage; not flower-size causal effects or target-species estimates.")

if __name__=="__main__":
    main()
