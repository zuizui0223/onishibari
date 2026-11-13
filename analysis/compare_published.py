#!/usr/bin/env python3
"""Recalculate descriptive within-study contrasts from published Daphne aggregates.

No onishibari observations, no statistical tests, no pooled meta-analysis.
Run: python3 analysis/compare_published.py --output docs/comparative_results.md
"""
import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived" / "published_aggregates.csv"
CONTRASTS = [
    ("jezoensis_2017", "flower_width_mm", "hermaphrodite", "female"),
    ("jezoensis_2017", "flower_mass_mg", "hermaphrodite", "female"),
    ("jezoensis_2017", "tube_length_mm", "hermaphrodite", "female"),
    ("jezoensis_2017", "hand_fruit_set_proportion", "herm_cross", "female_cross"),
    ("jezoensis_2017", "hand_fruit_set_proportion", "herm_self", "herm_cross"),
    ("miyabeana_2018", "natural_fruit_set_proportion", "open_canopy", "closed_canopy"),
    ("mezereum_2023", "branch_fruit_set_proportion", "hand_outcross", "ambient"),
    ("laureola_2007", "plant_stem_count", "female", "hermaphrodite"),
    ("laureola_2007", "flower_count", "female", "hermaphrodite"),
]

def load_data(path=DATA):
    with path.open(newline="", encoding="utf-8") as handle:
        data = list(csv.DictReader(handle))
    expected = {"species", "study", "metric", "group", "value", "doi"}
    if not data or not expected.issubset(data[0]):
        raise ValueError("Published data table is missing required columns")
    keyed = {}
    for row in data:
        if row["species"] == "Daphne pseudomezereum":
            raise ValueError("D. pseudomezereum has no verified measurements in this evidence table")
        if not row["doi"] or not row["study"] or not row["group"]:
            raise ValueError("Missing source metadata")
        k = (row["study"], row["metric"], row["group"])
        if k in keyed:
            raise ValueError(f"Duplicate metric: {k}")
        row["value"] = float(row["value"])
        if "proportion" in row["metric"] and not 0 <= row["value"] <= 1:
            raise ValueError(f"Invalid proportion: {k}")
        keyed[k] = row
    return keyed

def comparison(data, study, metric, numerator, denominator):
    a = data[(study, metric, numerator)]
    b = data[(study, metric, denominator)]
    if a["species"] != b["species"] or a["doi"] != b["doi"]:
        raise ValueError("Cannot compare different species or source studies")
    if b["value"] == 0:
        raise ValueError("Zero comparison denominator")
    ratio = a["value"] / b["value"]
    difference = a["value"] - b["value"]
    if "proportion" in metric:
        difference *= 100
        unit = "percentage points"
    else:
        unit = metric.rsplit("_", 1)[-1]
    return a["species"], study, metric, f"{numerator} / {denominator}", ratio, difference, unit, a["doi"]

def make_report(data):
    out = [
        "# Daphne published aggregate comparisons",
        "",
        "**Evidence only: these are NOT data from onishibari (D. pseudomezereum).**",
        "Each ratio uses two mean values from the same study and species; no effect-size uncertainty, pooled estimate, hypothesis test, or causal claim.",
        "",
        "| Species | Study | Trait | Contrast | Ratio | Difference | DOI |",
        "|---|---|---|---|---:|---:|---|",
    ]
    for args in CONTRASTS:
        species, study, metric, contrast, ratio, difference, unit, doi = comparison(data, *args)
        out.append(f"| *{species}* | {study} | {metric} | {contrast} | {ratio:.3f}x | {difference:+.2f} {unit} | [DOI](https://doi.org/{doi}) |")
    out += [
        "",
        "## Interpretation limits",
        "- No quantitative target-species observations are available in this compilation.",
        "- Self/cross fruit set ratios do not locate the pollen incompatibility or embryo failure mechanism.",
        "- Whole-shrub size, flower dimension, and flower number are distinct traits; their sexual biases may differ.",
        "- D. miyabeana open versus closed canopy is an observational habitat contrast.",
        "- D. mezereum is a 17-plant paired-branch comparison; ratio of means is not a paired effect estimate.",
        "- Quoted +/- errors and n details are retained in the CSV, but do not justify propagated CIs without the original sampling design.",
        "- Some source denominators are unextracted, never substitute floral-morphology n.",
        "- The 1989 Kikuzawa stigma loads are not part of this numerical comparison: use the primary-study abstract cautiously.",
        "",
        "## Next discriminating target experiment",
        "Randomize intact/sham versus emasculated flowers and standardize other-pollen input. Add self-first/cross-first pollen treatments; separately measure stigma grains, pollen tubes, embryos, ripe seeds and germination. A no-fruit bagging control cannot establish self-incompatibility.",
        "",
    ]
    return "\n".join(out)

def self_check(data):
    a = comparison(data, "jezoensis_2017", "flower_width_mm", "hermaphrodite", "female")
    b = comparison(data, "jezoensis_2017", "hand_fruit_set_proportion", "herm_cross", "female_cross")
    c = comparison(data, "jezoensis_2017", "hand_fruit_set_proportion", "herm_self", "herm_cross")
    assert abs(a[4] - 15.76 / 11.10) < 1e-12
    assert abs(b[4] - 0.211 / 0.604) < 1e-12
    assert abs(c[4] - 0.071 / 0.211) < 1e-12
    assert len(CONTRASTS) == 9
    assert all(row["species"] != "Daphne pseudomezereum" for row in data.values())
    print("9 within-study comparisons validated (published aggregates only)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    table = load_data()
    if args.check:
        self_check(table)
    result = make_report(table)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result, encoding="utf-8")
    elif not args.check:
        print(result)
