#!/usr/bin/env python3
"""v0.15: Sensitivity of inferred post-vector pollen loss to unmeasured
abiotic/non-vector removal in the published Johnson & Harder 2023 species table.

This is a bounded scenario analysis, not measured abiotic pollen loss.
It does NOT infer visit frequency, pollen compatibility, fertilization,
genetic siring, or floral biomass.
"""
from __future__ import annotations
import hashlib, io, json, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

URL = "https://ndownloader.figshare.com/files/42364722"
SHA256 = "e5d981d0a7b3804a7d2309b2df6e1418b86a2d9e52d850175fe69ab5a4144a9c"
OUT = Path("analysis_results_pollen_nonvector_sensitivity_v15")
GRID = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]
BOOT = 2500
RNG = np.random.default_rng(20261116)

def original():
    req = urllib.request.Request(URL, headers={"User-Agent": "Original-source-nonvector-pollen-sensitivity-v15"})
    with urllib.request.urlopen(req, timeout=70) as response:
        b = response.read(500001)
    if hashlib.sha256(b).hexdigest() != SHA256:
        raise ValueError("Original species-median file hash mismatch")
    d = pd.read_excel(io.BytesIO(b), sheet_name="Species medians")
    assert len(d) == 228
    assert {"removpcnt", "PTE", "Poltype", "familyGBIF", "Family"} <= set(d.columns)
    d["removpcnt"] = pd.to_numeric(d["removpcnt"], errors="coerce")
    d["PTE"] = pd.to_numeric(d["PTE"], errors="coerce")
    assert d.removpcnt.notna().sum() == 202
    assert d.PTE.notna().sum() == 224
    v = d.dropna(subset=["removpcnt", "PTE"]).copy()
    v["family"] = v.familyGBIF.fillna(v.Family).astype(str)
    v["RE"] = v.removpcnt / 100.0
    v["PTEf"] = v.PTE / 100.0
    v["RPE"] = v.RE * v.PTEf
    v["not_removed"] = 1.0 - v.RE
    v["observed_removed_not_received"] = v.RE - v.RPE
    assert len(v) == 198 and v.family.nunique() == 20
    assert np.all(v.RE.between(0, 1)) and np.all(v.PTEf.between(0, 1))
    assert np.allclose(v.not_removed + v.observed_removed_not_received + v.RPE, 1.0)
    return v

def round6(x):
    return None if not np.isfinite(x) else round(float(x), 6)

def scenarios(v):
    out = []
    for frac in GRID:
        z = v.copy()
        # q is hypothetical unobserved non-vector loss measured as fraction
        # of ALL observed removed pollen, capped by pollen not seen on stigmas.
        # This reclassifies attribution; it cannot change physical grain fates.
        q = np.minimum(frac * z.RE.to_numpy(), z.observed_removed_not_received.to_numpy())
        z["assigned_nonvector"] = q
        z["before_vector"] = z.not_removed + q
        z["after_vector"] = z.observed_removed_not_received - q
        assert np.all(z.after_vector >= -1e-10)
        assert np.allclose(z.before_vector + z.after_vector + z.RPE, 1)
        z["vector_loss_dominates"] = (z.after_vector > z.before_vector).astype(int)
        # Equal-family statistic is a sensitivity only, not a global angiosperm estimand.
        fam = z.groupby("family")["vector_loss_dominates"].mean().to_numpy()
        boot = np.mean(fam[RNG.integers(0, len(fam), (BOOT, len(fam)))], axis=1)
        out.append({
            "assumed_nonvector_fraction_of_reported_removal": frac,
            "species": int(len(z)),
            "families": int(z.family.nunique()),
            "n_vector_loss_dominates": int(z.vector_loss_dominates.sum()),
            "species_equal_fraction_vector_loss_dominates": round6(z.vector_loss_dominates.mean()),
            "family_equal_fraction_vector_loss_dominates": round6(fam.mean()),
            "family_resampling_ci95": [round6(x) for x in np.quantile(boot,[.025,.975])],
            "mean_percent_unremoved_plus_nonvector": round6(100*z.before_vector.mean()),
            "mean_percent_vector_loss": round6(100*z.after_vector.mean()),
            "mean_percent_conspecific_stigma_arrival": round6(100*z.RPE.mean()),
            "mean_assigned_nonvector_percent_of_produced": round6(100*z.assigned_nonvector.mean())
        })
    return out

def thresholds(v):
    # With assigned q=f*RE and before=1-RE+q, after=RE-RPE-q:
    # after>before iff f < (2 RE - 1 - RPE)/(2 RE).
    fcrit = ((2*v.RE.to_numpy() - 1 - v.RPE.to_numpy()) /
             (2*v.RE.to_numpy()))
    dominant = (v.observed_removed_not_received > v.not_removed).to_numpy()
    assert np.array_equal(dominant, fcrit > 0)
    a = fcrit[dominant]
    return {
        "n_baseline_dominant": int(dominant.sum()),
        "threshold_nonvector_fraction_of_observed_removal_if_initially_dominant": {
            "minimum": round6(np.min(a)) if len(a) else None,
            "q25": round6(np.quantile(a,.25)) if len(a) else None,
            "median": round6(np.median(a)) if len(a) else None,
            "q75": round6(np.quantile(a,.75)) if len(a) else None,
            "maximum": round6(np.max(a)) if len(a) else None,
        },
        "assumed_q_is_hypothetical_not_empirically_measured": True
    }

def main():
    OUT.mkdir(exist_ok=True)
    v = original()
    monads = v[v.Poltype.astype(str).str.lower().eq("monads")].copy()
    assert len(monads) == 34 and monads.family.nunique() == 18
    report = {
        "status": "source_SHA_verified_counterfactual_stage_attribution",
        "source_doi": "10.1098/rspb.2023.1148",
        "source_file_doi": "10.6084/m9.figshare.24146916.v1",
        "source_sha256": SHA256,
        "n_source_species": 228,
        "n_matched_RE_PTE": len(v),
        "n_monads": len(monads),
        "n_monad_families": monads.family.nunique(),
        "estimand": "unremoved_plus_hypothetical_nonvector versus vector_removed_not_delivered; both relative to produced pollen",
        "physical_post_anther_removal_fraction_is_unchanged_by_reallocation": True,
        "groups": {}
    }
    for k,group in [("all_matched", v), ("monads",monads),
                    ("non_monads",v[~v.index.isin(monads.index)])]:
        report["groups"][k] = {
            "scenarios": scenarios(group),
            "flip_thresholds": thresholds(group)
        }
    assert report["groups"]["all_matched"]["scenarios"][0]["n_vector_loss_dominates"] == 104
    assert report["groups"]["monads"]["scenarios"][0]["n_vector_loss_dominates"] == 32
    for group in report["groups"].values():
        seq = [p["n_vector_loss_dominates"] for p in group["scenarios"]]
        assert all(a>=b for a,b in zip(seq,seq[1:])), seq
        for s in group["scenarios"]:
            assert abs(sum(s[key] for key in [
                "mean_percent_unremoved_plus_nonvector",
                "mean_percent_vector_loss",
                "mean_percent_conspecific_stigma_arrival"
            ]) - 100.0) < .00001
    fp = OUT / "johnson_harder2023_nonvector_attribution_v15.json"
    fp.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    md = [
        "# Unobserved non-vector pollen removal: conditional source sensitivity (v0.15)",
        "",
        "Original author source: Johnson & Harder (2023), DOI 10.1098/rspb.2023.1148; "
        "species-median table DOI 10.6084/m9.figshare.24146916.v1. "
        "SHA256 verified. This is a **counterfactual attribution audit, not measured abiotic loss**.",
        "",
        "For each species, removal efficiency RE, transfer efficiency PTE, "
        "stigma delivery RPE=RE×PTE. Reclassify hypothetical non-vector removed pollen "
        "q=min(f×RE,RE−RPE), where f is a chosen share of the originally reported "
        "anther removal. Before vector entry =1−RE+q; after vector entry but "
        "not delivered=RE−RPE−q; stigma arrival=RPE. All sum to 100%.",
        "",
        "The physical location of pollen after leaving the anther is unchanged; "
        "this analysis asks **whether the loss can be attributed to pollinator transport**.",
        "",
        "| Assumed non-vector share of reported removal | Monads n=34: loss mainly after vector entry | All matched n=198: same |",
        "| ---: | ---: | ---: |"
    ]
    for i,f in enumerate(GRID):
        m = report["groups"]["monads"]["scenarios"][i]
        a = report["groups"]["all_matched"]["scenarios"][i]
        md.append(f"| {f:.0%} | {m['n_vector_loss_dominates']}/34 "
                  f"({m['species_equal_fraction_vector_loss_dominates']:.1%}) | "
                  f"{a['n_vector_loss_dominates']}/198 "
                  f"({a['species_equal_fraction_vector_loss_dominates']:.1%}) |")
    md.extend(["",
        "**Limits:** no independent measurement of abiotic pollen loss; all scenarios "
        "are assumptions. Distinct source taxonomic lineages are not a random global "
        "sample; orchid and milkweed dominance remains. The source stigma-pollen metric "
        "includes potentially self-derived or incompatible grains and is not "
        "fertilization, seeds, paternal siring or floral construction biomass.",
        "",
        "The threshold q*=max[0,(2RE−1−RPE)/2], in units of produced pollen, "
        "is the minimum reclassification that would eliminate an initially larger "
        "post-vector than pre-vector loss. The corresponding fraction of original "
        "anther removal is q*/RE. This is a classification tipping point, not "
        "an estimate of actual wind/rain/grooming losses.",
        "",
        "No *Daphne pseudomezereum* sex-function, pollen compatibility, flower-mass "
        "or reproductive fitness is measured here.",
        ""])
    (OUT/"johnson_harder2023_nonvector_attribution_v15.md").write_text("\n".join(md))
    print("MONAD_ABIOTIC_V15",json.dumps({
        "status":report["status"],
        "monad_threshold":report["groups"]["monads"]["flip_thresholds"],
        "monad_scenarios":report["groups"]["monads"]["scenarios"],
        "all_scenarios":report["groups"]["all_matched"]["scenarios"]
    }),flush=True)

if __name__ == "__main__":
    main()
