# Stage reversal in floral reproduction: genotype-matched Sagittaria original-data analysis (v0.10)

Date: 2026-11-16 JST. Status: Original published-source data secondary analysis **completed**, not new field data or a new causal treatment experiment.

## Main inference

Later flowering can either fail because development capacity has been preempted (Brassica rapa), or supply a second opportunity following early failure (Sagittaria trifolia). These two patterns are not a common causal coefficient and do not identify the universal optimal floral dry mass.

## Primary sources and sample definitions

Tang & Dai (2026) Ecology and Evolution, DOI https://doi.org/10.1002/ece3.73955 ; Figshare public original https://doi.org/10.6084/m9.figshare.27610539.v1 . Three published source files were downloaded from the official public Figshare API and uncompressed source SHA256 checked:

- inflorescences.csv, 655 inflorescence records, SHA256 72ff4762ef2b90455bdbb7c39be3c047ef7182e912964b92ae3b5c81f70e71fa;
- total.csv, 109 original individual total-seed records, SHA256 0b24262abe2dfe5520c5f5dd3141eb3c00bbf1acae29716c8ff871fbe5739551;
- plant information.csv, 159 source metadata rows, SHA256 962b712015b62b409c4f42f7306d87414c271982a517675cef5a9a605e7baf11.

The last source is Latin-1 rather than UTF-8 and contains duplicate / incomplete IDs; duplicates with consistent treatment assignment were deduplicated only to establish a treatment-code crosswalk. Treatment site codes: 1 Net, 2 Shade, 5 Open. **Each treatment is completely confounded with its physical location.** Net can change microclimate, light and pollinator access; phenotype changes are not isolated pollination intervention effects.

Four original plants with female inflorescence records lack matched plant-level seed records and were excluded from the final matched analysis. The original inflorescence ID2 label differed from the canonical plant total.csv genotype for one matched individual; the canonical total.csv genotype was used for genotype pairing (mismatch logged). Missing fruit counts are not zero fruit counts. No original plant IDs, sites or geolocations were exported.

Reproducible analysis source: scripts/sagittaria2026_original_sequential_reproduction_v09.py; successful GitHub Actions https://github.com/zuizui0223/onishibari/actions/runs/43793082841 (artifact ID 14193033985, archive ZIP SHA256 b00900b02e67a221449e5ba2b889dac51d99ca5f6eca3c56f83e8bef538d58b0). Resampling by distinct genotype with 5,000 bootstrap draws, paired-genotype sign flips 6,000, random seed 20261116. The paper already studied within-season reproduction; this is new source-verified sensitivity, not discovery of the biological pattern.

## 1. Within individual: first and later inflorescence

First inflorescence = original ordered index xu=1; later = xu>1. For each genet, fruit rate is number of fruits / original **female** flowers, not total male + female flowers, and latter windows are aggregated over later inflorescences. Same-individual paired late-minus-first rate only where both periods have female flowers.

| Site-defined treatment | Plants matched to total seeds | Plants with early+later stages | Mean first-stage fruit fraction, initial cohort | Mean later fraction, paired cohort | **Matched late-minus-first fraction**, genotype-bootstrap 95% |
| --- | ---: | ---: | ---: | ---: | --- |
| Net | 35 | 34 | .1566 | .5319 | **+.3707 [.2263,.5089]** |
| Shade | 34 | 29 | .4426 | .4985 | +.0140 [-.2097,.2324] |
| Open | 40 | 40 | .8685 | .8762 | +.0077 [-.0876,.1221] |

Warning: first-stage descriptive mean includes all 35/34/40 initial-cohort plants; paired late-minus-first means use only 34/29/40 with later opportunities. They cannot be recomputed by subtracting the two displayed columns.

Under Net, among 27 matched plants with zero fruit in first inflorescence, **21 (77.8%)** produced at least one later fruit. Shade 8/13 (61.5%); Open 3/3 (100%) but n=3 too small for broad inference. This is evidence that later floral opportunity can generate fruit after an earlier failure; it does not establish that lost seed numbers are fully replaced or that late flowering was an induced physiological response.

## 2. Within shared clone genotype across treatment sites

There are 47 distinct original source genotypes: 23 occur across all three site-defined treatment groups, 16 across two groups, 8 in one only. Comparing the same genotype between Net and Shade gives the following differences:

| Net minus Shade contrast | Shared genotypes | Estimated mean difference | Genotype-paired bootstrap 95% | Sign-flip p |
| --- | ---: | ---: | --- | ---: |
| Entire plant total seeds | 28 | −24.86 seeds | [−2279.56,+2269.45] | .9842 |
| Inflorescences per plant | 28 | +1.143 | [+.500,+1.821] | .0035 |
| First inflorescence fruit fraction | 28 | −.3061 | [−.5094,−.1040] | .0105 |
| Later inflorescence fruit fraction | 23 | +.1106 | [−.0943,+.3168] | .3136 |
| **Within-genet later-minus-first fraction improvement** | 23 | **+.4784** | **[+.2073,+.7641]** | **.0047** |

Important interpretation: a non-significant plant total seed difference between Net and Shade is **NOT statistical equivalence, full biological compensation or equal yield**; the CI is broad, and all site/environment effects are mixed with treatment. On same 29 clone genotypes, **Net minus Open** plant total seed count is instead −16,029.9 [−20,230.7, −12,272.5], showing second-chance fruiting does not restore high-resource/open-site reproductive output.

## 3. Independent opposite pattern in Brassica and other stage controls

Independent original 2026 Brassica rapa data (Harder & Hamnett-Mount, DOI 10.1002/ajb2.70204; repo v0.6) measured 2,201 deliberately hand-cross-pollinated individual flowers nested in 86 plants. Even when artificially pollinated, later opening flowers showed **24.9–41.4 percentage-point lower** fruit rate than earlier flowers, with a baseline plant-clustered fday logit coefficient −.4057 [−.6036,−.2078]. Opening order was not randomized, and resource competition, hormonal carpic dominance and later-flower physiology remain unresolved.

This helps discriminate two non-exclusive candidate mechanisms:
- **Developmental preemption**: earlier fruits or developmental programs restrict late fruit maturation even with adequate compatible pollen.
- **Second-chance opportunity**: if the early reproductive episode fails and later flowers remain receptive, later flowering can rescue some fruit production.

Silene dioica 2026 original two experiments (Joffard et al., DOI 10.1002/ajb2.70276) separately demonstrate steep loss of pollen deposited across sequential bumblebee visits and no detectable effect on father-diversity profiles from manipulating female floral display. Do not multiply independent experiment coefficients or infer the same pollen/seed individual causal chain. Clarkia 2025 4-species fruit-level pollen supplementation effects and Veratrum 2025 whole-plant-versus-flower supplementation effects further show that fitness conversion depends on recipient/plant and observation grain.

## 4. What remains unproved

The focal Daphne pseudomezereum's functional sex morph, pollen-compatibility status, floral construction mass, paternal genetic siring and lifetime maternal recruitment are **not** directly measured. None of the quoted published experiments estimates a universal mg of floral biomass → lifetime successful recruits relationship. A true causal test must jointly manipulate display attraction, independent compatible-pollen opportunity, early fruit developmental load and floral expenses at the genet level, monitor real female seed and paternal siring success plus future flowering/survival, and replicate across independent phylogenetic clades. Protect rare wild Daphne with the preexisting limit ≤50% confirmed buds and ≤10 manipulated/destructively sampled buds/genet.

**v0.10 conclusion:** later flowers should be viewed as conditional options, not invariably waste or invariably insurance. Which stage has positive returns can switch in time and environment, and maternal fruit-set, whole-genet seeds, actual floral dry-mass expense and future viability are distinct estimands.
