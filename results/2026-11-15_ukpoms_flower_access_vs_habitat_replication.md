# Independent UK PoMS evidence: flower units and landscape have complementary visitor-prediction information

Date 2026-11-15. Status **actual public records downloaded and analyzed**, source/quality, spatial heldout, year heldout and cluster bootstrap performed. **Not an NDVI–fitness or biomass–sex-system association.**

## Provenance and rights
UK Pollinator Monitoring Scheme (2026), *Flower-Insect Timed Count data from the UK Pollinator Monitoring Scheme, 2017–2023*. NERC EDS Environmental Information Data Centre. DOI https://doi.org/10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618
Official ZIP https://data-package.ceh.ac.uk/data/d706088c-9d48-4be5-b156-6ed6ad5e0618.zip
Archive SHA256 `0c0eb23f7d9df81ed112b00cf0568020a3865abed624a83256777ba84a04eafb`.
Use is under Open Government Licence with the source's required PoMS/UKCEH/JNCC partner acknowledgement. Public FIT Counts and allocated 1-km-square FIT Counts are **different survey streams**.
Survey grain: 10 min within a 50 cm × 50 cm quadrat observing a particular target flower type; **not** pollen grain delivery, seed set, population biomass or fitness.

Official archived tables:
- Public voluntary: 13,457 rows, 13,457 unique IDs, 2,971 grid squares, 148,350 recorded insect visitors, 774 zero-total (5.75%).
- Allocated 1km squares (primary inference): 3,053 rows, 3,053 unique IDs, 101 grid squares, 31,259 visitors, 271 zero-total (8.88%). Five with missing positive floral unit count were excluded from models → **3,048** model cases, 101 squares.
- In every row the sum of 10 insect-group columns equals all_insects_total. No invalid date/year mismatch. The 1km panel contains 2,932 British OSGB records and 121 Irish OSIE records. Reference precision 1 km, not a precise plant geolocation.
- Floral units differ: **flower head**, individual flower, spike, umbel. Do not collapse into one biological 'number of flowers' without target/structure controls.
- Categorical `flower_cover` indicates within-patch target flower share; `flower_context` says isolated, same-flower patch or mixed-flower patch. `land_cover` and `habitat_type` are source-recorded habitat classes, **not NDVI or woody AGB**.
- Useful primary data source audit: [run 43703098855](https://github.com/zuizui0223/onishibari/actions/runs/43703098855). Source content checksum file included, raw coordinates/recorder IDs not re-published.

## Primary model: out-of-location validity, 2017–2023
Outcome: log1p visitor count per 10-minute observation; targets: all insects, bumblebees, hoverflies, all bees (bumblebees+honeybees+solitary_bees).
Predictor backbone: target flower identity/type, flower structure/floral unit, year/month/hour, sunshine/cloud/wind.
Feature blocks:
- Habitat: `land_cover` + `habitat_type`
- Floral availability: `log(floral_unit_count)`, `flower_cover`, `flower_context`.
Model: Ridge `alpha=20` on log1p visits, fixed median imputation and one-hot encoding. 5-fold `GroupKFold` on **101 whole 1km-square identifiers**; no site appears in both train/test of a fold; all 2017–2023 in group.

| Outcome | baseline MAE(log1p) | + habitat | + floral | + both | %-improvement vs baseline (habitat, floral, both) |
|---|---:|---:|---:|---:|---|
| All insects | 0.696944 | 0.678016 | 0.665378 | **0.655902** | **2.72%, 4.53%, 5.89%** |
| Hoverflies | 0.603997 | 0.591083 | 0.589974 | **0.580894** | **2.14%, 2.32%, 3.83%** |
| All bees | 0.524255 | 0.511406 | 0.507701 | **0.500488** | **2.45%, 3.16%, 4.53%** |
| Bumblebees | 0.377610 | 0.374675 | 0.374014 | **0.371938** | **0.78%, 0.95%, 1.50%** |

These are **predictive error reductions**, NOT causal effect sizes, percentage changes in abundance, or adaptation. Nested features and differing number of parameters are evaluated only by held-out squares.

### Ablation: flower quantity versus patch availability
Relative MAE improvement over fixed backbone:

| New covariates | All visits | Hoverflies | All bees | Bumblebees |
|---|---:|---:|---:|---:|
| Floral-unit count only | +4.23% | +1.17% | +2.70% | +0.67% |
| Flower cover + surrounding flower context | +2.68% | +2.22% | +1.99% | +0.75% |
| Both floral components | +4.53% | +2.32% | +3.16% | +0.95% |

Important caution: 'number of floral units' can increase visitor **counts per quadrat** even when probability of a visit *per individual flower* does not increase. Units flower-head, individual-flower, spike and umbel vary across target flower types, which were controlled as categories but not proved commensurate.

### Uncertainty from cluster resampling (2000 replicates, conditional on OOF predictions)
Resampled whole **101 squares**, not each timed count. CI is a distribution of *differences in model errors*, conditional on fitted OOF predictions, **NOT** a refitted model confidence interval / causal inference:
- All insects, floral model MAE vs habitat-only MAE: improvement `+0.012639` log1p MAE units; 95% resampled interval **[-0.004257, 0.029439]**, so **DO NOT declare floral globally better than habitat on this evidence alone**.
- All insects, both vs habitat-only: `+0.022115` (CI **[0.014520, 0.030463]**) → complementarity is more reliable.
- Hoverflies, both vs habitat-only: `+0.010189` (CI **[0.005068, 0.015118]**); both vs floral-only `+0.009080` (CI **[0.002042, 0.016597]**).
- All bees, both vs habitat-only: `+0.010918` (CI **[0.006136, 0.015525]**). Bumblebees weaker and uncertain; never pool groups as if all respond equally.
Source **final bootstrap**: https://github.com/zuizui0223/onishibari/actions/runs/43704056334, artifact `ukpoms-1km-flower-vs-landcover-heldout-prediction`, artifact ID **14157521298**, ZIP SHA256 `d5d028a5cbc7d115ffb3fa607b80b91b63cdaa1a76f8bef95a114e87b8b1fe3f`.

## Independent-year validation: train 2017–2022; test all 2023
Training 2,497 records; 2023 testing 551 records, 90 different 1-km square IDs. Repeated square IDs **are allowed across train/test by year**; temporal holdout is NOT geographically independent.

| Outcome | +habitat MAE gain | +floral MAE gain | +both MAE gain (relative to fixed backbone) |
|---|---:|---:|---:|
| All insects | +2.01% | **+6.60%** | **+6.68%** |
| Hoverflies | +1.27% | +1.84% | **+2.22%** |
| All bees | +0.16% | +1.99% | +1.26% |
| Bumblebees | **-3.34%** | **-0.85%** | **-4.00%** |

For all visits, the combined value extends to a heldout new year. For bumblebees, MAE did not improve, showing guild-specific non-portability. Only one year held out; need sequential-year backtesting and geographic transfer across habitats.

Source [run 43703845929](https://github.com/zuizui0223/onishibari/actions/runs/43703845929), artifact `ukpoms-2023-year-holdout-resources-vs-habitat`, ID **14156559787**, ZIP SHA256 `f6cfe3494fd6c8a5f1f697d7ee6fd3bada2a9405fc139c248f14910533ef0374`.

## Cross-system synthesis with prior Zurich 2016
In Zurich 24 gardens, MOD13Q1 NDVI (2016 same flowering season) was initially associated with carrot hoverfly visits (rho+.659) but correlated very tightly with **imperviousness** (rho-.927). Adjusting urban cover and local floral richness eliminated the raw NDVI independent association (partial rho-.155); adding NDVI worsened LOOCV MAE. [Run 43702182669](https://github.com/zuizui0223/onishibari/actions/runs/43702182669).
**The new independent UK 1km PoMS result is about recorded flower quantity / patch availability and categorized land cover; it does NOT yet establish satellite NDVI predictive effects.** Contrasting Zurich and UK data suggests studying quantity × accessibility rather than substituting greenness for floral resources.

## Next valid remote-sensing bridge
PoMS `x1km_square` is an OSGB 2-letter+4-digit or Irish OSIE 1-letter+4-digit **1 km square**, not exact coordinates. Derive grid geometry in EPSG:27700 or EPSG:29903 and average NDVI/EVI across its full square or 1 km buffer rather than extract one 250-m pixel from the square center. Match survey **actual date** to image composite period, with QA flags; use ~2017–2023 MOD13Q1/AVHRR for visit period and explicitly test NDVI > habitat+flower units in spatial GroupKFold. Local floral abundance is directly observed; image greenness does not identify blossoms or biomass. GEDI L4B 2019–2023 landscape AGB is forest AGB and will not measure individual flowers or herbaceous pollinator resources; include only with compatible period and pixel uncertainty.

## Remaining causal and evolutionary gaps
- PoMS observes insect landings, not effective conspecific pollen deposition, mating compatibility, fruit/seed set, survival or sex function.
- No mass/biomass data in PoMS. Individual dry biomass analysis BAAD×Wang in the parent PR has only 91 height+mass+sex species across very sparse within-family sex-system contrasts; not enough to claim sexual system–biomass evolution.
- Zero observations depend strongly on sampling protocol: PoMS fixed-quadrat 10-minute zero rate ~6–9%, Zurich hourly garden windows ~94% zero; these cannot be interpreted as different country-level insect scarcity without protocol harmonization.
- Species/target flower source labels including "Other" mix taxonomic identities; sampling habitats remain observer-selected even in allocated squares.
- These results are **independent ecological prediction tests**, not novelty claims for pollinator scarcity theory, species diversification or adaptive strategy causation.

## Source code
- `scripts/ukpoms_fitcounts_inventory.py` and `scripts/ukpoms_fitcounts_profile.py`: actual EIDC ZIP and grain auditing
- `scripts/ukpoms_flower_vs_landcover_cv.py`: spatially blocked predictions, feature ablation, 101-square cluster error bootstrap
- `scripts/ukpoms_temporal_holdout.py`: future-year test
- `scripts/ukpoms_2023_modis_1km.py`: satellite season-matched NDVI **pilot script, status must be checked separately**; do not claim results if its workflow hasn't passed


## Fully executed remote-sensing bridge: NASA NDVI of 2023 UK PoMS 1km squares
The MOD13Q1 V061 seasonal context pilot **succeeded**, independent of 2016 Zurich, with actual NASA satellite imagery, field-year match and QA.

Source [GitHub Action run 43704350736](https://github.com/zuizui0223/onishibari/actions/runs/43704350736), artifact **14158485341** `ukpoms-2023-1km-NDVI-plus-floral-habitat-QA`, ZIP SHA256 `2eef04aa75f5a2ca60fc50cedfa6a4a8f2cb6bc9f5eeca26088027a0c8a94578`.
Script `scripts/ukpoms_2023_modis_1km.py`.
- Reconstructed 101 square geometries from OSGB/OSIE 1 km grid references (different national projections). For each cell, 4 interior quarter-square points were projected to MODIS sinusoidal coordinates, and NDVI values with pixel reliability QA 0/1 were averaged; **not one 250m point pretending to represent a precisely geocoded 1km survey**.
- 102 MOD13Q1 UK May–August 2023 STAC candidate scenes, chose 8 tile/month scenes, near May 9, June 10, July 12, August 13 (four seasonal representative composites).
- 99 grid squares had valid NDVI for ≥2 months, and 88 visited squares / **547 actual FIT records in 2023** passed model joins; the median per-site seasonal NDVI in these records was 0.75412. Original 2023 1km survey had 551 recorded FITs. No site coordinates, exact national-grid reference or observer identifiers are re-published.
- Model/validation: restrict to those 547 surveys, five-fold held-out 1 km grid squares; compare `recorded target flower + season/hour/weather + flower units/cover/context + landcover/habitat` against **exact same features plus** `2023 May–Aug mean NDVI` and `four-month representative-image NDVI amplitude`.

| Response | Baseline OOF MAE log1p | +NDVI OOF MAE log1p | Prediction MAE gain from adding NDVI |
|---|---:|---:|---:|
| All 10-minute insect visitors | **0.653702** | 0.656527 | **−0.43215% (worse)** |
| Hoverfly visits | **0.590383** | 0.597562 | **−1.21608% (worse)** |
| All bee visits | **0.516447** | 0.519174 | **−0.52804% (worse)** |
| Bumblebee visits | 0.398691 | **0.398011** | **+0.17052% (tiny improvement)** |

**Do not claim NDVI universally has no ecological effect**: this tests one year, small 1km-square sample and four QA-selected composites, not visit-date-matched bloom NDVI, flower abundance, GPP or woody AGB; candidate NDVI highly collinear with recorded habitat is possible but not directly estimated for UK. It does indicate that there was **no meaningful, consistent incremental predictive benefit** of these 2023 NDVI summaries once directly recorded floral access and habitat variables were known. The earlier 24 Zurich gardens independently showed NDVI covaried with urban sealed area (rho=−0.927), and it added no independent skill after that control; **do not generalize urban collinearity numerically to UK landscapes**.

Source of satellite processing is NASA `MOD13Q1.061`, 16-day 250m Terra NDVI (scale 0.0001; pixel QA 0/1), published catalogue https://doi.org/10.5067/MODIS/MOD13Q1.061 . Seasonal mean across a few composite scenes is environmental context, not biomass or flowering stage.

**What does explain variation (prediction, not cause)?** In 2017–2023 allocated 1km surveys, adding recorded flower availability and habitat jointly to target-flower/season/weather improved held-out visitors by ~5.9% MAE(log1p). Decompositions show target-floral-unit count more useful for total/bee visits; patch cover/context more useful for hoverflies in this dataset. The cluster bootstrap interval for floral-only vs habitat-only MAE improvement in all insects includes zero: the data do NOT support declaring floral information universally superior to habitat. Evidence supports *complementary*, group-conditional information and a limit on greenness-only proxy inference.

### Additional source-integrity point
UK PoMS `1km_square` metadata has 2,932 OSGB/121 OSIE records, whose grid references have distinct 2-letter+4-digit vs 1-letter+4-digit shapes. This is why a cell-average rather than fine point pixel is essential. Source audit run 43703595223.
