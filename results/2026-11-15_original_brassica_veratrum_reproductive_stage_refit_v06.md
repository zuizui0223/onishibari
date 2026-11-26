# Reproductive investment conversion: original-data two-study refit (v0.6, 2026-11-15)

## Scope, new data, limitations

**Status: NEW original-CSV reanalysis COMPLETE**, using only official published Dryad version ZIPs. Public individual-file API requests returned 401, but the documented version-download endpoint gave 200 for all five source studies. No credentials, workarounds, or private data used. The earlier v0.5 data-access blocker is resolved. These models are exploratory secondary analyses, not authors' confirmed causal models or new field work.

Sources:
- Harder & Hamnett-Mount 2026, *Brassica rapa*, https://doi.org/10.1002/ajb2.70204 ; Dryad https://doi.org/10.5061/dryad.hx3ffbgtq ; original published version ID 434166.
- Jiang et al. 2025, *Veratrum grandiflorum*, https://doi.org/10.1002/njb.04873 ; Dryad https://doi.org/10.5061/dryad.9s4mw6mw2 ; published version ID 405186.
- Source inventory also verified original *Silene latifolia* (208 individual-level phenotype/survival records), *S. dioica* (192 phenotypes and 1,721 offspring-to-candidate-father records), and *Clarkia* (430 flower/fruit rows across 4 species); **not** yet fitted in this v0.6 model.
- GitHub Actions original source inventory: https://github.com/zuizui0223/onishibari/actions/runs/43725181321 (SUCCESS). Primary two-study fit: https://github.com/zuizui0223/onishibari/actions/runs/43725981425 (SUCCESS).
- Model code: scripts/primary_brassica_veratrum_stage_refit.py ; workflow: .github/workflows/angiosperm-primary-stage-refit.yml ; source outputs are published as a time-limited workflow artifact.
- Source integrity: compare SHA256 **uncompressed original CSV member contents**, not dynamically regenerated ZIP envelope bytes. Example Brassica ovule-level file SHA 2ccafcb2c3ab73ec2c882c60f4047cbb34315287fcdc296388773dafaeba076a; Veratrum processed-pair file SHA d020dfa7fbdbe05ab822c64004b9de83b9840e07d02d04918678974f34848b7a. Source individual/flower IDs are not exported.

## Brassica: 86 plants, 5,658 total opened flowers, 2,346 fruits

Plant level reconstructed from 304 raceme rows; treatment assigned to **plants**, not to 5,658 statistically independent genets.

| Fraction of flowers deliberately hand outcrossed | Plants | Opened flowers | Mature fruits | Pooled fruits/flower | Average plant fruit fraction |
|---|---:|---:|---:|---:|---:|
| 25% | 21 | 1,541 | 323 | .210 | .223 |
| 50% | 23 | 1,533 | 656 | .428 | .428 |
| 75% | 20 | 1,321 | 639 | .484 | .502 |
| 100% | 22 | 1,263 | 728 | .576 | .605 |

**Correction to the previous narrative:** original authors reported **decelerating improvement**, not an exactly flat 75–100% plateau. In our different exploratory plant-level aggregated-binomial GLM adjusted for log vegetative mass, the contrast **75% versus 100%** remains negative (beta −.3688, HC3 95%CI −.6660 to −.0716, p=.0150). This is not an isolated randomized *flower-biomass* effect.

Vegetative mass coefficient in this model: beta **−.4183**, HC3 95%CI **[−.7064,−.1302]**, p=.00443, odds ratio per doubling **0.7483 [0.6128,0.9137]** for the *fraction* of flowers fruiting, holding pollen-assignment category fixed. Pearson overdispersion ≈3.93. The paper's differently specified mass effect beta −.321±.148 SE does not exactly equal this new model. A declining fruit fraction can coexist with more **total** fruits; vegetative mass is **not** true floral construction mass and mass was not randomized.

## Brassica: 2,201 artificially cross-pollinated individual flowers nested in 86 plants

For each individual plant, divide its *hand-pollinated flowers* into flowering on/before versus after its median anthesis day, and mean plant-specific fruit success (not pooled flower-level pseudo-replicates):

| Pollen treatment | Earlier hand-pollinated fruit fraction | Later hand-pollinated fruit fraction | Later minus earlier | Plant pairs |
|---|---:|---:|---:|---:|
| 25% | .8667 | .5659 | **−30.1 percentage points** | 21 |
| 50% | .9176 | .6690 | **−24.9 percentage points** | 22 |
| 75% | .8818 | .5074 | **−37.4 percentage points** | 20 |
| 100% | .8748 | .4612 | **−41.4 percentage points** | 22 |

Further exploratory binomial regression with **plant-clustered** standard errors and treatment-by-anthesis-day interaction gives base 25%-treatment fday beta **−.4057 [−.6036,−.2078], p≈0.000058**, adjusted for log vegetative mass. The *extra* interaction at 100% is −.2378 [−.4828,+.0071], p=.057; evidence is insufficient for a formal claim that the 100% treatment decline is steeper. Anthesis day was **not randomized**. These data locate frequent failure of **fruit formation after artificial pollen supply**, but cannot separate resource competition from hormonal carpic dominance, phenological/physiological quality, or other mechanisms.

## Brassica: 9,421 ovules from sampled developing fruits nested in 46 different Exp2 plants

Do not confuse **Experiment 2, developing-fruit conditional** data with all 86 plants / all flowers in Experiment 1.

| Ovule fate | Number | Fraction |
|---|---:|---:|
| Fertilized | 5,940 of 9,421 | 63.05% pooled ovule-weighted |
| Not fertilized | 3,481 of 9,421 | 36.95% |
| Developed into a seed *among fertilized* | 5,673 of 5,940 | 95.51% |
| Aborted after fertilization | 267 of 5,940 | 4.49% |

Original article reports about 61% fertilized and 96.3% development, based on its own summary/weighting; pooled ovules are a **different denominator**. These do not represent seeds per originally opened flower.

Resampling **independent plants** rather than independent ovules (4,000 bootstrap draws):
- Mean per-plant ovule fertilization: 25% treatment **0.5837** (n21) vs 75% **0.6507** (n25). Difference +0.0670, 95% bootstrap **[−0.0293,+0.1655]**; no clear difference.
- Mean per-plant development conditional on fertilization: **0.9508** vs **0.9516**, difference +0.00086 with bootstrap CI **[−0.0255,+0.0319]**; conditional development remains high.
- **Source QC flag:** 714/719 fruit-level raw rows satisfy total ovules = unfertilized + aborted + seed; **5 fully recorded rows conflict**, with total minus part sums **−16 (1 row), −4 (1), −2 (3)**. Do NOT silently fill or correct. The original 9,421 ovule-level source was separately audited. Resolve coding notes before a joint per-fruit likelihood.

## Veratrum: study manipulation scope changes measured pollen limitation

Original derived data contain **97 processed pair records**, divided 32 whole-plant, 32 inflorescence, 33 flower-level. Positive-pair log effect means reproduce the authors' exact published values, **without pseudocounting zero**. Flower and inflorescence levels include zero or missing seed ratios (7 and 8 respectively), so the positive-log samples are 26 and 24, versus 32 for whole plants.

| Artificial-pollen supplementation scope | Original pair rows | Positive seed-log pairs | Mean within-pair ln(Hand/Natural) | Independent-pair bootstrap 95% | Published mean |
|---|---:|---:|---:|---:|---:|
| Whole plant | 32 | 32 | **0.129219** | [0.011979,0.249027] | .129 |
| Inflorescence | 32 | 24 | **0.453325** | [0.162736,0.731887] | .453 |
| Flower | 33 | 26 | **0.787410** | [0.430806,1.168888] | .787 |

The difference between flower- and whole-plant means of log PL is **+0.658191** with exploratory between-group bootstrap 95% **[+0.270222,+1.051231]**. The ratio of log indices .787/.129≈6.10 is not six times **absolute seeds**. Crucially, whole-plant pairs are **size-rank matched but genetically different plants** (0/32 match plant IDs), while flower-level comparisons use the **same** plant (33/33 ID matches). Inflorescence pairs match ID in 31/32 cases. Thus between-scope difference is compatible with resource reallocation **but not causal within-genet identification** of the amount shifted.

**Aggregation trap:** mean of individual-pair ln(Hand/Natural) ≠ ln(mean Hand/mean Natural). Flower-level mean paired log PL=.7874, whereas log ratio of the two group arithmetic means=.5981. Study and meta-analysis estimates must preserve the original estimand, sampling exclusions, and manipulation scale. The flower-level log *fruit-set* contrast in matched positive pairs is **−0.0111**, despite strong seed/fruit contrast, confirming why fruit presence and seed quantity are distinct outcomes.

## Answer to original research question and next falsification

Plant resources or many flowers can fail to convert into **reproduction** at at least three distinct levels, all empirically exemplified but not linked into one universal causal chain:

1. **After flower production, before mature fruit**: even hand-cross-pollinated late Brassica flowers show markedly lower fruit formation. Not a lack-of-visiting-insects explanation.
2. **Inside fruits, between pollination and ovule fertilization**: more than one-third of ovules in sampled developing Brassica fruit were unfertilized despite pollen supplementation; >95% of fertilized ovules developed, though these are conditional on developing fruits.
3. **Between local flower-level and whole-genet seed benefit**: the magnitude of pollen limitation in Veratrum depends strongly on whether pollen is added to a few flowers or all flowers, consistent with fruit/seed resource redistribution among modules. Local flower seed gains are not automatically additional total viable offspring.

**Not yet answered:** how many extra *milligrams of floral construction/maintenance biomass* cause each transition, how parentage success behaves for the same male genets, and how present success trades against survival/future-year reproduction. Vegetative mass cannot stand in for flower biomass. Follow-up priority is original *Silene latifolia* individual phenotypic, repeated flowering, and 2007 survival data, while controlling the selected-line and physiological confounding and not treating repeated counts of open flowers as unique flower construction.

Scripts and source checks are durable in Draft PR #3; public ZIP contents are processed in-memory and **no source plant/site IDs are uploaded**. main branch remains unchanged.
