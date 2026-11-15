# Where floral investment stops converting into reproduction: stage contrast v0.3
Date: 2026-11-15 JST
Status: **new observational paired bootstrap completed on raw original Zurich source; independently published Brassica and Mimulus experiments re-audited; no new onishibari functional fertility observations.**

## Exact inference target
Do extra flowers and flower dry mass increase (i) effective encounter, (ii) compatible fertilization, (iii) a flower's probability of at least one seed, (iv) number of completed seeds per available ovule, (v) paternal genetic success and (vi) viable lifetime recruits? These responses need **separate denominators** and same-genet experimental changes. NDVI, vegetative mass and recorded insect landings are NOT proxies for measured floral allocation or seed-mediated fitness.

### 1. Independent source experiments: distinct conversion points, not one universal coefficient
**Harder & Hamnett-Mount 2026**, American Journal of Botany, doi:10.1002/ajb2.70204, study source https://pmc.ncbi.nlm.nih.gov/articles/PMC13206284/ and Dryad data description https://doi.org/10.5061/dryad.hx3ffbgtq .

- Experiment 1: 86 Brassica rapa plants assigned 25, 50, 75, or 100% of their flowers to hand outcrossing. Fruit formation fails late even in deliberately pollinated flowers. The reported mean temporal fruiting span is **9 days on poorly pollinated plants vs 5 on fully pollinated plants**, despite plants flowering for mean 11.6 ± 4.9 days. **This is an already-published causal experiment on pollination fraction, not a new analysis.**
- Exp 1 fruit-set relationship with ln vegetative mass b=-0.321 ± 0.148 SE; doubling vegetative mass leads to a logit odds ratio ~0.8005 (approx Wald CI 0.6547–0.9788) for **proportional** fruit set. Flower number and fruit number both increase with biomass but at different rates; does NOT imply absolute fruit number decreases, a randomized biomass treatment or real floral dry mass.
- Experiment 2: 46 plants in 25%/75% pollination groups; 719 developing fruit-bearing flowers; mean fraction of ovules fertilized **0.610 [0.567,0.658]**. Among **fertilized ovules**, the proportion developing seeds **0.963 [0.954,0.969]**. Simple arithmetic product 0.610×0.963=**0.5874** is a conditional illustrative seed/ovule ratio **among sampled developing fruits**, **not** seed yield across all initiated flowers or a joint-interval estimate. Mechanisms for 39% failure to fertilize were not separated into pollen incompatibility, pollen tube deficiency or ovule infertility.
- Authors discuss both development resource capacity and hormonal carpic dominance; direct mechanism unresolved.

**Mimulus ringens** randomized displays (2,4,8,16 flowers), paternity assigned for 1,935 progeny, doi:10.1093/aob/mcr193, full text https://pmc.ncbi.nlm.nih.gov/articles/PMC3278290/ .
- Original fitted male total siring elasticity with display **0.46**, female total seeds mothered elasticity **0.98** (model- and population-specific). An 8-fold display increase predicts respectively **2.60× total siring** and **7.67× total maternal seeds** from the published power fits. This corresponds to ~0.325× male siring/flower and ~0.959× seeds mothered/flower for the same display contrast. Effects are derived from published fits, NOT independent data fitting; dry floral mass was not manipulated. Plant total male and female fitness can both increase despite lower marginal returns.
- Published study also documented increased geitonogamy and pollen discounting; **male paternal success** cannot be substituted by visit count or maternal seed success.

### 2. New verified Zurich source reanalysis: same-garden paired difference
Original source: Reji Chacko et al. (2026) Journal of Applied Ecology doi:10.1111/1365-2664.70384; authors' original EnviDat doi:10.16904/envidat.676 and public https://github.com/mrejichacko/BetterBlooms . This project previously linked actual flower counts × minutes for 285 valid garden×plant-type×field-day entries and aggregated original study's end-season reproduction outcomes.

NEW script `scripts/zurich2016_paired_stage_contrast.py` joins **the identical 23 eligible gardens** across two comfrey reproductive outcomes:
- `any seeded flower fraction` = source n_flowers_with_seeds / (n_flowers_with_seeds + n_flowers_without_seeds)
- `ovule seed completion fraction` = source n_seeds / (n_seeds + n_unfertilized_ovules)
- x = source author's Bombus landings divided either by observer hours or approximate open flower-unit × observation hours (the latter counts flower number on days, not tracked exact flower lifetime).
- Unadjusted Spearman correlation, garden-label permutation p (10,000), paired garden bootstrap (10,000) of **difference in correlations**, leave-one-garden-out rho range, and rank partial correlations controlling original Urban_500 and flowering species richness.

| Exposure / statistic | `any seeded flower` | `ovule completion` | Paired difference |
|---|---:|---:|---:|
| Visits per **observation hour** rho | **0.5929** | 0.1591 | **0.4338** |
| Permutation two-sided p (10,000) | **0.00390** | 0.46915 | not a p-value |
| Paired bootstrap 95% of raw rho difference | - | - | **[0.0489, 0.8246]** |
| Leave-one-garden-out rho range | [0.5370,0.6601] | [0.0390,0.2309] | not assessed |
| Visits per **approx flower-unit hour** rho | **0.5375** | 0.2875 | 0.2500 |
| Permutation two-sided p | **0.00780** | 0.18168 | not a p-value |
| Paired bootstrap 95% of per-unit-hour rho difference | - | - | **[-0.0962, 0.5918]** |
| Urban + source floral-richness partial-rank rho, visits/hour | 0.4797 | 0.0274 | 0.4523, bootstrap 95 **[-0.0248, 0.9465]** |
| Urban + floral-richness partial-rank rho, visits/unit-hour | 0.3953 | 0.1135 | 0.2817, bootstrap 95 **[-0.2613, 0.7054]** |

**Meaning:** A reasonably stable garden-level association exists between Bombus observation frequency and the probability a comfrey flower bears ≥1 seed. Yet the **claim that an ovule-completion bottleneck has been located is NOT robust**: strength-of-association differences become compatible with zero after floral-unit normalization and urban/floral-richness adjustment. Correlation across gardens between flower count and outcomes is weak (flower endpoint -0.1571, ovule endpoint -0.1176). These are not flower biomass effects.

**Crucial qualification:** the two garden-level fractions arise from **different sampling/aggregation of flowers and ovules**, even though the gardens match. Recorded landing counts need not be compatible pollen, and no actual flower→ovule transition was tracked on identical flowers. Neither a significant vs nonsignificant pair of correlations nor a bootstrap difference proves a developmental bottleneck or causal mediation. The original authors already studied these data with guild-specific hierarchical models.

**Reproducibility:** GitHub Actions run https://github.com/zuizui0223/onishibari/actions/runs/43712775790 , **SUCCESS**, artifact `zurich-paired-seeded-flower-versus-ovule-completion` ID 14161775703, artifact ZIP SHA256 `a5d96fa43bfa72883d3973d0ec920b373b9cbcfef48e81b28598d317e1f7778f`. Source EnviDat ZIP SHA256 `b1ea0dc54d0c6e33bfacb44d7285377ec3e113807cc945fc2dfbfcc3e534bcc9`. Committed script and workflow `.github/workflows/angiosperm-zurich-paired-stage.yml`; log and JSON include full output. Raw source site IDs and coordinates are not released.

### 3. What is now supported vs unproven
**Supported:**
1. In previously published manipulations, floral-display size can differently affect maternal and paternal marginal success (*Mimulus*).
2. Increasing imposed hand-pollination can shift whole-plant fruit-set constraints and flowering-order success without making all ovules fertilize (*Brassica*).
3. Source-reanalyzed landings are relatively informative about a flower having ≥1 seed in some garden/community contexts, but less consistently about the number of seeds per potential ovule (Zurich comfrey); the gap is confounded and normalization-sensitive.

**NOT supported yet:** stage-switch boundaries estimated for wild long-lived Daphne; floral **construction dry-mass** conversion; direct pollen deposition, pollen compatibility, sex morph-specific genetic siring, seedlings, survival and lifetime performance. Differences in published experimental species and endpoints **cannot** be multiplied to make one global conversion rate.

### 4. Next tight experiments / analytics
- Primary target: repeated `genet × flowering_day × sex_morph × flower` records on two-or-more seasons, same plants tracking floral dry mass, flower count, floral duration, effective stigma compatible pollen, embryo fate, viable offspring and paternal assignments. Study size/power must account for within-genet correlation, not count flowers as independent plants.
- Randomized **visual display manipulation** with sham and unaltered receptive flowers to isolate advertisement, alongside a separately randomized **flower thinning / fruit removal** treatment to isolate resource/fruit-order feedback. Note thinning changes both advertisements and carbon budget; the two operations are not interchangeable.
- Split *Daphne pseudomezereum* morphs by independently verified functionality. The current focal target has **no own experimental pollen or seed function data**; retain the existing branch's protective max-flower manipulation and positive controls. No location-level private records in public artifacts.
- A further comparative test should ask whether **the limiting transition switches at a different flower number/flower-mass threshold depending on visit reliability and paternal/female mating pathway**. Confirm across independent clades, not by mixing garden/pollinator aggregate correlations into an artificial meta-effect.

### Raw data access note
The 2026 Brassica Dryad landing page lists eight original files including individual flower and ovule fate CSV. This environment could inspect file descriptions and primary paper but **could not download the raw Dryad CSV**, whose programmatic download endpoint requires an authenticated Dryad token under current service controls. Therefore no new individual-level Brassica model was fitted. No registration credentials or bypass attempted; published coefficient transformations and original authors' stage conclusions are labeled as such.
