# 2026-11-15: Direct biomass, same-season NDVI, measured visitor scarcity — evidence upgrade

**Research question:** Do plant mating/sex strategies, biomass and pollinator encounter rates track ambient vegetation production and environmental fluctuations? Is the apparent mismatch in Onishibari a case of evolutionary decoupling?

**Status:** empirical two-new-source audit and source-matched remote sensing are complete, but **no causal evolutionary model or individual flower-pollen-deposition data have been fitted**. Keep the broad PR #3 separate from Daphne P/C/K phylogenomic PR #2.

## A. Direct measured individual whole-plant biomass: BAAD × Wang 2021

Source: Falster et al. BAAD v1.0.0 https://github.com/dfalster/baad/releases/tag/v1.0.0 ; Wang et al. (2021) https://doi.org/10.5061/dryad.s1rn8pk7n.
BAAD `m.to` **is total dry biomass in kilograms (kg)**, `m.rt+m.so`, not remotely sensed AGB nor dry mass of reproductive tissues. `h.t` is total height in metres (m). The initially generated output mislabeled mass CSV column `*_mass_g`; **all stored numeric values were kilograms, not grams**. This is corrected in source script `scripts/baad_wang_biomass_pilot.py` and subsequent run; initial source audit run 43696999947.

- Source BAAD: 21,084 individual observations, 656 distinct binomial names recognized after normalization (678 species reported by source using its own taxonomy).
- BAAD 8,878 individuals have valid `m.to`.
- Of those, **4,562 individuals, 115 species** match a Wang sexual system; **4,133 individuals, 91 species** have both height+dry biomass.
- 91 species: hermaphroditism 60 spp/1,984 individuals, dioecy 12 spp/347 individuals, monoecy 19 spp/1,802 individuals. Species-sampled median `m.to`: hermaphroditism 0.00541 kg, dioecy 0.00107 kg, monoecy 2.945 kg.
- Allometric pooled `log(m.to) = -2.90366 + 2.06642 log(h.t)` from 4,133 individual observations, *not* phylogenetic or study-balanced, and not appropriate as universal allometric law.
- **Only 6 plant families** among the 91 spp contain more than one sexual-system class. A family-level contrast is severely underidentified: Aceraceae monoecy 2 spp vs dioecy 1; Myrsinaceae 1 hermaphrodite vs 1 dioecy; Theaceae 1 vs 1; Oleaceae 2 dioecy vs 1 hermaphrodite; Malvaceae 1 hermaphrodite vs 1 monoecy; Lauraceae 1 dioecy vs 1 hermaphrodite.
- Strongly different sampled stage/height distributions and biases in woody taxon/age/study composition mean group medians **do not establish reproductive strategy → individual biomass**. BAAD lacks floral allocation / lifetime fitness. Do not call `m.to` a sex-specific floral investment ratio.

### Decisions
Preserve all 4,133 actual mass+height observations grouped by species/study (not treat individuals as independent evolving strategies). For a causal allometry test require repeated sex-morph male/female plant dry masses within the same populations, body size and age, plus reproductive-tissue biomass. Independent phylogenetic replicates across sex-state transitions are especially lacking in this BAAD×Wang overlap.

## B. Actual 2016 vegetation satellite data × 24-garden pollinator observations

Source: Reji Chacko et al. (2025) EnviDat 10.16904/envidat.676; measured actual visits in **2016 Jun 15–Jul 20** 24 residential Zurich gardens, 4 standard phytometer plants; NASA MODIS/061/MOD13Q1 (Terra, 250m 16-day, NDVI and EVI, 0.0001 scale) accessed via Microsoft Planetary Computer STAC. Source garden coordinates remain in private runtime and were **NOT exported** to public GitHub artifacts.

A first correct-source 2016 run included two composites that did not overlap the field season. The final validation restricted MODIS to **three** 16-day composites starting 2016 Jun 09, Jun 25 and Jul 11 (the nominal composite windows overlap the field season), all 24 gardens passed `pixel_reliability` QA `0/1` and all were **24 distinct MODIS 250m pixels** for each of these scenes. This is contemporaneous ecological exposure, not a lifetime environmental trajectory or pollen deposition measurement.

Source final run https://github.com/zuizui0223/onishibari/actions/runs/43697110784; artifact `angiosperm-zurich2016-time-matched-MODIS-visitation`, ID 14154643154, ZIP digest `sha256:62bfb3e313e8b4a31f57403dd712aa95fbfee2273bc4692941555997031e9e76`.

**24-garden site-level Spearman correlations between 3-scene mean NDVI and mean insect visits per original observation window:**

| Plant (common/source name) | NDVI vs visitors rho | unadjusted p | q BH among 4 plants × 2 VI tests |
| --- | ---: | ---: | ---: |
| carrot / Daucus carota | +0.46161 | 0.02317 | 0.07248 |
| sainfoin / Onobrychis viciifolia | +0.49486 | 0.01395 | 0.07248 |
| radish / Raphanus sativus | +0.38547 | 0.06285 | 0.10057 |
| comfrey / Symphytum officinale | −0.16129 | 0.45149 | 0.45149 |

EVI correlations in the same 24 sites: carrot +0.4342 p=.0340 q=.07248; sainfoin +0.42944 p=.03624 q=.07248; radish +0.31020 p=.14015 q=.18687; comfrey −0.19268 p=.36704 q=.41948. Thus **no BH q < .05 across the eight comparisons**. Nominally positive associations for two or three plants need replication / source controls. NDVI is not "pollinator abundance" or direct flower count. All eight tests share 24 gardens (nonindependent); geographical autocorrelation, urbanization, specific flower resources and window effort remain unaccounted for.

Previous experiment data output: 21,120 observation windows, **94.31% with no visit** (~19,919 zero windows), **5,794 visitors in nonzero windows**. Site-level visits vs reproductive outcomes: carrot seed/umbel rho −0.100, radish fruit-set rho +0.206, sainfoin fruit-set +0.194, comfrey seeded-flower fraction rho +0.0458; none nominally p<.05 (N=24 gardens each). Source successful run https://github.com/zuizui0223/onishibari/actions/runs/43696492759 . Reproductive responses differ in units: cannot combine 4 species in a common unrestricted regression without standardizing/contrasts. We have NOT measured/estimated effective compatible pollen deposition or causal mediation by visitors.

### Conditional research hypothesis
Landscape greenness and its temporal covariation may be upstream context for partner encounters, while pollinator efficiency, compatibility, mate access and resource limitation gate effects on seed set. **NDVI → visit counts → realized reproduction is a hypothesis, not demonstrated mediation**. Specific plant-pollinator guild relationships may differ; zero visitor windows cannot alone identify effective pollination scarcity.

## C. Environmental R^2 and geographic prediction (a negative test)

Source GloPL×Wang×WorldClim v2.1 (1970-2000 10' rainfall and elevation), 580 species. Under pre-specified 5-fold geographic block (20 degree) and plant-family group holdouts, ridge regression (alpha10) using latitude/elevation baseline plus average rainfall, monthly rainfall CV and sexual system gave:
- spatial held-out MAE baseline 0.5834 versus full model 0.5817 (**+0.29% improvement**)
- family held-out MAE baseline 0.5863 versus full 0.5898 (**worse by ~0.60%**)
- all held-out R^2 were negative.
- Source validated run https://github.com/zuizui0223/onishibari/actions/runs/43695454189 .
Hence no reproducible broad **static rainfall+sex system** predictive signal for GloPL experimental pollen limitation in the tested holdouts. This does not negate flowering-season rainfall anomalies, NDVI phenology, forest canopy, effective pollinator guilds, or individual investment.

## Next model gates and remote sensing variables

**Layer 1 quantity**: true individual floral/vegetative dry mass; site MOD17 GPP annual carbon fixation; woody canopy AGB GEDI/ESA CCI with uncertainty. Note GEDI 2019+ is temporally incompatible with many historic GloPL studies; woody AGB not focal understory plant biomass.

**Layer 2 temporal reliability**: MOD13Q1 NDVI/EVI within-year amplitude and anomalies, CHIRPS daily rainfall anomalies/dry streaks, actual flowering/pollinator season, observation-year matched; controls for clouds and QA. Year of GloPL publication ≠ year of sampling; avoid treating 2013 snapshot as contemporaneous for most sites.

**Layer 3 accessibility**: available compatible mate density, nocturnal vs diurnal effective pollen transfer, garden local floral resource cover, open canopy/light, guild diversity, predator/nectar robber effects. Detectability-corrected visits per flower-hour and effective pollen deposition need distinct response models, even when using public observations.

**Layer 4 long-run outcomes**: fruit/seed quality, annual reproductive output, adult survival, paternity-mediated gene flow, changes in mating system/sexual dimorphism after lineage splits. None of these implied solely by RS associations.

Next important validation is to obtain garden-level flowering resources and urban landcover (or remote tree cover) **for the same 2016 season**, and separate visitor guilds from indiscriminate insects; a landscape NDVI effect may disappear after conditioning on local resources. Publish only anonymized garden-level summaries and model coefficients, not personal house/garden coordinates. For woody sex vs biomass, seek direct reproductive-tissue mass data in matched sexes, not solely total plant mass across lineages.

## Repository and source paths
- `scripts/zurich2016_modis_visits.py` (field-year matched MODIS and BH+pixel QC)
- `scripts/visitor_scarcity_real_data_pilot.py` (measured visits vs seed/fruit outcome)
- `scripts/baad_wang_biomass_pilot.py` (direct kg dry biomass and allometric residuals)
- `scripts/glopl_spatial_cv.py` (explicit out-of-region holdout)
- `docs/remote_sensing_pollination_biomass_protocol.md` (scaling and source policy)
- These are reproducible exploratory analyses on Draft PR #3, not peer-reviewed conclusions.

## Within-family sex-system comparability check (corrected source run 43697541650)
Only **six named families** in the 91 matched biomass+height species contain at least two sex-system categories, each with one to three sampled species in this data slice:
- Aceraceae: two monoecious vs one dioecious species
- Myrsinaceae: one hermaphroditic vs one dioecious
- Theaceae: one hermaphroditic vs one dioecious
- Oleaceae: two dioecious vs one hermaphroditic
- Malvaceae: one hermaphroditic vs one monoecious
- Lauraceae: one dioecious vs one hermaphroditic
The sign of **height-adjusted** dioecious vs hermaphroditic mass residuals is mixed across the four applicable families (Myrsinaceae & Theaceae positive for dioecy; Oleaceae & Lauraceae negative for dioecy). These are sparse, unmatched species-level contrasts, not a robust phylogenetic analysis. Therefore this dataset does **not** support a universal sex-system dry-biomass difference after conditioning on size and lineage. Experimental flower tissue dry-mass comparisons from within populations remain the key missing variable.
