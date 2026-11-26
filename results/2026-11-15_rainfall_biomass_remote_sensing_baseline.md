# Do plant reproductive strategies covary with biomass, remote sensing productivity and rainfall? First evidence audit

2026-11-15. Question inspired by sex dimorphism and intermittent pollination in summer-deciduous woody *Daphne pseudomezereum*. This is a **broad comparative ecological research program**, NOT P/C/K introgression research.

## Clear definition of outcomes and causal scales
- **Individual dry biomass** (g/plant) and **floral biomass allocation** (g floral tissue/g plant) are direct functional investment traits. Not measured by either PROTEUS height or remote-sensing NDVI.
- **Landscape woody AGB** (Mg ha^-1; ESA CCI Biomass or NASA GEDI), **GPP** (kg C m^-2 year^-1; MOD17), **NDVI/EVI**, **canopy LAI** are distinct ecosystem variables. High woody canopy AGB can coincide with low understory light. NDVI may saturate; vegetation greenness is not flower resource availability or pollinator abundance.
- **RF** provisionally means rainfall (if Random Forest was intended, that is a modeling algorithm, not an RS covariate): use CHIRPS v3 for satellite/gauge rain (60S–60N) when actual observation dates are known; WorldClim static normals for first climate background.
- **True partner scarcity** must be decomposed into visits per flower-hour (with detection effort, pollinator guild and time of day), effective compatible pollen deposition, hand cross-pollination response (pollen limitation), compatible mate abundance, and temporal reliability/cross-year variance. They are not equivalent.
- **Evolutionary strategy vs current environment** cannot be causal without ancestral trait reconstruction, geographic/phylogenetic controls, or manipulative/panel evidence; do not claim that contemporary NDVI selected a deep ancestral sex system.

## Actual analysis A: world precipitation vs pollen limitation and sex systems
Source GloPL published ES from https://github.com/idiv-biodiversity/pollen-limitation-data-descriptor/blob/master/Output/GloPL.csv (2,969 cases), Wang et al. 2021 deposited sexual system list https://doi.org/10.5061/dryad.s1rn8pk7n, and WorldClim v2.1 10 arcminute 1970–2000 precipitation normals/elevation from https://geodata.ucdavis.edu/climate/worldclim/2_1/base/.
- Source-matched selected/numeric sex system cases **1,583** across **611 species**, at **736 exact distinct coordinate pairs**. These are not 1,583 independently evolved species.
- Raster has land values for **580 species**. 45 experimental cases fell in no-data precipitation raster cells.
- **Unadjusted species-level group median climatological rainfall** (mm/year): hermaphroditism n=503 => **974**; dioecy n=34 => **1,202**; monoecy n=43 => **883**.
- **Monthly rainfall CV** medians: hermaphroditism 0.5158; dioecy 0.4296; monoecy 0.5915.
- **Spearman rho** against species-level median experimental PL effect:
  - Annual RF ρ=+0.04071 (naive p=0.328)
  - Mean monthly RF CV ρ=+0.08403 (naive p=0.0431)
  - Elevation ρ=-0.04531 (naive p=0.276)
  In dioecious subset n=34, annual RF rho=-0.160, RF CV rho=+0.107; both unadjusted p>0.3. These correlations are small, and nominal p=0.043 for seasonality is not evidence after multiplicity, non-independence, and unmodeled phylogeny/habitat.
- **Critical error guard:** GloPL Year is PUBLICATION year, not necessarily the field observation year; 2001–2015 MODIS is **not contemporaneous exposure** for all historic GloPL studies, and WorldClim is a climate norm. Do not label this "remote sensing NDVI test". WorldClim 10' pixels are coarse, and site coordinates may be approximate.
- Reproducible source/quantitative work: https://github.com/zuizui0223/onishibari/actions/runs/43693752526 ; artifact `angiosperm-GloPL-geographic-rainfall-climatology-pilot` ID 14154185500, SHA256 `b1233499c898115f6123ef90c8a341d0f874a14a5d347090da9d84f9b109faf3`.
- Artifact provides `species_median_PL_worldclim_10m.csv`, `site_coordinates_for_ndvi_sampling.csv` with **736** sites, and `pilot_summary.json`. The first failed pipeline run was corrected for zero-padded month filenames and later successful runs verified.

## Actual analysis B: independent-plant size proxies, NOT biomass
PROTEUS/DiveRS original authors' cleaned public trait data: https://github.com/divers-it/rs-traitspace/blob/main/supplementary_data/cleaned_data.csv, 364 angiosperms. Height known 342, flowerSize known 323, both 305.
- Rank correlation of maximum height with flower size among 305 species: **rho = -0.0019** (p=0.974).
- Woody subset n=142: **rho = -0.0355** (p=0.675).
- Herbaceous subset n=158: **rho = +0.2151** (naive p=0.00664). These are **unadjusted** cross-taxon descriptions; five species are woody/herb categorization missing/other with both traits.
- Selfing 40 species: median maximum height=0.3775 in source units, woody **5/40 (12.5%)**. Outcrossing 195 species: median height=5.0, woody **127/195 (65.1%)**. This difference is strongly entangled with growth form, taxonomic composition, habitat and data selection; it is **not biomass measurement and does not prove a selection process or relative efficiency**.
- Morphologically dimorphic sexual system is not verified functional dioecy; analyze sex expression separately from male/female actual flower production and fitness.
- Successful quantitative run: https://github.com/zuizui0223/onishibari/actions/runs/43694073229 ; artifact ID 14154235207, SHA256 `91c6699fba4d6eaa6f0c9420b8a0b05c01552d38bbb4481138c7978930d68759`.

## NDVI/EVI/GPP and true woody biomass: extraction status & design
Earth Engine code editor source template: `scripts/gee_glopl_modis_context.js`.
Input: **736 site coordinates** extracted from the successful GloPL/WorldClim real-data run.
- Source: `MODIS/061/MOD13Q1`, 250m/16 days, QA mask `SummaryQA<=1`, scale by 0.0001, 2001–2015 monthly and annual/overall context NDVI/EVI.
- Source: `MODIS/061/MOD17A2HGF`, GPP 500m/8 days; sum by calendar year and then average across 2001–2015. Do not call GPP aboveground biomass.
- Buffers: 250 m, 1 km, 5 km as exploratory spatial-sensitivity bands, not universal pollinator foraging distances.
- **NOT RUN:** Earth Engine credentials and uploaded site asset are not connected here, and source study dates cannot be inferred from publication year. The source script contains a placeholder for an EE asset and must not be reported as a measured NDVI association.
- Biomass: GEDI L4B version 2.1 maps woody aboveground biomass density and uncertainties for 2019–2023 at 1km; historical GloPL pollen experiments 1981–2015 are temporally mismatched. ESA Biomass CCI has annual woody AGB and uncertainty at 100m for multiple years (2005–2012, 2015–2024); only link to date-verified study years and forest-dominant habitat; does not measure understory Daphne biomass.
- Long NDVI alternative: NOAA AVHRR NDVI CDR daily 0.05-degree ~1981–present, useful for broad vegetation phenology if actual field year and geographic sampling uncertainty are known.

## The ecological prediction worth prioritizing
**Quantity–reliability–accessibility decomposition**:
- Quantity: annual GPP, standing AGB, mean NDVI.
- Reliability: within-year NDVI amplitude, interannual NDVI anomalies, consecutive drought years, rainfall seasonality and RF autocorrelation.
- Accessibility: canopy closure, understory light, landscape heterogeneity, flowering timing and effective available partner guilds.
- Mechanistic endpoints: effective pollen deposition, hand-supplementation PL, viable seed set, male pollen export and long-distance pollen paternity, individual dry-mass allocation and adult survival.
- Response of sex system depends on ancestry, habit, mating compatibility, clonality, adult longevity. A dense forest may exhibit **high canopy biomass yet low understory pollinator accessibility**. Do not assume an always-positive effect.

Use hierarchical models with study/year/sex morph/population nested properly, phylogeny, spatial dependence, location-precision checks, outcome type and experimental level. A test of a single NDVI/AGB correlation does not demonstrate divergence, maladaptation, or pollinator scarcity.

## Hypothesis falsification and design gates
H1: Sex dimorphism reflects **individual size and floral investment allometry**, not necessarily total landscape biomass: compare flower dry mass ratio and plant dry biomass within sex × population, common garden when feasible.
H2: Effective pollination is better predicted by **flower-season** vegetation and rainfall reliability + pollinator phenology than mean NDVI/annual RF. Counter: lag correlations disappear when actual flower and pollinator observations are supplied.
H3: At fixed landscape GPP, canopy AGB interacts with understory light/accessibility. Counter: no sex-specific or guild-specific interaction.
H4: Rare effective cross-pollen transfer may maintain genetic connectivity but not high annual seed output: compare paternity and seed set from same marked plants.
H5: Apparent richness deficit from GBIF may reflect detection; validate hierarchical occurrence and timed visits.

## Repositories
Broad workstream PR #3 (research/angiosperm-reproductive-decoupling-20261115) is separate from P/C/K cytonuclear PR #2; main remains unchanged. This document is analysis evidence and a preregistration-style direction, not a published finding.
