# Remote-sensing test of reproductive-strategy, biomass, rainfall and pollination

Status (2026-11-15): **initial static precipitation pilot run submitted to GitHub Actions; NDVI/GPP/GEDI extraction NOT executed**. It requires an authenticated Earth Engine project, proper site geometry import and validation of source observation dates.

## Question
Does reproductive mode (dioecy, monoecy, hermaphroditism; autonomous selfing; sex-function allocation) covary with (a) ambient **vegetation production**, (b) **local plant biomass**, (c) **resource seasonality**, and (d) **pollinator availability and effective pollen delivery**? Are any patterns shared across plant lineages or conditioned on canopy, climate and individual life history?

**RF** is provisionally interpreted as **rainfall**, not Random Forest (machine-learning algorithm). WorldClim 2.1 rain is meteorological interpolation; CHIRPS v3 is gauge + satellite precipitation. NDVI ≠ direct biomass, NDVI ≠ floral abundance, NDVI ≠ pollinator numbers.

## Distinguish scales of biomass, evidence and units
1. **Individual biomass**: dry mass of plant vegetative organs, male/female flowers, fruit and seeds (g or g/individual). This tests allometric sex allocation, and requires field data / trusted TRY or BIEN size models. Remote sensing is inadequate for understory *Daphne* individual biomass.
2. **Vegetation productivity**: MOD17 GPP (kg C m^-2 yr^-1), NIRv/spectral indices, vegetation greenness. May influence resource competition and pollinator activity, *not identical to standing plant mass*.
3. **Landscape above-ground woody biomass**: GEDI L4B (Mg ha^-1) with pixel-level uncertainty (1km, observed mainly 2019-2023) or ESA CCI Biomass (~100m woody AGB, selected years 2005-2024 with documented uncertainty). Dense canopy can increase AGB while decreasing understory irradiance; forest AGB is not understory focal plant biomass.

## Satellite and meteorological candidate data
| Variable | Source | Native grain | Valid time | Primary role | Limitation |
|---|---|---|---|---|---|
| NDVI/EVI, annual and monthly amplitude | MOD13Q1 v6.1 | 250m, 16-day | 2000–present | landscape greenness and seasonality | saturated at high LAI, not floral biomass |
| LAI/FPAR | MCD15A3H v6.1 | 500m, 4-day | 2002–present | canopy growth / shade proxy | modelled |
| GPP | MOD17A2HGF v6.1 | 500m, 8-day | 2000–present | carbon fixation/productivity | modelled |
| AGB mean/SE | GEDI L4B v2.1 | 1km, one-period summary | 2019–2023 | woody canopy aboveground biomass | 2018-and-earlier GloPL experiments temporally mismatched |
| woody AGB/uncertainty | ESA CCI Biomass | 100m + coarser | 2005–2012,2015–2024 | true woody AGB context | does not measure nonwoody understory biomass |
| RF precipitation normals | WorldClim 2.1 | 10' in fast pilot | 1970–2000 baseline | broad climate context | static climatology, not experiment-year rainfall |
| RF precipitation time series | CHIRPS v3 | 0.05° daily, ~1981–present | year-matched only with true experiment year | rainfall anomalies, dry spells | quasi-global 60°S–60°N; no direct biomass |

Official sources:
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD17A2HGF
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD15A3H
- https://doi.org/10.3334/ORNLDAAC/2299
- https://climate.esa.int/en/projects/biomass/
- https://www.chc.ucsb.edu/data/chirps3
- https://geodata.ucdavis.edu/climate/worldclim/2_1/base/

## Important temporal caveat for the existing GloPL corpus
GloPL cases are experiments published between 1981–2015. The CSV `Year` field is **publication year**, not the actual field-observation year for every case. **Never match 2015 MODIS imagery to a site solely because its paper was published in 2015.** The GloPL georeferenced data can be assigned long-term baseline climate (WorldClim 1970–2000) and broad 2001–2015 vegetation climatology as a **retrospective habitat proxy**, not evidence of contemporaneous environmental exposure.

Historical NDVI time series back to 1981 can be obtained from NOAA AVHRR NDVI CDR at 0.05°, but precise matching still needs a verified study observation date. Official data: https://www.ncei.noaa.gov/products/climate-data-records/normalized-difference-vegetation-index .

## Steps now
1. `scripts/reproductive_worldclim_rainfall_pilot.py` downloads 10-minute monthly rainfall and elevation rasters from WorldClim, extracts values for matching GloPL lat/lon, and estimates annual rainfall normal, monthly precipitation CV and elevation. It exports `site_coordinates_for_ndvi_sampling.csv` keyed by a hash of reported site coordinates, **not by an inferred field year**.
2. Upload that CSV to Earth Engine as a point table (longitude and latitude columns), set your valid Earth Engine project asset ID in `scripts/gee_glopl_modis_context.js`, run three buffer summaries (250m/1km/5km). In 2001–2015 climatological composite, use QA-masked MOD13Q1 NDVI/EVI, seasonal amplitudes of monthly climatological VI, annual-mean MOD17 GPP. This script is a **ready-to-execute template, not a claimed completed export**.
3. Only after site time matching, evaluate climate/NDVI anomalies, flowering DOY overlap, and seasonal resource gaps. GloPL alone is insufficient for temporal alignment.
4. Keep the separate canopy biomass products for newer, time-compatible field sites. For *Daphne* understory, validate canopy openness, understory light and floral display with on-site measurements.

## True pollinator scarcity as a measurement problem
- **Encounter rate:** species-corrected visits per flower-hour, time of day, pollinator taxon, temperature and weather, observer effort; use zeros and hierarchical detectability.
- **Functional pollen transfer:** conspecific pollen deposition or pollen tubes per visit, standardized against compatible-mate availability. A visit is not equivalent to pollination.
- **Pollen limitation:** additive hand cross-pollination effects on seed production, separated from bagged-self, compatibility and inbreeding depression; GloPL is a compilation of experimental outcomes, NOT pollinator abundance.
- **Temporal reliability:** coefficient of variation and autocorrelation of effective visit rate, frequency of consecutive bad years/seasons, overlap of focal flower phenology with pollinator activity, alternative guild compensation (diurnal/nocturnal).
- **Landscape capacity**: pollinator richness/occupancy models corrected for record effort plus habitat diversity; GBIF pollinator occurrence is not abundance.

## Pre-registered competing predictions
H1 **Resources-to-biomass**: Higher GPP/AGB corresponds to greater landscape resource production, but within dense forest understory pollen transfer may decline with canopy closure despite high AGB. Test forest habitat × light × flower abundance rather than NDVI-only slope.

H2 **Productivity variability versus strategy**: Temporal NDVI amplitude and rainfall anomaly persistence, controlling climate, phylogeny, land cover and sampling, predict prevalence of autonomous assurance, partner-route redundancy, or chronic pollen limitation. The effect can depend on ancestry and lifespan rather than being monotonic.

H3 **Phenological mismatch**: The circular distance between focal flowering peak and *community* vegetation green-up/peak (measured NDVI), plus measured flower availability and pollinator phenology, predicts functional pollen transfer better than annual NDVI alone. NDVI is not flower phenology; ground truth required.

H4 **Biomass allocation**: Male/female floral display allometry relates to individual mass and resource reserves, rather than remote landscape canopy AGB. Test standardized sex morph dry-mass and flower-size measures in common environments.

H5 **Evolutionary consequences**: A low visit rate can coexist with rare long-distance pollen export sustaining gene flow while local seed set remains pollen-limited. Use paternity, population genomics and coalescent models; F_ST or cp/ITS discordance alone does not establish speciation.

## Essential confounders, scales and safeguards
- Source `study × species × experimental case × site × actual year`, with species and study random intercepts, family/phylogeny and spatial correlation; do not falsely treat 2,969 GloPL cases as independent populations.
- Rainfall, temperature, elevation, biomes, latitude, canopy, plant habit (woody/herbaceous), plant size, disturbance, urban agriculture, soil and landcover.
- Satellite sampling buffers 250m, 1km and 5km are sensitivity choices, **not automatic bee flight ranges**. Foraging traits should determine meaningful scales.
- Report GloPL latitude/longitude accuracy; site coordinate inferred from articles can be coarse enough that 250m pixel extraction is misleading.
- Select minimum-coverage or inverse-probability weighting for photographic, field and study selection; avoid full-data-case restriction.
- For seed production effect-size meta analysis use measurement SD and sample size when available, random effects and phylogeny; no claim of causality from cross-sectional NDVI/biomass correlation.

## Verified follow-up: same-year NDVI and published urbanization controls (2026-11-15)

The original top-of-document status ("NDVI not executed") applies only to **the Earth Engine 736-site long-term GloPL composite workflow**. A **separate, genuine MOD13Q1 satellite analysis has now completed** for 2016, same-year 24 Zurich garden sites using Planetary Computer STAC and pixel reliability QA. Do not describe all NDVI extraction as pending.

1. True 2016 MODIS NDVI/EVI: 24/24 gardens, three 250m 16-day scenes; original visitor matrix 21,120 site × plant × time windows and 5,794 observed insect visitors. <https://github.com/zuizui0223/onishibari/actions/runs/43697781055>.
2. Plant–insect guild response: 100% of 5,794 visits mapped to 167 source taxa; carrot hoverfly visits versus NDVI rho=+.659 (q_BH=.0124 over 27 nondegenerate tests), while carrot bee visits rho=+.251. <https://github.com/zuizui0223/onishibari/actions/runs/43701737786>.
3. **Essential negative control:** the same 24 gardens' `Urban_500` impervious cover and local flower richness are already published by Chacko et al. (2026), DOI 10.1111/1365-2664.70384, source <https://github.com/mrejichacko/BetterBlooms>. NDVI was strongly inversely associated with imperviousness (rho=-.927), but not local flowering species richness (rho=-.035). After adjusting for both, NDVI–carrot-hoverfly partial-rank rho fell to -.155. An independent NDVI covariate worsened 24-garden leave-one-out MAE 2.19% relative to imperviousness+local flowers. <https://github.com/zuizui0223/onishibari/actions/runs/43702182669>.
4. Direct woody individual biomass from BAAD v1.0.0 joined to Wang sexual system labels: 4,133 mass+height+sex individuals, 91 species. Only 6 families have mixed sexual-system states, none with 2+ species in both contrasted groups. **Biomass and sexual system causal/evolutionary claims are not established**. <https://github.com/zuizui0223/onishibari/actions/runs/43697541650>.
5. Full evidence: `results/2026-11-15_ndvi_guild_urban_biomass_integrated.md`.

**Updated ecological priority:** distinguish high ecosystem productivity, low sealed-surface fraction, dense but shaded canopy, local flower resources, and effective pollen delivery. The original 2026 Zurich authors already demonstrated site-level urban densification and flower-richness effects: our new remote-sensing association should not be presented as a novel demonstration of those effects. A new broad study must replicate across habitat types and species, with true biomass allocation and more than one year.



## 2026-11-15 completed independent UK PoMS 1km remote sensing result (2023)
The 736-site GloPL GEE historical-climatology script is **still not executed**, but another satellite study **has been executed and validated**. Using UK PoMS 2017–2023 timed surveys (EIDC 10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618), we retrieved 2023 MOD13Q1 NDVI for 101 OSGB/OSIE 1km survey squares from 8 May–August QA-selected MODIS tile-month images (4 subpixel quadrants per square). For 547 2023 recorded FITs in 88 squares with 2+ valid months, NDVI mean and seasonal range added **no consistent predictive improvement** over models already given target flower, weather, local flower units/cover/context, and source landcover/habitat. Independent 1km-square heldout MAE relative change: all insects −0.43%, hoverflies −1.22%, all bees −0.53%, bumblebees +0.17%. This is **not** a failure of remote sensing generally, and NDVI ≠ actual floral biomass.

Successful run: https://github.com/zuizui0223/onishibari/actions/runs/43704350736 . For full provenance, temporal limits and comparisons with Zürich see `results/2026-11-15_ukpoms_flower_access_vs_habitat_replication.md`. Next key gap is *actual flowering dates* and functional pollen-delivery observations at independent sites, not another simple average-NDVI regression.


## New 2026-11-15 completed: flower-count elasticity at fixed effort
UK PoMS observed **10-minute fixed 50cm quadrat landings** scale **sublinearly** with counted target floral units. Structured 1km-square subset 3,048 surveys, 101 square clusters, after plant target, year/month/weather, floral type and square fixed effects:
- All insect visitor counts beta(log F)=**0.169 [cluster-robust 95% CI 0.095,0.244]**; doubling observed floral units predicts total landings ×1.124, equivalent to **recorded insect landings per floral unit ×0.562**.
- All bees beta=.217 [.134,.301], hoverflies .090 [.003,.177], bumblebees .144 [.024,.265].
- Within UK PoMS's **different public-volunteer sampling stream** 13,372 valid observations / 2,978 grid cells, source-harmonized adjusted all-insect beta=.092 [.049,.134]; allocated-square beta=.220 [.159,.281]. Both source streams show sublinear count relationships, but volunteer sample selection differs strongly.
- Any-visitor event in a 10-minute quad increases with F: all-insect adjusted OR **1.223** per doubling; all-bees OR **1.241**, while per-structural-unit total landing counts drop. These are complementary definitions of resource access.
- Respects the known prior *Digitalis purpurea* dilution literature (Grindeland et al. 2005, DOI 10.1111/j.1365-2435.2005.00988.x): **not a new ecological principle**, and cannot be called actual pollen-deposition/seed-set dilution. Floral units differ among flower heads/spikes/umbels/individual flowers.
- Open-versus-closed floral shape interaction appears strong but **not identifiable independent of plant taxon**: only the heterogeneous `Other` source plant group contains both categories. No within-taxon morphological selection conclusion.
- Floral patch isolated/same flower/mixed flower did not significantly modify slope in this analyzed sample (all-insect joint p=.849); not evidence of a universal null.

Completed report: `results/2026-11-15_floral_density_sublinear_visitors.md`; effect table `data/2026-11-15_poms_floral_unit_scaling_effects.csv`.
Reproducible runs:
- main fixed site + sensitivities https://github.com/zuizui0223/onishibari/actions/runs/43705791623
- two PoMS streams https://github.com/zuizui0223/onishibari/actions/runs/43706325276
- zero-visitor GLM https://github.com/zuizui0223/onishibari/actions/runs/43706071915
- hoverfly vs bee composition https://github.com/zuizui0223/onishibari/actions/runs/43705978912
- floral shape/patch interaction and identity audit https://github.com/zuizui0223/onishibari/actions/runs/43706798707

Updated conceptual bottleneck: **biomass investment (individual dry mass) → visible floral display (F) → total visitor attraction (Y; sublinear) → per-flower effective compatible pollen → seed viability and male pollen export → lifetime genetic connectivity**. Each arrow needs different data and cannot be proved from FIT Counts alone. Future direct experiments should vary display *within species*, sex morph and physical habitat over years, measure per-flower pollen deposition and reproductive allocation biomass, and evaluate partner scarcity by observation time+functional delivery rather than occurrence records/NDVI.
