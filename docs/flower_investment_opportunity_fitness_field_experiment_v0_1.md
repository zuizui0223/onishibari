# Reproductive investment → effective encounter → compatible pollen → lifetime fitness
## Multi-clade comparative field experiment v0.1 (2026-11-15)

The overarching question is: **when do flowers invested per unit biomass increase effective pollen deposition and male pollen export, when do they only help escape zero-visitor periods, and when do additional flowers experience diminishing returns?** *Daphne pseudomezereum* and its close relatives are focal cases, but the causal prediction should be replicated across distinct clades and environments.

### Evidence motivating a test, not already demonstrated causal effects
1. UK PoMS 2017–2023 FIT Counts: 3,048 allocated 1km observations + 13,372 volunteer stream show strong between-quadrat **sublinear visitor counts per structural floral unit** (all insects beta ~.17 in allocated-square site-controlled model, beta ~.09 in volunteers). No direct effective pollen or seed data; flower units vary across flower-head and umbel structure.
2. Reji Chacko et al. 2016 Zurich actual site-date flowers and visitors (published 2026 JAE) yield 285 active site×date×species observations: **within garden bloom changes are closer to proportional** in fixed-alpha NB/date adjusted models (beta ~.87–1.13 across four species). Poisson / seasonal adjustments disagree for carrot. Therefore **no universal visitor saturation coefficient**.
3. End-of-season reproductive associations are trait/group-specific: comfrey bumblebee visitation vs fraction flowers setting seeds rho=.538 after approximate bloom-minute normalization; ovule-level seed completion correlation weak. The source authors already studied group-specific reproductive success. This is **not new experimental evidence**.
4. Same garden/day conditional null: rare sainfoin 23 zero-visitor days versus ~19.9 expected under fixed garden total and measured flower×observer-minute exposure (p_excess≈.13). Carrot 7 zeros versus .5 under linear F exposure but 7.54 under exponent F²; **not identifiable whether low-display threshold vs phenology/weather causes it**.
5. Same-garden zero/nonzero label permutation (10,000) with common floral DOY detrending: **sainfoin** visits favored positive-bloom residuals (p_two-sided=.0108); **carrot** apparent F effect disappeared after detrending (p=.178). Supports testing **encounter-occurrence** separately from **visitor intensity** without inventing fitness impacts.

### Experimental causal structure (not merely species trait correlations)
Unit: `genet × population × date/period × sex morph × display manipulation × flower`. The ideal contrasts involve several genuinely independent gynodioecious/dioecious/hermaphrodite clades with a shared, harmonized species tree; sympatric paired populations and repeated seasons are preferred. *Daphne* examples cannot substitute for independent evolutionary replicates.

**Intervention**: Randomly choose modest floral-display reduction or augmented standardized artificial floral visual displays / supplemented intact flowers within individuals, with appropriate sham controls. Do not confound flower removal with plant energy allocation: the manipulation must separate **visual display**, **nectar/odor reward**, and **resource cost**. Record disturbance and damaged flower physiology. Compare paired matched plants or repeated cycles when carryover effects are modeled. Genet identity is essential for clonal plants.

**True floral investment**: Dry mass/flower (non-destructive subsample), tissue nitrogen/sugar/nectar, number of open flowers, days of flower exposure, whole-plant size (stem/allometry or calibrated destructive subset), leaf/branch/fruit mass allocation by sex. Woody canopy AGB from GEDI/ESA and NDVI are **not** focal *Daphne* individual biomass.

**Exposure and pollinator process**: Continuous instrumented time windows + human observations by day and night, visitor guild, weather, 0-visitor intervals, exact flower-hours, first-visitor waiting time, shared floral patch context and conspecific alternative flowers. Separate a two-part observation/Poisson/hazard model:
  - `P(any effective visit | open flower, effort, compatable partner)`
  - `E[n effective visits | >=1, flower exposure, guild]`
  - a zero-inflated or hurdle model only when justified against a properly modeled negative-binomial count baseline with detection.
Use per-window weather and effort offsets and random effects for individual, population and calendar date.

**Functional pollen**: Count stigma conspecific pollen grains / pollen tubes, distinguish incompatible vs compatible pollen by controlled crosses and donor genetics, identify male pollen export and sired seed paternity where feasible. A bee/hoverfly landing is not a fertilization.

**Fitness**: female seed set, 4-ovule per-flower completion if applicable, viable seedlings, individual survival, male siring, repeated reproductive seasons. Evaluate lifetime reproductive success or integrate demographic projection with measured uncertainty; never assume that woody perennial longevity fully compensates a failed year.

**Landscape resources**: remote MODIS EVI/NDVI phenology (year/date matched) plus MOD17 GPP and measured understory photosynthetic irradiance, canopy closure and local floral census; environmental quantity × reliability × accessibility. Use 250 m–5 km buffers as foraging-scale sensitivity, not universal ecological scale. Use real biomass only when instrument footprint/date/forest fraction is compatible.

### Distinct testable outcomes
- **H0 encounter assurance**: additional flowering reduces the chance of a zero-visitor interval at a low-F threshold without changing the conditional positive-visit rate much (as suggested by sainfoin observational pattern). A randomized display modification that changes first-arrival hazard at fixed timing/environment would support a causal encounter effect.
- **H1 display saturation**: at larger F, total effective pollen delivery grows subproportionally; estimate the *fitness relevant* elasticity of compatible pollen/flower, not just visitor counts/structural floral unit. Distinguish pooled landscape differences from within-genet changes.
- **H2 sex-specific returns**: male siring and female seed viability have different relationships with display, resource allocation and compatible mate abundance; morphological male/female differences need not imply opposite display optima.
- **H3 environmental reliability**: the display-to-pollen and pollen-to-seed conversion changes with bloom phenology and temporal rainfall/visitor variability; landscapes with the same annual NDVI can differ in short-term effective pollen transfer.
- **H4 reproductive system and speciation**: quantify genetic mating system, population pollen connectivity and verified reproductive barriers using independent low-copy nuclear loci rather than inferring speciation from current flower size or cp/ITS discordance.

### Preregistered inference boundaries
Do **not** treat 285 site-day observations as 285 independent populations, or 3,048 PoMS fixed quadrat observations as 3,048 plants. Keep the study/site/genet/flower hierarchy; within- vs between-individual decomposition and study-year split are explicit. Fitting multiple correlated outcomes requires error/multiplicity control and negative controls. Reproductive strategy may reflect historical constraints, not current phenotypic selection. Set formal sample numbers by simulation based on a pilot's variance and transition probabilities, not arbitrary fixed N.
