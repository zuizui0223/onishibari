# From Onishibari to a cross-angiosperm evolutionary ecology question (v0.1, 2026-11-15)

## Core question (beyond Daphne)
**When can an organism retain phenotypes that appear mismatched to currently observed partners and interactions, and when does intermittent partner availability instead drive reproductive assurance, local extinction, or lineage divergence?**

The focal phenomenon is *not* simply that dioecy and male-biased floral display persist despite rare visits. The deep distinction is between:
1. **encounter rate** (pollinators/other partners observed in short survey windows);
2. **functional success** (effective conspecific pollen transfer per encounter, reproductive compatibility, pollen export, seed production);
3. **temporal integration** (survival and repeated reproduction over many seasons, whether failures covary in space/time);
4. **evolutionary connectivity** (cross-population gene flow and reproductive isolation, which can be nonlinearly related to visitation rate);
5. **evolutionary legacy** (phylogenetic inertia, historical partner community, plastid/nuclear histories, constraints).
A large flower, dioecious morphology, or selfing syndrome is not diagnostic of its functional success, past selection, or speciation mechanism.

## Correct prior simplifications
- *Daphne pseudomezereum* is a **summer-deciduous woody shrub**, not a herbaceous perennial.
- Its apparently low visitor frequency at a field site is not established as low *effective lifetime pollen transfer*; observation effort, daytime/nighttime gaps, paternity and seed viability matter.
- Long-lived plants do **not** necessarily compensate single-year pollen limitation in later years: Tye, Dahlgren & Sletvold (2020) experimentally rejected compensation in two orchids, doi:10.1007/s00442-020-04623-x.
- Xu (2021), doi:10.1111/jeb.13936, explicitly modeled the effect of temporal vs spatial variation on flower longevity, an important precedent. Thus stochastic pollination by itself is **not** a novel hypothesis; the new test must distinguish temporal structure, functional sex allocation and long-term population/phylogenetic outcomes.
- Helmstetter et al. (2026), doi:10.1111/nph.70894, already analyzed 21 reproductive/flower/life-history traits in 361 species and found that dioecy primarily associated with floral investment; simple lifespan correlation may be spurious. Do not claim novelty for re-running sexual system–longevity correlations.
- Wang et al. (2021), doi:10.1111/ele.13815, already investigated global climate and sex-system transitions in 68,453 angiosperms. Avoid repeating basic latitudinal sexual-system regressions.

## Empirical cross-dataset audit (sources fetched, executable scripts)
1. **Wang sexual systems** deposited `SplistWithSexSyst.csv`: 72,459 unique names and labels (55,270 hermaphroditism; 11,546 monoecy; 5,643 dioecy) downloaded from authors' Dryad-replica Zenodo record 4895147, doi:10.5061/dryad.s1rn8pk7n. Original 2021 publication reported 68,453 **analysis species**; DO NOT equate the deposited list with the filtered analyzed set without reconstructing its exclusions and taxonomy. `GloSexSyst.csv` has 464 **spatial grid rows** and is NOT a species table.
2. **PROTEUS/DiveRS** authors' cleaned CSV (`divers-it/rs-traitspace/supplementary_data/cleaned_data.csv`): 364 species, 363 with sex morph classification, 310 with mating category, 330 pollination mode, 354 lifespan, 228 floral reward. Data are 2025 public DiveRS and include recodings/possible missingness, and are NOT equivalent to functional male/female fitness estimates.
3. **GloPL**: 2,969 **experiment cases**, 1,265 unique species, 459 species with multiple cases; authors' source `idiv-biodiversity/pollen-limitation-data-descriptor/Output/GloPL.csv` includes published effect sizes (Dryad doi:10.5061/dryad.dt437). Experiment cases can differ in site, year, morph, treatment or outcome, so 2969 cases are NOT independent species/year trajectories.
4. **Exact species overlap** (casefold, underscore to whitespace, no fuzzy synonyms): Wang∩PROTEUS=251; Wang∩GloPL=625; PROTEUS∩GloPL=45; triple=**28 species**. Therefore full-data-case restriction is strongly biased and insufficient for the intended broad question. Additional taxonomic synonym matching should be separately audited (e.g. POWO) with retained original names.
5. **GloPL × Wang descriptor pilot**: standardized PL_Effect_Size with `Included_Bennett_et_al_meta_analysis=yes` and finite numeric effect sizes gave 611 matched species, 1,583 experiment cases. Species-median log pollen-limitation effects in the overall mixed-protocol pool: hermaphroditism n=526, median 0.24360; dioecy n=40, median 0.11595; monoecy n=45, median 0.15133. These are **unadjusted descriptive values, not evidence of adaptive group differences**. Mix of Sup/Bagout response types, fruit/seed outcomes, incomplete study design, and phylogenetic/geographic sampling bias.
6. **Whole-plant supplementation subset:** among matched species, only **4 dioecious species**, 3 monoecious, 95 hermaphroditic. It is not appropriate to infer how dioecy relates to lifetime pollen limitation from such an unbalanced subset. Whole-plant supplementation during a limited year is also NOT automatically a lifetime experiment. Source paper Rosenheim/2024 doi:10.1002/ece3.10911 warns about whole-lifetime versus flower-level comparisons.

## Central falsifiable model: reliable functional exposure matters, not observed visitation alone
Let `lambda` be effective successful encounter rate per breeding year and `s` the independent post-breeding survival probability (0<=s<1). Under annual Poisson encounters, constant q and geometric survival, the probability of >=1 effective event before death is:
`P_event=(1-exp(-lambda))/(1-s*exp(-lambda))`.
This is an **illustrative probability**, not seed set, fitness or an empirical Daphne estimate.
Compare equal-average-rate regimes: (i) constant lambda for every individual, (ii) 50% of individuals permanently zero, 50% permanently 2*lambda. For lambda=0.1, s=0.95: case (i) P=0.677775, case (ii) P=0.407886. Thus mean encounter rate alone is insufficient; persistent variation and longevity interact. Source code `scripts/lifetime_rare_interaction_model.py`, run 43692850742. This is not proof of demographic compensation: increasing opportunity over lifespan cannot make a lost year of seed production disappear.

## Four new empirical predictions / alternative explanations
**H1, temporal integration:** Differences in annual visits poorly predict differences in *lifetime* pollen delivery and offspring recruitment once adult survival, flower duration, timing, and interannual encounter autocorrelation are measured. Prediction: stronger temporal correlation/long drought of partner events lowers effective lifetime success at fixed mean rate. Alternative: chronic absence + no adult/seed compensation -> reproductive failure and demographic extinction. Needs individual × year × pollen deposition and vital rates, not GloPL pooled cases.

**H2, response non-uniqueness:** Under the same scarcity gradient, alternative solutions (autogamy, day/night partner redundancy, floral investment changes, functional sex bias, long flowering duration, or decline/extinction) emerge depending on ancestral sexual system, resource constraints and encounter covariance. Prediction: phylogenetic effects and prior morphology × reliability interactions explain more than 'scarcity implies selfing'. Distinguish mating self-compatibility, autonomous seed set, and inbreeding depression.

**H3, phenotypic–functional lag:** Clades with similar morphology but different current partners sometimes retain ancestral sexual dimorphism because of historical selection/constraints, *or* because function remains maintained under rare effective encounters. Counterfactual: within a lineage compare common-environment morphs, multiple populations and functional experiments against non-phylogenetic convergence. Avoid calling a morph maladaptive from one year of visitation.

**H4, demographic vs evolutionary connectivity:** Rare successful cross-pollination may have small demographic contribution but disproportionately large genetic consequences (paternity across populations). Prediction: visit count/fruit set patterns and realized cross-population pollen-mediated gene flow need not covary. Test genetically assigned pollen fathers with plastid/seed dispersal separately; no speciation claim from FST alone.

## Multi-scale work package (not a single complete-case data frame)
- **Macroevolutionary background:** Wang 68k/72k species label backbone; audit coverage, taxonomic ambiguity, source, geographic/phylogenetic bias. Infer sexual-state transitions only with an independently justified species tree and observation coding uncertainty; not already published basic climate regression.
- **Experimental ecological bridge:** Wang×GloPL 625 spp for evidence map and hierarchical meta, controlling study, outcome, treatment (Sup vs Bagout), supplementation level, family/phylogeny, species and nonnative status. Whole-plant weak cross-system overlap => *do not fit mechanistic life-history effect* without new data.
- **Multitrait functional comparison:** DiveRS 364-species evidence, missingness-aware, distinguish morphological sexual system from mating system and visitation. Avoid equating 364 species with comprehensive functional sex estimates.
- **New data product:** `species × population × year × sex morph × individual × interaction event` structured tables, target 10–20 clades and multiply sampled populations. Use voucher-confirmed female/male/hermaphrodite floral morphology + conspecific deposition/paternity, seed viability, survival and alternate routes (autogamy/day/night pollination). Record site-level observation effort, weather, partner identities and uncertainty.
- **Species-divergence layer:** sister-lineage contrasts with independent nuclear loci, flowering phenology, overlap, gene flow, reproductive isolation, not inferred from the sexual system map alone.
- **Ecological generalization:** rare seed dispersal, mutualistic animals, host/herbivore and microbial symbiont encounters can be tested under the same interaction-frequency × longevity × temporal-structure logic, but the plant datasets must not be interpreted as universal evidence outside plants.

## Immediately executable next activities
1. Reproduce the **68,453 versus 72,459** Wang release filter; preserve input/analysis-set distinctions and taxon provenance. Dryad metadata notes some variable sexual systems excluded from the published analysis.
2. Preserve **three-level analysis outputs** rather than forced triple-complete cases: `Wang+GloPL`, `PROTEUS`, `new image+field functional sex table`.
3. Add automated QC to protect against denominators 2969 cases vs 1265 species, variable granularity and 4/3/95 whole-plant subset, and against double-counting GloPL site-year cases.
4. Pre-register a sampling design in which pollen/deposition, pollinator observation, fitness and genetic connectivity are measured for the same sex morph and location across seasons.
5. Use new experiments to distinguish long-lived 'buffering' from a false assumption of complete compensation (Tye et al. 2020).

## Source and CI pointers
- Global Wang sexual system source: https://doi.org/10.5061/dryad.s1rn8pk7n
- GloPL descriptor: https://doi.org/10.5061/dryad.dt437
- DiveRS data: https://doi.org/10.48579/PRO/HU3QYM, original analyses: https://github.com/divers-it/rs-traitspace
- Verified source-coverage run: https://github.com/zuizui0223/onishibari/actions/runs/43692209469
- Three-way exact-join run: https://github.com/zuizui0223/onishibari/actions/runs/43692403634
- Pollen limitation exploratory descriptive run: https://github.com/zuizui0223/onishibari/actions/runs/43692755764
- Time-integrated rare interaction model run: https://github.com/zuizui0223/onishibari/actions/runs/43692850742
- This is a standalone **broad** workstream; the original narrow Daphne phylogenomic PR #2 remains separate and unchanged.

## Update: remote-sensed greenness is strongly confounded with landscape imperviousness

A source-reproducible 2016 MODIS × Zürich 24-garden validation now shows **why landscape greenness alone is not a satisfactory pollinator scarcity metric**. NDVI correlates with carrot hoverfly counts (rho+.659, BH q=.0124 across 27 guild tests), but also with the already-published 500m impervious-cover fraction (rho -.927). After adjusting for original urban cover and local floral richness, NDVI partial rank rho is -.155 and adds no leave-one-garden-out prediction (MAE worsens 2.19%). This is a negative result for **independent NDVI explanatory value in this one case study**, not a universal null effect of remotely sensed productivity. The original 2026 paper already tested densification, local flowering resources and pollination success on these exact sites; proper novelty lies in broader cross-landscape quantity × accessibility × reliability comparisons, functional pollen transfer and size-normalized reproductive biomass.

Direct woody biomass test BAAD×Wang: 91 mass/height/sex-system species (12 dioecious, 19 monoecious, 60 hermaphroditic), but only 6 families mix sex-system labels and no family has replicated species in two compared states. Avoid claiming evolutionary size allocation from these descriptive medians.

Verified runs 43701737786 (guild), 43701970972 (sensitivity), 43702182669 (urban flower controls), 43697541650 (dry biomass); integrated evidence `results/2026-11-15_ndvi_guild_urban_biomass_integrated.md`.



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
