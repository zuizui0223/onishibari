# New analysis: floral resource scaling and ecological dilution, PoMS 2017–2023

**Date:** 2026-11-15 JST. **Status:** observed source analysis complete, with named-source quality controls, two PoMS sampling streams and four linked outcomes. Ecological correlational model **only**. Not an experiment and not reproductive success.

## General question
Does greater local flower investment provide proportionally greater pollinator service, or are insect landings a saturating resource? If the latter, a plant may retain large floral display despite highly diluted per-floral-unit visits, while an individual's lifetime reproduction is a separate outcome. This distinguishes floral biomass invested, pollen visitor encounter opportunity, effective pollen transfer, and lifetime seed success.

An existing predecessor must be acknowledged: Grindeland, Sletvold & Ims (2005), *Functional Ecology*, DOI https://doi.org/10.1111/j.1365-2435.2005.00988.x, found declining visitation rate per flower at large floral displays in *Digitalis purpurea*. **Do not claim flower dilution as a newly discovered principle**. The new analysis extends the *pattern check* to community-level standardized FIT Counts across many target floral forms, **without claiming a mechanism is proven**.

## Data and grain
UK Pollinator Monitoring Scheme (2026), EIDC DOI: https://doi.org/10.5285/d706088c-9d48-4be5-b156-6ed6ad5e0618, Open Government Licence with required partnership acknowledgement.
- Allocated-square FITs: 3,053 raw surveys, 3,048 with numeric positive `floral_unit_count`, from 101 1-km squares; one count is 10 minutes observing landings at a selected target flower in 50cm × 50cm.
- Voluntary public FITs: 13,457 raw surveys, 13,372 usable (85 lacking positive floral-unit data), from 2,978 1-km references. **Not a probability sample**; these are two *streams from the same program*, not independent research studies.
- Floral-unit counts are *per target floral unit*: individual flowers, flower heads, umbels or spikes. Median standardized survey F=18; 95th percentile=200, 99th=563, largest recorded count=3000.
- Visitor categories used here: all insects, all bees (bumblebees + honeybees + solitary bees), hoverflies, bumblebees. Insect landings are not conspecific effective pollen deposits, seed set or true per-flower visitation probabilities.
- Fixed selection of positive-flower survey quadrats itself biases away from nonflowering sites: do not infer population-wide absence of pollinators from the fraction of observed 10min zero counts.
- None of these fields measures plant or floral dry biomass, NDVI, reproductive fitness or species-level mating system.

## Statistical parameter
Poisson pseudo-ML for mean visitor counts:
`log E(10min insect counts) = intercept + beta * ln(floral_unit_count) + target flower identity/unit, floral shape, patch cover/context, weather, season + habitat OR square site fixed effects`.
Cluster robust standard errors **by 1-km square** (101 clusters). Check fixed-alpha negative binomial, exclude ambiguous 'Other' target flower (980/3048 records), and drop top 1% of F. beta=0 tests any increase in total landings, beta=1 tests proportionality.
Per-floral-unit *descriptive landing-count ratio* changes by `2^(beta-1)` if F doubles; this is not actual pollen transfer/flower.

### Main structured-square estimates
| Response | Within-square beta [cluster 95% CI] | Total insect landing multiplier for doubling F | Landing counts per floral unit multiplier for doubling F |
| --- | --- | --- | --- |
| All insects | **0.1693 [0.0945, 0.2442]** | 1.124 | **0.562** |
| All bees | 0.2173 [0.1339, 0.3007] | 1.163 | 0.581 |
| Hoverflies | 0.0897 [0.0029, 0.1766] | 1.064 | 0.532 |
| Bumblebees | 0.1443 [0.0239, 0.2647] | 1.105 | 0.553 |

P(beta=1) extremely small for each model, though Poisson mean-variance assumptions are violated (Pearson overdispersion main ~8.7 for total, 2.9 for all bees, 3.6 hoverflies). **Cluster-robust** SE and negative binomial sensitivity matter. P(beta=0) for within-square hoverflies=.0429, bumblebees=.0188, all insects=9.17e-6, all bees=3.3e-7. Interpret borderline small-group p values cautiously.

Other specification checks:
- **All insects** with source habitat rather than square fixed effects: beta=.2081 [0.1537,0.2626].
- Excluding unknown `Other` flower target (no square FE due separation): beta=.2481 [0.1835,0.3127].
- Excluding highest 1% flower-unit cases with square FE: beta=.1595 [0.0834,0.2357].
- Negative-binomial GLM alpha=1 (habitat model): beta=.2093 [0.1579,0.2608].
Sublinear association persists, although standard errors/model and group selection depend on model assumptions.

Source CI and logs: successful CI [run 43705791623](https://github.com/zuizui0223/onishibari/actions/runs/43705791623), artifact ID **14158084034**, ZIP SHA256 `a75a976ed1b199dee627edbb43d8192788d949c6611111636b5a52e623a9bc3b`. Primary script `scripts/ukpoms_flower_visit_elasticity.py`.

### Two PoMS survey stream sensitivity
For a harmonized habitat-adjusted comparison, **same common covariates** on each stream: target flower, flower structure, unit, floral cover/context, weather/season, habitat type, with clustered SE per national 1km grid:

| Response | Allocated-square beta [CI] | Public volunteers beta [CI] |
| --- | --- | --- |
| All insect landings | **.2202 [.1590,.2813]** | **.0916 [.0491,.1340]** |
| All bees | .3063 [.2154,.3972] | .0993 [.0411,.1574] |
| Hoverflies | .1251 [.0381,.2121] | .0625 [.0140,.1109] |
| Bumblebees | .1924 [.0728,.3119] | .0825 [.0295,.1355] |

All estimates suggest *total landings rise sublinearly* with F, but volunteer choice/square/habitat composition differs (public 38.2% `Other` target vs allocated 32.2%). Both survey streams remain correlated by common methodology and country; no claim of independent cross-country validation.

Fixed-alpha negative binomial yields public/allocated all-insect beta .0831/.2101, similar direction. Excluding unknown-target observations also preserves positivity and strong sublinearity (public beta .0997, allocated .2629).

Source [run 43706325276](https://github.com/zuizui0223/onishibari/actions/runs/43706325276), artifact **14158168265**, SHA256 `49b08d73290fda743a103dd9022c08ad3f4ba423e34b4a380c28026ca2110498`.
Script `scripts/ukpoms_public_vs_allocated_flower_scaling.py`.

### At least one visitor versus per-structural-unit service
Additional logistic GLM: P(any visitor during a fixed 10min quadrat) with target+habitat+weather/season and square-clustered robust standard errors.
- All insects: zero-window share 8.86%, log F beta .2899 [.1203,.4596], **odds ratio when F doubles = 1.223**, p=.00081.
- All bees: zero share 58.24%, beta .3113 [.1926,.4301], OR=1.241, p=2.7e-7.
- Bumblebees: zero share 71.95%, beta .1651 [.0463,.2839], OR=1.121, p=.00645.
- Hoverflies: zero share 49.21%, beta .0748 [-.0298,.1793], OR=1.053, p=.161 (inconclusive).
Thus **floral density can make a 10min observation less likely to be entirely zero** while at the same time showing strongly **lower landing counts per structural floral unit**. A probability of at least one observer-recorded landing is not a pollination assurance probability.

Source [run 43706071915](https://github.com/zuizui0223/onishibari/actions/runs/43706071915), artifact **14157788359**, SHA256 `b467567a0745aa61806f455d0d9c6fa66fd5cb30ce1554e811fdcaa5f295a37e`. `scripts/ukpoms_zero_visits_logistic.py`.

### Guild composition difference remains uncertain
Conditional paired binomial counts (hoverfly vs all bees, same FIT windows, only 2056 with either observed):
- Habitat-adjusted `log F` beta −.1577 [-.2831,−.0323], p=.0137, hoverfly-to-bee odds factor .896 per doubling F.
- Within-square FE: −.1329 [-.2818,+.0161], p=.0805.
- Top1% F trim: −.0972 [-.2447,+.0502], p=.196.
Do not claim that bees strongly outperform hoverflies at greater flower density after site/outlier controls. Site and taxon/flower-patch covariates are not experimental manipulations.

Source [run 43705978912](https://github.com/zuizui0223/onishibari/actions/runs/43705978912), artifact **14158686123**, SHA256 `e7748d169832ed4f8010c0de52424e299be831d953f65ff8f43b80c67567e42f`.

## Floral structure & patch-context moderator audit: important identification failure
An exploratory six-test family (Poisson square FE, square cluster robust):
- Open structural flowers: all-insect beta .2032 [.1244,.2819]; closed structural flowers beta −.0242 [-.1249,.0765]; interaction BH q=6.44e-5. Similar point contrast for bees and hoverflies.
- Source floral patch context (isolated vs larger same-flower patch vs larger mixed-flower patch): **no statistically clear effect-slope differences**, joint p=.849 for total insects, p=.438 all bees, p=.389 hoverflies.
**BUT:** The target flower × open/closed cross-tab reveals that 14 of 15 coded target taxon categories are *entirely one flower structure*. The only category containing both is **"Other - please describe below"**, with 671 open and 309 closed. Therefore the apparent interaction by open/closed is **confounded almost perfectly with target taxon identity and source aggregation**; it is NOT evidence that flower opening *causally modulates* crowding/dilution within a plant species. Named-target-only sensitivity fits for all insects and all bees became unidentified/unstable, with only hoverflies converging (wide closed-structure CI). The shape claim must be downgraded to a **comparative hypothesis requiring within-taxon replication or experimental floral manipulation**.

Source [run 43706798707](https://github.com/zuizui0223/onishibari/actions/runs/43706798707), artifact **14158876159**, SHA256 `8100814366ba3904f6e6e7d8f3bf0466d77dcd175dd97c6888b35e1694552387`. Script `scripts/ukpoms_flower_accessibility_interactions.py`.

## What this does and does not say about *Daphne* and global evolution
The empirical result is **nonproportional insect landings per counted floral unit**, where plant taxa, floral unit type and site habitat differ. It cannot show that a large male flower is favored, that a female flower remains self-incompatible, that long-lived plants buffer individual-year pollen failures, or that historical introgression/speciation arose from resource limits.

A testable *mechanism* connecting them:
1. Individual biomass → flower display and lifetime floral investment (need actual dry-mass and size/flower number data).
2. Flower display → total insect visits and distribution of visits among floral units (PoMS shows this is strongly nonlinear across floral quadrats).
3. Insect landings × effective pollen delivery × compatible-mate availability → pollen receipt/export and viable seed set (**not measured here**).
4. Temporal reliability × adult survival → lifetime breeding and pollen-mediated gene flow (**not measured here**).
5. Different lineages/sexual strategies may survive on different sides of these conversion bottlenecks. Whether that generated selection/sexual dimorphism requires new comparative field/phylogenetic evidence.

Prediction: Male display may increase pollen export in absolute terms while yielding diminishing returns per floral unit. Female fitness can be pollen-limited for different reasons (compatibility, fruit/seed resource constraints), so do not equate observed male/female display contrast with measured pollen export/seed set. Experiments should manipulate flowers within both sex morphs at same population, hold area and time, and measure per-flower stigma loads, pollen export, fruit set, and adult survival over years. Stratify remote-sensing NDVI/GPP versus actual **local flower availability**, canopy/urban accessibility and pollinator guild.

### Immediate next designs
- Within-target flower species (several independent clades), flowering phenology/date and individual-genet-level allocation, multiple years, controlled variation in display (natural/experimental) and flowering neighbor density.
- Estimate PoMS `floral_unit_count` elasticity for visitors per structural unit by target flower, but don't infer biology from heterogeneous `Other`.
- Link modern monitoring to biomass/landscape EVI/NDVI and vegetation seasonality only at actual observation dates; previous 2023 UK NDVI heldout gave almost zero incremental predictive power beyond floral units and landcover (run 43704350736).
- Move from counts to actual deposited conspecific pollen/seed output and cross-population paternity for evolutionary consequences.
