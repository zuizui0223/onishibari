# Two meanings of pollinator scarcity: number of visitors vs interrupted floral exposure

2026-11-15. Source: original 2016 Zurich field experiment, published Reji Chacko et al. (2025 data, 2026 *J. Applied Ecology*, doi:10.1111/1365-2664.70384). Categorical ecological process audit, **not new field evidence** of adaptations, climatic causation or species-level reproduction.

## Observed data and denominator
Exactly **285** flowering garden×actual-date×plant observations, each with original per-round flowers/flowerheads and positive observer minutes. Dates June 15–July 20 2016. Total 5,788 directly matched visitor landings in those active cells. Do **NOT** interpret nonsurveyed site×date matrix zeros as actual biological zeros; 21,120 source event-window rows include all site/date/plant combinations and many structurally unobserved combinations.

The question is whether visit events are more temporally clustered than expected **after conditioning on the site and amount of currently observed floral opportunity**. For each garden×plant with >=2 observed dates and >0 total visitors:
- Keep its total observed visitor count `N_{gp}` fixed.
- For each active sampled date `d`, probability `p_d = F_{gpd} × effort_minutes_{gpd} / sum_d F_{gpd} × effort_minutes_{gpd}`.
- Simulate 10,000 independent `Multinomial(N_{gp}, {p_d})` allocations. Count site×plant×date zero-visit observations and visit concentration HHI across observed dates; aggregate observed vs expected by plant.
- This is a conditional **constant rate per floral-unit-minute** null. It is not a time-series model incorporating site/day weather, insect flight seasonality, rewards, flower attractiveness, pollen efficiency or behavioral synchrony. Rejections are not directly evidence of ecosystem collapse.

## Actual observed vs conditional-on-flowering null

| Plant | Nonzero-total site-plant groups | Active date observations | Total recorded visitors | Observed zero-visit garden-dates | Expected zero garden-dates (F×effort null) | One-sided p for EXCESS zeroes | Observed mean within-site visit HHI | Null mean HHI | p for HHI excess |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Carrot | 24 | 72 | 2,569 | **7** | **0.4963** | **0.0001** | .5261 | .4965 | .0003 |
| Radish | 24 | 71 | 2,595 | 1 | 1.8882 | .9034 | .5745 | .5760 | .5703 |
| Sainfoin | 22 | 65 | 121 | 23 | 19.9052 | .1316 | .6845 | .6551 | .1804 |
| Comfrey | 24 | 75 | 502 | 6 | 5.9938 | .5956 | .4660 | .4919 | .9179 |

**Scientific interpretation:** the species with the fewest visits, sainfoin, did NOT have statistically evident excess day-level interruptions *beyond a simple null that already accounts for its floral exposure and the observed site total*. Carrot, by contrast, has many visits overall but far more complete no-visit garden-date observations than this floral-exposure model predicts.

This is a **study-specific, conditionally defined discrepancy**, not a universal rank of pollinator scarcity. Carrot open umbels and a largely generalist guild differ from sainfoin and comfrey.

## Why measuring floral exposure matters for null models
As a negative sensitivity, condition on observation effort MINUTES only, ignoring flowering display. This version shows much more apparent temporal overclustering:
- Carrot 7 observed zeros vs **0.0318 expected**, p=.0001.
- Radish 1 vs **0.0154 expected**, p=.0154.
- Sainfoin 23 vs **16.0755 expected**, p=.0021.
- Comfrey 6 vs **2.4605 expected**, p=.0206.

Thus **ignoring flower availability can turn ordinary floral phenology-driven absence into an apparent "pollinator crash" or very rare interaction process.** Conditional-on-F and effort removes the surplus-zero signal for three of four plants in this particular simple null. Carrot remains a candidate for additional weather, activity and visitor-guild/target analyses.

Source run: https://github.com/zuizui0223/onishibari/actions/runs/43708697778 — artifact **14159278802**, ZIP sha256 `e15efcbb9fef61d19cc05ea2a4aa8e6249d1f43c556e6bf3040f28be51a549f1`. Script `scripts/zurich2016_temporal_visit_concentration.py`, source checksum and date-join QC `scripts/zurich2016_flower_date_link.py`.

## Importance to NDVI/biomass and the general theory
NDVI, EVI and GPP are **whole-community vegetation greenness and productivity**, not instant number of compatible flowers accessible to a focal visitor guild. Even accurate forest AGB or plant standing mass cannot identify a **zero-flower day** or the conditional encounter opportunity of a focal plant in 2016, nor separate capture/hour from stigma pollen deposition.

Any cross-species/forest/urban/island comparison of 'scarcity' should distinguish at least:
1. **Floral exposure**: active flowers or flower structures × true observation minutes × sex expression/phenophase.
2. **Mean landing rate conditional on exposure**: # insect arrivals over available floral units and effort, not raw observer zeros.
3. **Temporal intermittency beyond expected**: conditional zero-day fraction / concentration or interannual reliability, explicitly controlling flower availability, day weather and species-specific guild calendars.
4. **Pollination effectiveness and compatibility**: stigma conspecific viable pollen, mating partner/sex phase, seed viability or paternity, absent here.
5. **Life-history integration**: number of breeding years and total reproductive success, absent here.

Each stage can vary independently. A mechanically higher `visit/flower/hour` is not necessarily higher male pollen export or female seed set. Forest canopy biomass can affect insect/flower spatial accessibility differently from landscape greenness. Long life expectancy can accumulate repeated chances without guaranteeing reproduction.

## Critical limits and falsification
- 2016 Zurich experimental phytometer garden network **is the same source as Chacko et al. 2026**: there is NO independent validation of functional visitor scarcity/speciation.
- If the carrot over-zero signature disappears after adjusting measured site-day weather, volunteer identity, nectar status, garden management, visitor groups, or focal flower quality, the stronger "behavioral interruption" interpretation would be falsified. Until then it is residual overclustering under a deliberately restrictive null.
- The null conditions on observed total, which removes between-garden variation and does not estimate adult pollinator abundance or whether plants obtain enough compatible pollen.
- Sainfoin contains some gardens with zero visitors for all observed dates; they are excluded from a **conditional allocation** test because it is mathematically impossible to allocate zero events. Excluded sites must be kept in population-level encounter or lifetime failure models.
- Only summer 2016 and 2–5 observed rounds per garden. Cannot estimate **interannual** variance, persistent multiyear interaction losses or lifetime fitness.

## Next required data collection
Mark male/female flower cohorts and measure plant identity and dry mass (and actual flowers), day + night effective pollen transport, seed set and genotypes over several years. Keep the flower-exposure-adjusted no-visit probability and actual pollination deposition as distinct outcomes. Include within-site floral season and habitat/NDVI metrics, but avoid using NDVI as the primary measure of floral availability.


## New diagnostic: non-linear floral opportunity vs simple day-level rarity
A follow-up sensitivity checked the conditional observed-zero test by assigning each within-site day's visit probability according to `p_d ∝ observed_minutes_d × flower_units_d^alpha` rather than only `alpha=1`.
Carrot has seven complete zero-visitor garden×date observations, **all at 1–3 umbels**; four at exactly one umbel, six at 1–2. Every carrot-zero day still had visitors to at least one **other** flowering phytometer in the same garden/day, and radish received insects on all seven. Median recorded visits across other phytometers were **50** on a carrot-zero day, so those sites were not devoid of insect activity. This does **not** prove compatible carrot pollinators were available or effective.
- Carrot zero days expected if alpha=0.5: .0308, p_excess=.0001.
- If alpha=1: .5038, p_excess=.0001.
- If alpha=1.5: 2.5899, p_excess=.0012.
- If alpha=2: **7.542**, p_excess=.7163.
- If alpha=3: 17.557, p_excess=1.
Thus an imposed **superlinear floral-display attraction weight** could reproduce the *number* of carrot zero days, but alpha=2 **overestimates visit concentration** (simulated mean HHI .633 vs observed .526), so it is not a jointly fitting ecological process model. No direct evidence alpha=2 is a biological mechanism.
- The apparent carrot tiny-display/nonzero association is also confounded by common flowering DOY, as the conditional garden permutation check below makes clear.
Sources: nonlinear-exposure null [run 43709205090](https://github.com/zuizui0223/onishibari/actions/runs/43709205090) artifact **14160205286** ZIP SHA256 `537dfcf588a9ae60ec8c26996bdb2e7260b58764decf7a271db100aa2dd9d57d`; zero-site alternative flowers [run 43709226965](https://github.com/zuizui0223/onishibari/actions/runs/43709226965) artifact **14159239656** SHA256 `0f7f1f669871d16ef3cee57a553501e0f08f0e3d5c4515c3d14abb0b971786ce`.

## Two distinct models: an observation day receiving any visitor versus visit intensity after arrival
Exploratory split models at **same flower-bearing site/day**:
- Part 1: `Pr(at least one captured visitor)` binomial logit on ln flower abundance, recorded effort + DOY, 24-garden clustered SE; **NOT** true pollination assurance.
- Part 2: positive-only Poisson QMLE for total captured visitors, with observation-minutes offset; **NOT** a fitted zero-truncated hurdle likelihood. Do not draw formal hurdle mechanism claims.

**Sainfoin** (67 site-days, 24 zero visitor site-days):
- F=2–3 source units: 0/3 dates with recorded visits;
- F=4–9: **10/22 (45.5%)** any visitor;
- F≥10: **33/42 (78.6%)** any visitor.
- DOY-adjusted part-1 logit beta=1.8699 [cluster95 .8675,2.8723], OR for doubling flower units **3.655** (p=.000256, exploratory). Part-2 positive-only count elasticity beta=.2145 [−.1449,.5739]; the latter is **not distinguishable from zero** under this model.
- This fits a plausible **encounter occurrence bottleneck**: having more flowers may improve chances of at least one landing without guaranteeing more landings after a site/day is visited.

**Comfrey** (75 date/site, 6 zeros): adjusted part-1 OR 2.236 per double display; positive-only count elasticity .4419 [.1191,.7647]. But only six zero days across all sites and limited within-site contrasts.
**Carrot** (72 site/day, 7 zeros): unadjusted part-1 OR 3.327 per double F, but after DOY adjustment it changed sign and was uncertain (beta −2.865, CI −6.905 to +1.175), so **do not claim a robust display threshold**.
**Radish** (71 site/days, only one zero): logit coefficient is **UNIDENTIFIABLE due quasi/completely separated outcomes**; the script explicitly invalidates huge estimates instead of showing them as statistical results.
Source [run 43709547786](https://github.com/zuizui0223/onishibari/actions/runs/43709547786), artifact **14159379741** ZIP SHA256 `7bc9e7ca8def241f454d9501aa0465c8444d7af27718f9256be0492007680423`, `scripts/zurich2016_display_two_part.py`.

### Within-garden randomization of observed zero-visit day labels
Condition within each garden/plant on its observed number of zero vs nonzero days; randomly permute the labels 10,000 times. Equal-weight informative gardens. Time-detrended sensitivity removes a **global plant-specific quadratic DOY trend** from log flower counts before permutation. No garden-specific weather or observer correction.
- **Sainfoin**: 17 informative gardens, 53 site/date points. On observed positive vs zero-visitor days, raw log-display difference +.5843 (exponential ratio **1.794**), two-sided p=.0072. After quadratic DOY de-trending, residual log-display difference +.5348, **two-sided p=.0108**. The latter is *residual seasonal-adjusted log F*, not a literal between-days flower-count ratio.
- **Carrot**: 7 informative gardens, 25 site/date points. Raw difference +1.4636 (raw flower ratio **4.32**), two-sided p=.0009. After quadratic DOY detrend difference **−.2994** (two-sided p=.1785): the apparent raw signal is not robust to phenological timing.
- **Comfrey**: 5 informative gardens, 17 points, raw two-sided p=.1066, DOY detrended p=.0637; too few informative site contrasts for a firm interpretation.
- **Radish**: only 1 informative garden (almost every active sampling day had ≥1 recorded landing): cannot run meaningful within-garden label permutation.
Source [run 43709729983](https://github.com/zuizui0223/onishibari/actions/runs/43709729983), artifact **14159443455** ZIP SHA256 `637262f18839754c9325c5290312a9fe294fbb0c50b44dd88f96dae84bdee159`, `scripts/zurich2016_zero_display_permutation.py`.

### Reconciled interpretation
- Under the *constant rate per floral-unit-minute* conditional null, sainfoin's **23 zero days** vs expected **19.9** were not excessive (p=.13). The same species' **relative flower display** still predicts which days receive a visitor in within-garden permutation. These facts **are not contradictory**: flower availability can explain many observed zeroes, while seasonal/study site interactions still determine when visits occur.
- The original 2026 paper already found trait-dependent visitor communities and seed output. Our contribution is **method-sensitive decomposition** into display, occurrence, and conditional intensity; not a novel experiment or strong evolutionary causal claim.
- The appropriate research hypothesis is that **biological partners may be rare because flowering opportunity is absent/low, because interaction arrival is stochastic beyond exposure, or because effective compatible pollination fails even after arrival**. These must be modeled separately across temporal grain and floral morphology. For sex-system evolution, add floral dry biomass, male/female functional fitness, and repeated years in independent lineages.
