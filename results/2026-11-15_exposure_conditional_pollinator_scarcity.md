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
