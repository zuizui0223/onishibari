# Same-day floral display, visit scaling, and reproductive output — Zurich 2016 (2026-11-15)

## Research question and novelty boundary
**Does a larger floral display reliably translate to more insect landings and then to viable reproduction?**

Prior study **already** answered many local landscape and pollination-success aspects: Reji Chacko et al. (2026) *Journal of Applied Ecology*, DOI https://doi.org/10.1111/1365-2664.70384, *No one-size-fits-all: trait-dependent effects of local plant diversity on pollinators and pollination service in a densifying city*. Same original EnviDat data https://doi.org/10.16904/envidat.676. This analysis **adds explicit within-garden date-matched flower count × observer time, robustness to season / overdispersion, and normalization of site guild visits by flower-unit exposure**; it is a sensitivity reanalysis, not an independent field dataset or proof of adaptation.

## Source integrity and biological grain
- Authors' published source: https://github.com/mrejichacko/BetterBlooms/blob/main/cleaned_data/sampling_effort_gardens_plants_aggregated_per_date.txt (per garden × field day flower units and observer sampling minutes) and original EnviDat source `06_trait_data/species_temporal_flower_visitation_matrix.csv` (garden × actual observation date × phytometer × visitor taxon × window).
- Original design: 24 Zurich gardens, **19 experimental pots of four species per garden**; planted June 9 2016. Flowers or inflorescences were counted in each sampling round; carrot floral units are **umbels**, the other taxa use their own floral units. Visits observed on dates June15–July20 2016.
- Public source has **91 source garden×date observation rounds**, 21,120 source matrix rows for taxon/window sampling grid and **5,794 total insect landings** before source filtering, including some ineffective/nectar robbing visitors. With positive recorded flowers and observing effort, **285 plant×garden×date** cells (24 gardens; four phytometers), totaling 5,788 matched landings: carrot 72 cells/2,569 landings; radish 71/2,595; sainfoin 67/122; comfrey 75/502.
- **Do not interpret a blank/nonsurveyed date×site×species row as an ecological zero**. Matrix contains 2,112 possible site×date×plant combinations but only 285 recorded positive flower×effort exposures with corresponding observation data.
- Source QA successful [run 43707896582](https://github.com/zuizui0223/onishibari/actions/runs/43707896582), artifact **14158738776**, ZIP SHA256 `7dec704727fba542285004e42173322b04665d349a1c5ae9614b2a5efecef004`.

## Temporal *within-garden-plant* visit–flower elasticity
Observational count model:
`log E(landings_{garden,plant,date}) = log(observed_minutes) + effects(site,plant,date/phenology) + beta_within[plant] × {log floralunits - mean log floralunits(garden,plant)} + beta_between[plant] × mean log floralunits(garden,plant)`.

Because flower counts were sampled once each field observation round, using `floral units × recorded observer minutes` is an **approximate floral-unit time exposure**, not a direct individually tracked flower-hour count. 24-garden clustered sandwich SE; offset observer minutes (not fabricated flower-hour exposure) in primary model. Overdispersion large in Poisson models, motivating fixed-alpha NB sensitivity. Plant- and site-mean centering is not a random experiment: flower number and visitor activity both respond to seasonal weather and development.

| Model | Carrot within β [95%CI] | Radish within β [95%CI] | Sainfoin within β [95%CI] | Comfrey within β [95%CI] |
|---|---|---|---|---|
| Poisson, no date | .837 [.668,1.007] | .596 [.434,.758] | .864 [.275,1.453] | .623 [.318,.927] |
| Poisson, **calendar-date FE** | .670 [.449,.891] | .763 [.499,1.027] | .957 [.252,1.662] | .821 [.512,1.130] |
| Poisson, plant-specific seasonal smooths | **−.238 [−.592,.116]** | .477 [.146,.808] | .906 [.188,1.624] | .628 [.247,1.009] |
| Fixed-alpha NB GLM, calendar-date FE | **1.010 [.679,1.341]** | **.936 [.667,1.205]** | **1.134 [.474,1.793]** | **.871 [.469,1.273]** |

NB day-FE slope CIs all contain the **proportional visitation** value β=1, unlike the UK PoMS simple between-quadrat 2017–23 slopes (~.17 for all insects under within-square controls). **Thus the previous universal “strong per-flower dilution” story is NOT robust to sampling unit and model choice**. In the same experimental gardens, proportional landings may be compatible with observed floral variation, particularly once overdispersion is acknowledged. The carrot beta even changes sign after plant-specific flexible seasonal control, highlighting phenology confounding.

Successful source action [run 43708099704](https://github.com/zuizui0223/onishibari/actions/runs/43708099704), artifact **14158804207**, ZIP SHA256 `3144d1462739ae8ad2bfd4a71577f3c4bbb568ca2fabd9f00b28c9990e6f6fea`.
Script `scripts/zurich2016_within_garden_display.py`.

## End-of-season fruit/seed response after visitor normalization
Original authors published garden-specific pollinator group totals and individual/umbel level reproductive data. Chosen source group for each phytometer follows their functional guild analyses: carrot hoverflies; radish social bees; sainfoin all bees; comfrey bumblebees. Garden 39 omitted (original study's insufficient observation period). Sainfoin also omits gardens 19, 28, 52 for insufficient flowering/sampling per original analysis.

For each garden×plant, `hours=sum(day observer minutes)/60` and `approx_flower_unit_hours=sum(day flower units × observer minutes)/60`. Compare site guild counts/hour and counts/approx floral-unit-hour to an **aggregate seasonal reproductive endpoint**, each plant analyzed separately; no pooling incompatible endpoints. All values below Spearman correlation of **23 gardens each** except sainfoin **20**.

| Species/end point | Guild count per observing hour, ρ | Guild count per approximate flower-unit hour, ρ |
|---|---:|---:|
| Carrot mean seeds per umbel | +.302 (p=.161) | +.059 (p=.788) |
| Radish fruit set fraction | +.447 (p=.0326) | +.377 (p=.0758) |
| Sainfoin fruit set fraction | +.348 (p=.132) | +.350 (p=.130) |
| Comfrey fraction of flowers setting ≥1 seed | **+.593 (p=.00287, BH q=.0287)** | **+.538 (p=.00816, BH q=.0408)** |
| Comfrey fertilized ovules / potential ovules (4 per flower) | +.159 (p=.468) | +.288 (p=.183) |

BH correction across the **ten correlated descriptive comparisons**; no taxon-specific genetic phylogeny, garden spatial covariance or source individual/micro-plant mixed effects. On comfrey fruit/flower seed-setting, partial **rank correlation** (conditioning on published Urban_500 and local insect-pollinated flowering species richness) is +.480 for visitor rate/hour and +.395 for visitor rate/floral-unit-hour. These are neither causal mediation estimates nor a replicate of the authors' hierarchical GLMM.

The flower-hour normalization did **not systematically strengthen** the association with reproduction: carrot became weaker; radish somewhat weaker; sainfoin similar; comfrey fruit/seeded-flower remains positive but declines; comfrey ovule-fill increases but is uncertain. **Visitation saturation proxies are not a universal measure of how much offspring is produced per biomass allocated.**

Successful source action [run 43708479678](https://github.com/zuizui0223/onishibari/actions/runs/43708479678), artifact **14159053751**, ZIP SHA256 `dc17cbd3e6757598f4d90c8dba794ea84e3638a832488d6ff1f33867828623e6`. Script `scripts/zurich2016_display_to_seed.py`.

## Interpretation: ecological opportunity versus resource investment
A useful conceptual series is
1. vegetative individual dry biomass / available energy (direct plant mass; **not NDVI or AGB at landscape**);
2. total floral biomass allocation, flower-unit count, sex-morph and phenology (only flower unit count measured here);
3. **visitor landings per observation time** and visitor arrival clustering (date-matched data);
4. actual deposition of *compatible conspecific pollen* per unit flower and floral reward (not measured);
5. fruit-setting probability, ovule fertilization, seed quality, survival and fecundity (only some endpoints measured);
6. male fitness from pollen export/paternity and long-term population gene flow (not measured).

**There is no single proportionality or saturation exponent that connects all six steps across systems.** The British PoMS between-quadrat F–V sublinearity and Zurich within-garden day response can coexist without contradiction because they answer different statistical and biological questions. The Zürich original paper already showed that pollination success varies among visitor groups and floral specialization, so do not label this broad “trait-dependent responses” concept as newly discovered.

### Specific falsifiable next test
Compare the same flowering genet's flower display (number and dry mass) and effective day/night visitor deposition (compatible grains/stigma per visitor), with male pollen export and female viable seeds under resource-matched growth, at 2+ years and multiple sites. Test whether **the slope of effective pollen delivery/flower with display size** is lower in low-reliability sites (interannual failure and longer no-visit runs), and whether life expectancy changes the demographic impact. This is an experimental and longitudinal question beyond both PoMS and Zürich.

Full code, source SHAs, sample counts and uncertainty are saved on Draft PR #3, main unchanged. Exact private garden coordinates/IDs are excluded from re-published results.
