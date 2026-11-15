# Vegetation greenness is not independent of urbanization: new 2016 RS × visitor-guild + biomass audit

**2026-11-15 — externally sourced, reproducible exploratory analyses.**
Broad research question: how much of the apparent reproduction/visitor-environment association is (i) landscape structure / urbanization versus (ii) independent remotely sensed vegetation quantity, temporal reliability, (iii) direct local floral resources, and (iv) individual biomass allocation?

## A. Same-year, same-garden MODIS × flower visitor data — actual imagery retrieved
Source: Reji Chacko et al., EnviDat 10.16904/envidat.676 (2016 Zurich four standardized phytometers, 24 gardens, 21,120 visitor windows, 5,794 recorded visitors); NASA MOD13Q1.061 250-m 16-day Terra composites from June 9, June 25, July 11 2016 through Planetary Computer STAC. All 24 gardens occupy distinct 250m pixels and 3 composites passed pixel reliability QA at each garden.

- Initial four plants / 2 indices 24-garden test: positive NDVI association with carrot total visitation (rho=+.4616, unadjusted p=.0232), sainfoin (+.4949, p=.0140), radish (+.3855, p=.0629), comfrey negative (-.1613, p=.4515); **none FDR q<.05 for eight tests**.
- Group-specific 167 taxon mapping from the original author-supplied checklist: **100% of 5,794 recorded visits successfully matched**.
- Carrot guild counts: 67.95% bees, 27.55% hoverflies, 4.51% wasps. Sainfoin 94.31% bees, only 123 recorded visitor events across 5280 windows. Zero-visitor windows by plant: carrot .9201, comfrey .9487, radish .9205, sainfoin .9833.
- **Positive exploratory group-specific signal:** Carrot hoverfly visits per observation window vs NDVI: Spearman rho=+.65913, p=.000460; FDR q=.01242 across 27 nondegenerate plant × index/guild-diversity comparisons. Carrot bees vs NDVI rho=+.2506, p=.238. This does not imply plant-specific effective pollination by hoverflies.
- Geographical robustness of the raw carrot hoverfly association: leave-one-garden-out rho range +.619 to +.754; garden-bootstrap percentile 95% interval +.305 to +.846; rank-based linear latitude/longitude residual correlation +.714. This **does not** address observed garden urban cover or spatial autocorrelation (and is post hoc).
- Source runs: https://github.com/zuizui0223/onishibari/actions/runs/43701737786 and https://github.com/zuizui0223/onishibari/actions/runs/43701970972. Artifacts 14157385078 and 14157300861. No garden coordinates or IDs in uploaded results.

## B. Decisive confounder check: published urban/sealed ground versus NDVI
The 2026 *Journal of Applied Ecology* paper by Chacko et al., DOI **10.1111/1365-2664.70384**, already analyzed **these same 24 gardens**, the 500-m sealed-surface (Urban_500) gradient and local floral richness, finding densification decreased visitation/richness and floral enrichment affected certain guilds, but not hoverflies. **Do not claim our raw NDVI–hoverfly association is a new causal mechanism.** The authors' reproducible sources:
- https://github.com/mrejichacko/BetterBlooms/blob/main/raw_data/explanatory_variables.txt
- https://github.com/mrejichacko/BetterBlooms/blob/main/raw_data/plant_floristic_data.txt
Both provide a one-row-per-garden source covariate table; all 24 IDs matched original EnviDat in memory, private garden coordinates excluded from released outputs.

New joint analysis, identical gardens, original 2016 site-season MODIS, fresh preexisting published covariates:
| Descriptive statistic | Value |
|---|---:|
| NDVI vs sealed-surface proportion (Urban_500), Spearman rho | **-0.92696** |
| NDVI vs local May–August flowering plant species richness, rho | **-0.03533** |
| carrot hoverfly visits vs NDVI, raw rho | **+0.65913** |
| carrot hoverfly visits vs Urban_500, raw rho | **-0.70307** |
| carrot hoverfly visits vs local flowering species richness, raw rho | **-0.07355** |
| NDVI ↔ carrot hoverfly visits, *partial rank rho* controlling Urban_500 and flowering richness | **-0.15451** |
| EVI ↔ carrot hoverfly visits, *partial rank rho* controlling same | **-0.04575** |

Predictive skill (leave-one-garden-out, ridge alpha=5, outcome log1p(carrot hoverfly visits per recorded window)):
| Features | MAE | LOOCV R² |
|---|---:|---:|
| Urban_500 + local flowering richness | **0.047826** | **0.32994** |
| Baseline + NDVI | 0.048875 | 0.29560 |
| Baseline + EVI | 0.049244 | 0.29323 |

Adding NDVI increases LOOCV MAE by **2.19%**; EVI increases by **2.96%**. This is consistent with MODIS greenness being primarily a proxy for less sealed cityscape, **not an independently predictive floral resource effect at these 24 sites**. Near-collinearity rho=-.927 limits any attempt to isolate a unique NDVI causal slope.

For radish and comfrey the added NDVI also did not improve LOOCV MAE; no unique NDVI effect proven in this case study. Densification and flower richness were already explicit and quasi-orthogonal in source design, which explains very weak NDVI–local flower richness correlation.

Run: https://github.com/zuizui0223/onishibari/actions/runs/43702182669 — artifact **14156434500**, ZIP SHA256 `16cfb5476a7292c58f45205b48f1ab5782f2a0865f2b9408ae07addafabefb2f`. Source: `scripts/zurich2016_ndvi_added_value.py`.

## C. Individual biomass is a DIFFERENT ecological outcome than greenness
Direct measured woody mass, BAAD v1.0.0 Falster et al. 2015 DOI 10.1890/14-1889.1, m.to is in **kg** (whole-plant dry mass). Joined species-level Wang 2021 sexual system labels via exact binomial names:
- 21,084 woody individuals original; 8,878 with valid total dry mass; 4,562 mass+sex matched across 115 species; 4,133 mass+height+sex individuals across **91 species** (60 hermaphroditic; 12 dioecious; 19 monoecious).
- Species-median measured mass, medians of these species: hermaphroditism .00541 kg, dioecy .00107 kg, monoecy 2.945 kg. **Do not attribute this to sex systems:** individuals from different ages/sizes and studies; species are not representative of a typical mature plant size; monoecious group strongly enriched in much larger woody individuals.
- Pooled log mass~log height slope=2.0664, but a pooled height-adjusted residual does not control study, family, wood density, habitat or ontogeny.
- Within-family *trait contrasts* are almost absent: only **6 families** have two different sex systems and **zero** have two or more species in each contrasted class. Therefore no robust within-family estimate of sexual-system mass effect can be made.
- Original data do **not measure floral dry mass or male/female reproductive investment**. Need direct flowers/fruits/seed mass, adjusted for plant biomass, sex morph and phenophase before claiming reproductive allocation patterns.
- Run https://github.com/zuizui0223/onishibari/actions/runs/43697541650, artifact **14154593704**, ZIP SHA256 `f861c4172e189367df2f117a1ffc70d2719fa7aeffe17fd10ec39e90689df4fb`.

## D. What this means for the broad ecological question
**Apparent biomass/greenness–reproductive-strategy correlations cannot establish biomass allocation or ecological specialization.**
Defensible decomposition:
1. **Quantity:** individual dry mass, canopy AGB, productivity GPP, greenness NDVI/EVI (not equivalent metrics).
2. **Accessibility:** impermeable cover, understory light/canopy openness, local flowering abundance, bloom–pollinator trait matching and alternate visitor guilds.
3. **Reliability:** within flowering-season and interannual variability, visits per effort, effective pollen per visit, repeated poor seasons, mate compatibility.
4. **Evolution:** phylogenetic origin of sexual systems; measured male/female floral investment vs potential gene flow/reproductive isolation.

Independent conclusions available **now**: (a) 24-garden NDVI and hoverfly visits correlate but collapse after garden urban and flower controls; (b) guild-response heterogeneity; (c) existing direct woody biomass overlap insufficient to attribute sexual system effect. Neither the 24-garden observational study nor BAAD/Wang join establishes selection, evolutionary transitions, speciation or pollination efficiency.

## Next data action
Stop treating 736 historical GloPL coordinates as time-matched to MODIS merely through publication year. Existing weather / rainfall WorldClim sites can provide broad habitat background. To quantify **genuine common environment**, use multi-year NDVI/GPP/RF variance and true flower phenology (with field year), pollinator evidence including pollination efficiency, and flowering male/female biomass across multiple lineages and matched sites. Include canopy AGB × understory light and land-cover × floral-resource interactions. Distinguish independent out-of-region prediction from within-urban garden prediction; a 24-garden NDVI effect should not generalize to whole angiosperms.

## Source licenses/ethics
Garden coordinates are available to researchers under the source archive, but no exact private garden positions or site-level IDs are re-published in this repository's output artifacts. All printed results aggregate across gardens. Endpoints/field years from original published sample files. MODIS via public image services; 2016 scenes and pixel QA recorded in successful workflows. This note is research evidence, not a published study or peer-reviewed causal estimate.
