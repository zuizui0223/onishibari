# v0.21 | Spatial advertisement versus temporal insurance in Daphne flower investment

2026-11-16 JST. This is a **secondary analysis of published Daphne jezoensis source observations**, not a new experiment or new Daphne pseudomezereum field evidence. Original study Shibata & Kudo (2023), *American Journal of Botany* DOI https://doi.org/10.1002/ajb2.16239. Author's publicly deposited source https://zenodo.org/records/8253921.

## Reproducible source integrity

Original Zenodo ZIP SHA256: 97a403f17ac4025e0a5834485cc9ceb820fa760ed5c38452d4a866abc7d8eec8.
Original FlowerLongevity.csv SHA256: 140aeae0064ff9e13ff4199794abe7dc533d5ac4f3fa0e929fc815eb0d4df017.
Exactly 408 dated records (51 inflorescences × 8 census times) on 14 distinct ramets: 5 female F, 9 male-biased hermaphrodite H. Date strings are month/day with **no explicit calendar year** and were analyzed only as offsets 0,2,5,7,12,15,20,24 days from April 16; the older generic date parser's year 0001 was a software artifact, not a documented observation year.

Nine source FloweringNo entries are the literal string 'na': 5 entries within ONE F outcross inflorescence and 4 within ONE H emasculate_outcross inflorescence. They are **missing**, never observed zero. FlowerNo (initial inflorescence total, 3–14 flowers) is constant at each of the eight censuses. No original plant IDs, coordinates or genotype rows are emitted.

## Operational estimand

Normalized exposure = trapezoidal integral from first to last census of number of open flowers / number of flowers in the inflorescence. Unit: **observed-window flower-days per initial inflorescence flower**. This measures realized display time within the common 24-day calendar interval, not individual marked flower lifespan or real maintenance biomass expenditure.

For literal missing data, present a lower bound by substituting zero open flowers and a conservative upper bound by substituting *all* original flowers open. Interpolation is only an explicitly labeled sensitivity; none of these substitutions are reported as measured outcomes.

## Reanalysis

| Sex-morph / source treatment | Independent ramets/inflorescences | Missing counts | Mean flower-days within 24-d window |
|---|---:|---:|---:|
| F bag | 5/5 | 0 | 22.2533 |
| F natural | 5/5 | 0 | 22.7333 |
| F outcross | 5/5 | 5 in 1 inflorescence | 18.5400 to 22.1400 |
| H emasculate_bag | 9/9 | 0 | 15.1884 |
| H natural | 9/9 | 0 | 16.0847 |
| H emasculate_outcross | 9/9 | 4 in 1 inflorescence | 13.7471 to 15.3582 |
| H self | 9/9 | 0 | 13.5703 |

F bag vs H emasculate_bag, original paper's operational potential-longevity contrast:
- Difference +7.0649 days, ratio 1.4652.
- Plant-group bootstrap 95% interval on difference [+5.6247,+8.6644] days (n=5 F and n=9 H; small-sample, nonexperimental sex contrast).
- Original author paper already found approximately eight additional female flowering days, so this **quantitative normalization is a sensitivity check, not a novel independently discovered sex difference**.
- Bagging F versus emasculation-plus-bagging H differs in handling, so intrinsic genetic sex effect is not isolated.

Exploratory within-ramet treatment label contrasts, not preregistered and not source-verified randomized causal comparisons:
- F outcross minus bag, 5 paired ramets: mean bound **[-3.7133,-0.1133] flower-days**, but one F outcross inflorescence has 5 missing censuses; most plausible fill assumptions must be checked before interpreting a flower-closure mechanism.
- H emasculate_outcross minus emasculate_bag, 9 paired ramets: mean bound [-1.4413,+0.1698] flower-days, 4 missing censuses in one H outcross inflorescence.
- H self minus H emasculate_bag, 9 paired ramets: −1.6181 flower-days, 8 of 9 paired ramets lower. **This comparison confounds emasculation and self-pollen treatment**, so it cannot demonstrate self-pollen-induced early closure or self-incompatibility.
- Within F, natural minus bag +0.4800 flower-days; within H, natural minus emasculate_bag +0.8963. Do not infer higher success in nature from longer open flowers: fruit/seed and pollen reception are not followed within these inflorescences.

## Prior literature and non-novelty boundaries

- Alonso (2004), *Annals of Botany* DOI https://doi.org/10.1093/aob/mch004 manipulated whole-season exposure in **Daphne laureola**: ~50% shorter opportunity reduced F fruit set by ~50% and H by ~25%. This is a real temporal-opportunity intervention, but its time scale is flowering season, not individual flower duration.
- Gavini (2023), *AoB PLANTS* DOI https://doi.org/10.1093/aobpla/plad084: in **Dianthus plumarius**, larger H flowers contrast with longer F stigma receptivity and higher female seed fitness. A similar display-versus-duration idea has already been tested in other genera.
- Wang et al. (2021), *Evolution* DOI https://doi.org/10.1111/evo.14113: in **Cyananthus delavayi**, limited pollen removal prolongs male phase and shortens female phase without changing total flower life. This is a mechanistically distinct, experimentally demonstrated sex-function interference.
- Flower longevity and gynodioecy are established research topics. A manuscript cannot claim sex-specific temporal investment was first discovered here.

## Focal onishibari test

Hypothesis: for sex-function morphs in early spring, male fitness may return chiefly from real pollen-loaded new arrivals/paternal siring, whereas female fitness may benefit from persistent receptive stigmas spanning rare successful visits. Which is more valuable depends on compatible mate and pollinator opportunities, not binary day versus night identity or apparent lack of visual signal.

Measure within identified genetic individuals (not just morphological stems):
1. Actual fresh floral dry mass and opening dates/census-specific alive and receptive flowers for both morphs; cost per surviving flower-day requires a separate empirical respiration/nectar/carbon assay.
2. Newly arriving vs re-visiting identified pollinators, pollen removed from anthers and verified insect-carried grains; total stigma grains versus genetically compatible germinated pollen and mature seed.
3. Paternal offspring identity and distance separately from maternal viable offspring and next-year survival.
4. Bag, sham, hand-outcross and hand-self controls to distinguish developmental retention of pistils from functional female fertility. Do not presuppose physiological SI.
5. Conserve wild rare plants: cumulative destructive or manipulated bud ceiling remains <=50% of verified accessible buds and <=10 confirmed buds per genet, with common garden if experimental sample power insufficient.

**Main scientific boundary:** source 2018 flower-days and 2018/2020 published fruit/paternity measurements are from separate cohorts. They cannot be multiplied into one same-flower causal display-duration-to-genetic-fitness chain. The final new claim needs truly joint measurement or a randomized manipulation in focal species/abundant comparator.

## Execution

- Source script: scripts/daphne_jezoensis_2018_flower_days_v21.py
- Source/sex treatment audit: scripts/daphne_jezoensis_2023_longevity_source_v20.py
- Source original result: https://github.com/zuizui0223/onishibari/actions/runs/43877491700 (successful)
- Machine JSON: analysis_results_daphne_flower_days_v21/daphne_2018_flower_opportunity_auc_source_v21.json (GitHub Actions artifact)
- Draft PR only; no change to main. Results are secondary author-data fits and are not additional independent species for the prior cumulative evidence ledger.
