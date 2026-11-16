# v0.14 | Where does male floral investment disappear? Original 228-species pollen-fate accounting

Date: 2026-11-16 JST. Status: source-original, **secondary analysis** of already published 2023 comparative data; it does not estimate floral dry biomass, genetic siring, an evolutionary optimum, or any experimental reproductive outcome in *Daphne pseudomezereum*.

## Primary-source integrity and valid estimand

Johnson & Harder (2023), *Proceedings of the Royal Society B*, "The economy of pollen dispersal in flowering plants", DOI https://doi.org/10.1098/rspb.2023.1148. Authors' source Figshare collection DOI https://doi.org/10.6084/m9.figshare.c.6837249 ; individual original 228-species Excel DOI https://doi.org/10.6084/m9.figshare.24146916.v1 , file id 42364722, pinned uncompressed raw-file SHA256 e5d981d0a7b3804a7d2309b2df6e1418b86a2d9e52d850175fe69ab5a4144a9c, worksheet "Species medians." This is a species-median synthesis of published studies and authors' own observations. The associated DOCX file is additional methodological supplementary material, not a second species replication.

The source contains **228 species**, of which **202** have measured percent pollen removal and **224** have measured percent pollen transfer efficiency. Exactly **198** have **both** measurements. All calculations below use those SAME 198 species unless explicitly stated otherwise. Data concern grain delivery to a **conspecific stigma**, which can include genetically self-derived, incompatible, nonviable or otherwise infertile pollen. It is NOT confirmed genetically assigned fathers or viable male offspring.

Let RE = (pollen removed from anther)/(pollen produced), PTE = (pollen reaching conspecific stigma)/(pollen removed). Then RPE = RE × PTE = (pollen reaching conspecific stigma)/(pollen produced). For monad pollen, unmeasured loss by wind and rain can contribute to apparent RE, spuriously lowering inferred PTE, as the original authors caution.

## Two non-equivalent definitions of a reproduction-stage bottleneck

**A. Physical produced-grain accounting** (the appropriate measure for "where did the pollen grains go?"):
- Unremoved fraction = 1 − RE;
- Removed but not delivered to conspecific stigma = RE(1 − PTE);
- Delivered fraction = RE × PTE.
They sum exactly to 1 **within each original species**.

On an equal-species (NOT pollen-production-weighted) average over **198 matched species**, per 100 produced pollen units:
- **44.59** remain unremoved;
- **49.11** are removed but never reach conspecific stigmas;
- **6.30** reach a conspecific stigma.

Only **104/198 (52.53%)** species have larger physical post-removal grain loss than pre-removal loss. The physical share is NOT 84.8%.

**B. Logarithmic proportional shortfall accounting** (a different measure of *multiplicative efficiency*):
- Loss intensity from removal = −log(RE);
- Loss intensity from transport = −log(PTE);
- Total loss intensity = −log(RPE) = sum of the two.
The median share of **log loss intensity** attributable to transport is **80.50%** (family bootstrap 95% interval **65.9–97.3%**); **168/198 species (84.85%)** have transport log loss greater than removal log loss. This is a mathematical decomposition, NOT a claim that 84.85% of real pollen grains disappear after removal or a "80.5% selection pressure" estimate. The two bottleneck definitions answer different biological questions and must never be conflated.

Matched 198-species median efficiencies: RE **57.0%**, PTE conditional on removal **9.405%**, RPE **4.227%**. Equal-species mean RPE **6.296%**, geometric mean **2.909%**. Medians of RE and PTE cannot simply be multiplied to obtain a true median RPE, because the quantities co-vary and the median is non-linear.

## Dispersal packaging determines the location of physical losses

These are SECONDARY raw-source group means of produced-grain fates; columns within a row sum to 100%, and N is only the subset with both RE and PTE.

| Original pollen packaging category | Complete species | Families | Unremoved % | Removed, not delivered % | Conspecific stigma received % | Species where post-removal losses exceed pre-removal losses |
|---|---:|---:|---:|---:|---:|---:|
| monads (loose grains) | **34** | **18** | **19.95** | **78.52** | **1.53** | **94.1%** |
| Massulate (orchid segmented pollinia) | 81 | **1** | 31.18 | 62.02 | 6.80 | 75.3% |
| Solid (orchid solid pollinia) | 30 | **1** | 58.46 | 29.24 | 12.31 | 16.7% |
| Asclepiad (milkweed pollinia) | 53 | **1** | 73.05 | 21.77 | 5.19 | 11.3% |

**Critical taxonomy warning:** 111/198 species belong to Orchidaceae and 53/198 to Apocynaceae: **164/198 (82.8%) from TWO families**. The entire Asclepiad column is one family; both Massulate and Solid are from Orchidaceae. There can be no robust phylogenetically independent claim that "pollen packaging causes the pattern" from these three categories alone. 34 monads span 18 families and constitute a more diverse, but modest and incompletely representative, sample.

The original paper's published **2.4% average PTE for monads and 27% for solid orchid pollinia** are estimated marginal means based on different availability/adjustment sets, whereas the **matched-species subgroup medians** in this reanalysis are respectively **0.545%** and **22.69%**. Neither is a mistake; don't compare conditional subgroup medians with the authors' marginal mean estimates as though they had identical estimands.

At the 198-species level, the family-cluster bootstrap 95% interval on physical post-removal *percentage points* ranges about **32.0–80.9%**, reflecting extreme taxonomic imbalance. A family-equal weighting yields mean RPE of **1.67%**, rather than the naive species-equal mean **6.30%**. That is NOT a better universal mean in any taxonomic sense, but a strong sensitivity that "angiosperm overall pollen economy" is not identified by this taxonomically skewed convenience sample.

## Is there a real removal-versus-transfer tradeoff?

Across the 198 source-matched species, **Spearman ρ = −0.37056** between RE and PTE; family-block bootstrap interval **[−0.660,−0.175]**. Crucially this is NOT the algebraic RPE-PTE correlation; it is a genuine association of two measured stages, though a causal tradeoff cannot be inferred.

Within pollen packaging categories, the corresponding descriptive correlations differ:
- Asclepiad n=53: **+0.075**;
- Massulate n=81: **+0.063**;
- Solid n=30: **+0.139**;
- monads n=34: **−0.343** (uncertain in this modest subset).

Further an exploratory family-clustered OLS model of logit(PTE) versus logit(RE), adjusted for original pollen packaging, source grooming indicator, and source geographic region, leaves **n=172** species across **16** families. The logit(RE) coefficient is **β −0.00369**, family-cluster 95% interval **[−0.13065,+0.12328]**, **p=0.955**. Thus the observed broad negative RE-PTE association is compatible with **a strong between-packaging/lineage-composition effect and no within-adjustment linear relation**. It would be incorrect to claim a universal within-species allocation tradeoff or that this analysis disproves the original paper's published marginal treatment results. No full PGLS / experimental manipulation was performed; selection from missing Region/groom is additional uncertainty.

## External independent phenotype bridge and duplication risk

Harder & Johnson (2023), *American Journal of Botany*, DOI https://doi.org/10.1002/ajb2.16177, original Dryad https://doi.org/10.5061/dryad.w6m905qtq published version 250874 includes **51 taxa** with species-level PTE, pollen number per flower and ovule number per flower (a separate article by the **same authors** and potentially overlapping source records). That source was successfully downloaded/verified; its separate exploratory regression and exact GBIF binomial overlap audit are stored in scripts/johnson_harder2023_pollen_ovules_51_overlap_v14.py and an independent workflow. DO NOT treat the 228 and 51 as two independent experiments, nor claim all 228 species have own flower number/flower dry-mass measurements. An exact name match does not establish that PTE came from the same population/year, so inspect any overlap and source disagreement.

**Cross-publication source-overlap result (corrected GBIF taxonomy):** Original accepted binomials use an underscore, e.g. Genus_species, not Genus species. After normalizing, the exact one-to-one accepted binomial join identifies **49 overlapping taxa out of the 51** PTE + pollen/ovule taxa and 228 pollen-fate species. The median absolute log difference of reported PTE percentages on these matched taxa is **0**; just **one** shared species has more than a 10% relative PTE discrepancy. **These ARE NOT two independent datasets containing an additive 228+51 species.** The two unmatched names may be differences in source inclusion or taxonomy; no evidence for two new independent experimental observations.

**Unadjusted 51-species sensitivity (21 named families):** ln(Pollen/flower) vs ln(PTE percentage) family-clustered OLS slope **−0.0965** [−0.6189,+0.4258], p=.717; ln(Ovules/flower) vs ln(PTE) **+0.8831** [+0.3548,+1.4113], p=.00105; ln(Pollen/Ovules) vs ln(PTE) **−0.9796** [−1.1613,−0.7979]. But family bootstrap for the ovule slope spans zero **[−0.0334,+1.3164]**, and the P:O ratio slope is exactly the difference of the two underlying log slopes, not a third independent response. The source's original 2023 AJB authors used phylogenetic comparative methods and other controls, not this secondary raw unadjusted family-cluster fit. The different conditional significance pattern is therefore a **sensitivity/denominator warning**, not a refutation of original published results. No real floral dry mass or actual genetically sired offspring are provided by these data.

The original 2023 AJB article already found a negative evolutionary association between pollen transfer efficiency and Pollen:Ovule ratio with phylogenetic models. Unadjusted raw-species regressions are **not a direct replication of phylogenetically informed estimates**, especially if family composition is skewed.

## Why this advances the onishibari question

Flowers can have four separable direct and indirect male/female benefits:
1. Visual/reward display attracts **independently pollen-loaded visitors** or causes repeated probes;
2. Anthers load pollen onto vectors (**removal**), whose success depends on packaging and reward;
3. Pollen reaches appropriate stigma, may include geitonogamous self pollen, and must be **viable/compatible** to fertilize;
4. Fertilized ovules must become mature viable offspring; male siring and female seeds can saturate or face resource limits independently, with late flowers providing either backup or extra expense.

A simple *flowers × pollen per flower* calculation is not genetically realized male fitness. In the focal *Daphne pseudomezereum*, neither direct paternal successes, floral dry construction mass, compatible pollen tubes nor sexual function of apparently retained female structures has yet been measured in our source records. The mechanism of stomatal interception by male-like flowers established for andromonoecious Besleria triflora is known prior work, not proof of the same mechanism in Daphne. Major empirical next step: **same-genet true vs fake flowers, direct measurement of amount removed vs donor-resolved stigma transfer vs genetically assigned seeds**, with wildplant conservation bounds already defined (≤50% confirmed buds, ≤10 buds/genet total manipulative/destructive).

## Status, reproducibility and novelty

- Manuscript source: Johnson & Harder 2023 doi:10.1098/rspb.2023.1148, https://pmc.ncbi.nlm.nih.gov/articles/PMC10547555/
- Authors' public 228 original dataset doi:10.6084/m9.figshare.24146916.v1
- Public 51 species PTE + pollen/ovule dataset doi:10.5061/dryad.w6m905qtq; original work doi:10.1002/ajb2.16177
- Source audit: scripts/johnson_harder2023_global_pollen_source_audit_v14.py, [Actions #43844539853](https://github.com/zuizui0223/onishibari/actions/runs/43844539853) SUCCESS.
- Source-matched 228 refit: scripts/johnson_harder2023_pollen_fates_v14.py, [Actions #43846329500](https://github.com/zuizui0223/onishibari/actions/runs/43846329500) SUCCESS.
- Separate pollen/ovule 51 species and original taxon-overlap audit: scripts/johnson_harder2023_pollen_ovules_51_overlap_v14.py, through .github/workflows/angiosperm-pollen-ovules-overlap-v14.yml.
- This report is a secondary **data-quality/identifiability** contribution, not a causal claim that pollen packaging was experimentally changed, nor a new species-level 228-sample ecological experiment. No original species rows or localities released in public outputs.
