# Does floral investment fail because pollinators exhaust carried pollen? (v0.11)

Date: 2026-11-16 JST. Status: **new model-shape and withheld-independent-bumblebee analyses completed on original published 2026 experimental data**, not a new independent field study or a fitted universal life-history cost.

## Biological question and causal boundary

Total observed insect contacts are not the same as independent arrivals of pollen-loaded insects or genetically effective pollen transfer. The same bumblebee can visit numerous neighboring flowers as its carried pollen declines. In contrast, flower production across multiple reproductive windows can afford second chances, as shown by the independently published 2026 Sagittaria study. These are distinct experiments, not a single observed treatment chain.

Joffard et al. 2026, *American Journal of Botany*, DOI https://doi.org/10.1002/ajb2.70276; original author-controlled French Dataverse https://doi.org/10.57745/CF1VNW, publicly downloadable tab file ID 765423. Raw tab-delimited pollen counts on sequentially visited female stigmas; original source byte SHA256 **8291aea33a6c2894edccc55dd911f52a5c9c014049c20fb04ce7397f7e6dee68**. The primary paper itself had already established strong pollen carryover / male siring-diversity consequences. Reanalysis here compares alternative *functional curves* out of sample.

## Source quality / independent grain

- **186 original female-stigma observations, 25 independently observed bumblebees, ten possible sequential positions**. A single bumblebee-position combination occurs twice; those two measured grains were averaged, giving **185 unique bee×sequence-position records**. Record missingness was not interpreted as zero pollen.
- 51 female recipient plants occur in the source. The present model uses bee as cluster; does not fully model recipient-genet crossed effects, stigma compatibility, or pollen germination.
- The 2026 separate floral display experiment assigns 1,150 source-offspring fathers among 21 maternal plants; this **does not share the same individual flowers or pollen donors** with the sequential-stigma controlled experiment. Do not fit a cross-experiment causal mediation product.

## What does the carried-pollen decline look like within the same bumblebee?

Outcome **log(1 + counted grains deposited on stigma)**; each predictor is centered within the same bumblebee, thus removing bee-specific baseline pollen load. Models assessed by **25-fold leave-one-bumblebee-out** cross validation: train model on the other 24 bees, then evaluate *within-bee centered deviations* in the excluded bee. This is an external-within-bee *shape* check, **not absolute grain prediction for a completely new bee** (the held-out bee's own response mean is used for centering).

| Predeclared candidate curve | Number of slope terms | Mean held-out-bumblebee MSE | Within-bee R² |
|---|---:|---:|---:|
| Exponential decrease in carried pollen (linear in visit ordinal) | 1 | **2.4414** | .2389 |
| Power-like curve (linear in log visit ordinal) | 1 | **2.3704** | .2579 |
| Initial steep drop at 2nd flower, then separate linear later slope | 2 | **2.3718** | .2690 |
| Continuous piecewise curve with bend at 3rd flower | 2 | **2.4082** | .2602 |

**Exploratory best** is the power-like curve, but the models are not distinguishable given 25 pollinator groups:
- Exponential minus power mean bee-MSE = **+.07092**, paired bootstrap 95% **[−.05692,+.17936]**.
- Power minus second-flower-kink mean bee-MSE = **−.00132**, paired bootstrap 95% **[−.09533,+.10076]**.
- The latter early kink model has first→second change in log(1+grains) **−1.4096 [−2.1898,−.6294]**, and each subsequent step **−.24635 [−.34843,−.14427]**, fitted with bee-cluster robust uncertainty. It is a viable *description*, **not** a uniquely inferred mechanistic threshold.
- Power-like curve has within-bee slope for log(sequence ordinal) **−1.31392 [−1.66879,−.95906]**. This is about the transformed outcome, **not** proof the arithmetic mean pollen count follows a literal power law of exponent 1.314.

All these are analyses of a single source study, **not four phylogenetically independent replications**.

### Zero-grain risk, independent supporting endpoint

Among the 185 unique bee×ordinal observations, **18.38%** had zero grains. A separate exploratory GEE binomial model with exchangeable **within-bee** correlation estimates odds of zero deposition multiplying by **OR=1.2654 per later ordinal position** (95% robust **[1.1265,1.4213]**, p≈7.21×10^−5). The original first-position mean grains (~529 over 19 measured stigmas) falls to position 10 (~22 over 18), but those are **not 19/18 identical bees**. The zero-grain model supports decreasing delivery opportunity, **not** direct loss of self-compatible/compatible fertilizing pollen.

## Interpretation: replenish-versus-revisit prediction

Let B be independent pollen-loaded pollinator bouts, and m the number of successive flower contacts made within each bout, so observed contacts V=B×m. Let g(j) denote deposited grains at the j-th flower within a bout, under the explicitly theoretical assumption that fresh bouts restart the same decreasing g sequence. The delivered grains over the plant's flowers would scale as:

D(B,m) = B × sum_{j=1}^m g(j).

When g(j) is decreasing, **at the same V**, additional independent fresh bouts (higher B, smaller m) yield more pollen than longer within-bout series. This is a **mathematical conditional prediction**, not a causal effect of manipulating flower number or a measurement of g in Daphne. The reproductive consequence depends further on compatible pollen tubes, fruit set, viable offspring and paternal genetic assignment.

Field discriminator for a future onishibari experiment: track (i) new independent pollinator arrivals versus repeated within-genet probes, (ii) actual stigma-compatible pollen deposition per visit/bee position, (iii) total genetic fathering vs seed development. **Total visits and floral construction mass alone do not identify the gate.**

## Contrast with original observed temporal recovery

*Sagittaria trifolia* (Tang & Dai 2026, DOI https://doi.org/10.1002/ece3.73955; repo v0.10) original genotype-matched array showed **later flower/inflorescence fruit-set recovery of ~37.1 percentage points under netting** relative to initial inflorescence within the same plants, while *Brassica rapa* original 2026 hand-crossed late flowers **declined** in fruit success by 24.9–41.4 points across assigned pollination treatments. These contrast distinct *time scales and developmental situations*; the Sagittaria netting intervention is fully confounded with site and microclimate, and no independent bee-bout counts were collected across the earlier/later stages. Thus original evidence supports the existence of both trends but **does not identify the renewal mechanism as their cause**.

## Reproduction / validation

- Source script: scripts/joffard2026_pollen_carryover_shape_cv_v11.py
- CI: .github/workflows/angiosperm-joffard-carryover-shape-v11.yml
- Successful source-refit workflow: https://github.com/zuizui0223/onishibari/actions/runs/43804095718 , artifact 14197329012 SHA256 d8fab3e66bcbe96dc17035343d63687eb3ffa72527674f690d4cd04aa454f0ad.
- Source original 2026 Joffard paternity+carryover previous refits: https://github.com/zuizui0223/onishibari/actions/runs/43792243360 . These show source-labeled male mating-diversity invariance, not actual lifetime siring.
- Five data limitations: bee cluster n25, 51 recipient plants' crossed effects not modeled, position is a controlled sequence not flower-display treatment, deposited grain counts do not indicate compatible pollen, apparent best curve selected among several exploratory candidates. Uncertainty on CV differences includes zero.

## Next decisive original empirical step

Do not infer flower-biomass efficiency from visitation effects or cross-study pollen slopes. In the *Daphne pseudomezereum* focal system, first establish verified male and female **functional** morphs, pollen tube/germination controls and independent genet IDs. Under the preexisting wild-plant intervention caps (at most 50% and at most ten confirmed buds per genet), compare adding a fake advertising signal versus providing real pollen-bearing flowers, while recording *fresh individual pollinator arrivals*, stigma-compatible pollen, genetic fathering, mature viable seeds, and future plant survival. Direct sexual function and pollen compatibility of the apparently retained pistil remain unmeasured.

**Decision:** Stage-gate switching is a viable falsifiable framework, not a demonstrated universal law. The only robust conclusion of this v0.11 reanalysis is clear within-bee pollen delivery loss with repeated contacts and heightened zero-grain risk, *not* the uniquely identified shape of a biological exhaustion threshold.
