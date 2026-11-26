# Floral investment conversion gates: pollen transfer versus genetic parentage (v0.9)

**Date:** 2026-11-16 JST. **Status:** completed original-public-data reanalysis and validation; not new field observations.

## Publication and source
Joffard et al. 2026, *American Journal of Botany*, DOI https://doi.org/10.1002/ajb2.70276. Public data DOI https://doi.org/10.57745/CF1VNW, v1.0. Official, anonymous Dataverse access with original file member SHA256 pinned:
- pollen_counts_on_stigmas.tab, file ID 765423, 186 rows, SHA256 8291aea33a6c2894edccc55dd911f52a5c9c014049c20fb04ce7397f7e6dee68
- paternity_assignment_results.tab, file ID 765422, 1,150 rows, SHA256 0ae809e19d4b19f60e80ed6e81c2985ecdeb0882b7b857c449b567b325462a50

The two source tables are **different experiments**. No original bee, plant, offspring, father, flower or geographic identifier was published into our project. Inspectable code:
- scripts/joffard2026_dataverse_source_inventory_v09.py
- scripts/joffard2026_pollen_paternity_primary_refit_v09.py
- CI angiosperm-joffard2026-original-v09.yml and angiosperm-joffard2026-original-stage-v09.yml.
- Source access audit https://github.com/zuizui0223/onishibari/actions/runs/43791048111 SUCCESS.
- Completed model plus sensitivity https://github.com/zuizui0223/onishibari/actions/runs/43791675028 SUCCESS.

## A. Experimental sequenced pollen transport: 186 stigmas, 25 bumblebees

Controlled pollen-laden bumblebees sequentially visited flowers of female Silene dioica; pollen grain counts were measured microscopically. Source original: **186 stigma samples from 25 distinct pollinator IDs and 51 female recipient IDs**, visits 1 through 10, 35 zero-grain stigmas and 23,800 counted grains. A pollinator/sequence position repeats once, and individual bumblebee trajectories have only 2–10 observations, so adjacent-position data may not follow the same bees.

| Sequence position | Recorded stigmas | Mean grains | Median grains | Stigmas with zero grains |
|---:|---:|---:|---:|---:|
| 1 | 19 | 529.21 | 435.0 | 1 |
| 2 | 24 | 168.96 | 72.5 | 3 |
| 3 | 22 | 103.64 | 32.5 | 2 |
| 4 | 17 | 98.24 | 25.0 | 1 |
| 5 | 19 | 83.68 | 20.0 | 6 |
| 6 | 17 | 76.47 | 30.0 | 4 |
| 7 | 18 | 67.50 | 42.5 | 1 |
| 8 | 16 | 48.44 | 22.5 | 5 |
| 9 | 16 | 29.06 | 20.0 | 5 |
| 10 | 18 | 21.94 | 5.0 | 7 |

**Bee-fixed-effect, bee-cluster-robust OLS on log(1+grains)** per additional visit sequence position: **beta −0.300076, 95%CI [−0.401101,−0.199051]**, multiplicative *1+grains* factor **0.740762** per position. Bumblebee cluster bootstrap 95% [−0.395248,−0.205108], 6,000 within-bee visitation-position permutations p **0.000167**. Excluding top ~1% high-pollen counts leaves beta **−0.283725 [−0.377423,−0.190026]**. Repeating after removing first-position data (n167) leaves **−0.236983 [−0.356014,−0.117952]**. Leave-one-bee-out beta range **−0.321125 to −0.280016**. Collapsing the duplicate bee-position cell gives **−0.302692 [−0.401757,−0.203626]** across 185 unique cells. Original count-scale Poisson bee-FE IRR 0.706968 [0.622768,0.802552] but **Pearson overdispersion 153.4**, so Poisson serves only as a secondary sensitivity check, not the primary calibrated error model.

**Key limit:** Measured grains on stigma are NOT measured compatible grains, pollen tubes, fertilization, genetically sired offspring or whole-plant lifetime fitness. Visit order within a bee's sequential pollen depletion is NOT an independent manipulation of real flower number or floral biomass. First-vs-tenth pooled means differ with potentially different pollinator cohorts; don't convert ratio of raw means into an individual pollen-loss coefficient. Some flower and bee trajectories are incomplete.

## B. Separate paternal diversity experiment: 1,150 assigned offspring and 21 mothers

The original 2026 journal paper genotyped **1,320 seeds**, but the accessible file has **1,150 confidently assigned** seeds (all distinct IDs). Those records belong to **55 fruits from 21 individual mothers**, of which 10 had female floral display experimentally enhanced and 11 were controls; 77 candidate fathers appear in assigned offspring. There are 575 *assigned records* in each treatment, which is a source sampling property, NOT an estimate of total number of produced viable seeds.

Methods: each fruit's distinct assigned sire count, unbiased Simpson genetic father diversity, effective father number and expected 12-assigned-seed rarefied father richness are computed. Among-fruit father profiles in each mother are compared using Bray–Curtis dissimilarity. Aggregate all fruit outcomes to ONE value per independent mother. Two-sided mother-label permutation n=6,000 and independent mother bootstrap n=3,500; no treating seeds or fruit as independent treatment plants.

| Mother-aggregated paternal metric | Control mean n11 | Manipulated mean n10 | Manipulated–control effect | Mother bootstrap 95% |
|---|---:|---:|---:|---:|
| Distinct fathers per fruit | 6.727 | 5.950 | −0.777 | [−2.043,+0.636] |
| Unbiased Simpson father diversity | .7093 | .6515 | −.0578 | [−.1901,+.0664] |
| Effective father number | 3.7636 | 3.4994 | −.2641 | [−1.3993,+1.0562] |
| Exact 12-seed-rarefied expected fathers | 5.1179 | 4.6908 | −.4272 | [−1.4060,+.5787] |
| Between-fruit sire-profile dissimilarity within mother | .7416 | .7229 | −.0188 | [−.1509,+.1138] |

Mother-level permutation p values, in row order: .311, .400, .700, .439, .803. These five comparisons show **no detectable treatment difference**, consistent with the original publication. But 10 vs 11 mothers, selection/assignment uncertainty and multiple metrics limit power. Genetic father **diversity and composition** are NOT fathered offspring total, plant fertility or actual new-genotype count in the population. No assertion that the effect must be exactly zero.

## C. Generalization across previously audited source experiments

- 2022 Silene dioica randomized **dummy flower visual signal**: observed insects increased (IRR 1.272, 95%CI 1.004–1.610), but maternal viable seed and male genetic siring outcomes had no detected improvement; source hand-pollen control showed no detected female pollen limitation. Thus lack of seed response can be pollen sufficiency, not necessarily ineffective pollinators.
- 2026 independent *S. dioica* pollen carryover: directly counted grains per sequential pollinator visit fall steeply, despite the visual-display manipulation showing no detectable increase in within-fruit genetic father diversity. These are two independent experiments inside the 2026 paper, not a joint visit-to-parentage fate linkage.
- 2025 Clarkia matched fruits: extra hand pollen increased per-fruit seed fill in three of four congeners, but the observational real flower count was not an identified moderator in any of four. One supplementary flower per mother does not identify total plant seed gains.
- 2026 Brassica: many late manually cross-pollinated flowers still fail fruit maturation; within sampled developing fruits ovules often stay unfertilized but most fertilized mature.
- 2025 Veratrum: local flower-level pollen limitation log ratio ~.787 versus whole-genet level ~.129 across different treatment/plant sets.
- 2011–2012 Silene latifolia: original 208 genets link floral census and 2007 survival, but final negative survival association is outcome-time, missingness and leaf physiology sensitive.

**Important novelty limit:** The scientific relationships were already hypothesized or discussed by original source authors; this repository now has raw verification, cluster-appropriate analyses and a common *stage evidence map*, **not an estimate of a single causal numerical flower-biomass-to-fitness pathway**. Source organisms, manipulations and units cannot be mechanically multiplied into one proportion or treated as independent clade replicates.

## D. Distinguishing a real bottleneck in onishibari

The first mechanistically necessary measurements are actual number and dry mass of unique flowers per independently verified genet; real flower-hours; donor/recipient sex morph; within-plant floral advertisement; visited flower ID; counted stigma *conspecific compatible* grains, pollen tube and embryo fate; genetically assigned paternal outcross success and viable maternal seed; future recruitment and adult survival. A pollen-laden bumblebee's sequential deposit decay is not necessarily the decay in *genet's flower-number return* in Daphne. For rare wild *Daphne pseudomezereum*, first verify male-type retained pistil female fertility/self-compatibility with positive cross-pollen controls; retain existing 2027 field protection limit of at most 50% of verified buds and 10 buds per genet across all assays. Independent clade, site, season and within-genet designs are needed to identify the threshold where floral biomass ceases to increase genetic fitness.

**Falsifiable prediction (not observed yet):** At fixed floral construction cost and mean visitor numbers, a plant whose frequent visits mainly redistribute a quickly depleted pollen load among genetically identical flowers may gain fewer unrelated outcross fathers or paternal offspring from an expanded display than a plant with fresh pollen delivery across compatible genets. Test with measured within-genet vs between-genet visits, genetic pollen donors and downstream seed viability; don't call 2026 sequential stigma grains a test of this causal effect.
