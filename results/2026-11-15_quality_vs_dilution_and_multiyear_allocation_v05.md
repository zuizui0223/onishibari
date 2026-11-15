# When does a larger floral display cease to generate fitness? (v0.5, 2026-11-15)

**Status:** 2026-11-15 JST. Independent primary-paper synthesis + *new explicitly hypothetical mathematical counterexample*, checked by GitHub Actions. A measured 3.5-year **Silene latifolia** study provides real floral-mass, pollen-source, siring and survival information. This **does not** supply a fitted joint cross-clade conversion coefficient or new data from *Daphne pseudomezereum*.

## 1. The missing mathematical distinction: flower-level dilution is not loss of whole-plant pollen service

Suppose a single genet produces F open receptive flowers at a fixed time. Let its total **observed landings** be V(F)=a F^beta. Let the probability that each landing is an *independent effective conspecific, compatible pollen delivery* decay as q(F)=q0 F^-gamma. This q is hypothetical, currently NOT measured by PoMS or Zurich (it bundles compatible pollinator behavior, between- vs within-genet movement, deposition and pollen viability). Effective visits then have expected intensity lambda F^(beta-gamma), where lambda=a*q0.

Assume effective events are independently uniformly allocated across the F flowers. In that special Poisson model:
* Mean effective arrivals *per flower* = lambda F^(beta-gamma-1).
* Chance **any effective arrival to a flower** = 1-exp(-lambda F^(beta-gamma-1)).
* Expected total distinct flowers receiving ≥1 effective arrival: 
  G(F) = F [1-exp(-lambda F^(beta-gamma-1))].

**Derivative:** writing k=beta-gamma and x=lambda F^(k-1),
dG/dF = 1-exp(-x) + (k-1)x exp(-x).

For k>=0, dG/dF >=0: since exp(x)>=1+x, we have exp(x)>=1+(1-k)x. Thus increasing F and diluting each flower's *visit rate* does not cause a reduction in total effectively serviced flowers if the number of **effective** visits is nondecreasing. When k=0, G(F) approaches lambda as F becomes large: total success saturates instead of falling. When k<0, G(F) can peak and decline even if **observed visit count** increases with F, because per-visit *effective quality* declines faster than visitor attraction increases.

This is an identity only under the **uniform independent visit assignment** model. Repeated contacts to the same inflorescence, visit-order effects, flower size, phenology, pollen carryover, pollen competition, clonal geitonogamy and incompatible grain deposits violate its interpretation as actual flowers fertilized. A "one effective arrival" event is not a known fraction of fertilized ovules or seed set.

### Completed sensitivity, none of these numerical settings are measured
Use beta=.17 as a hypothetical value similar to the published prior UK PoMS **between-quadrat** β≈.169, not a transferable within-genet effect; set a=4 initial insect landings, q0=.5 **invented** compatibility, lambda=2 and arbitrary development resources B=120, flower cost .4, fruit cost 2, 10 ovules, conditional z=.6, seed development cost .5. The seed index in the model is zoo × min(G(F), positive available development capacity) and is an abstract index, NOT actual seeds.

| gamma (hypothesized decline of per-landing compatibility) | k=beta-gamma | Peak F for G within F=1..300 | Peak F for synthetic maternal-index including floral expense | Mechanism |
|--:|--:|--:|--:|--|
| 0 | +.17 | 300 (still rising at upper bound) | 237 | G rises even though visits per flower decline; resource budget eventually binds |
| .17 | 0 | 300 (approaches plateau λ=2) | 275 | Effective contacts stay constant; finite budget determines eventual decline |
| .27 | -.10 | **9** | **9** | Effective pollen service itself declines beyond 9 because per-visit compatibility deteriorates faster than attraction |
| .47 | -.30 | **3** | **3** | More dramatic quality bottleneck |

These "optimal F" values depend completely on arbitrary scenario choices and are **NOT predictions for onishibari or Silene**. If only observed visitors are counted, beta can be positive throughout, so one cannot diagnose the quality collapse without **compatible pollen delivery per visitor**.

Reproducible script: \`scripts/effective_visits_flower_coverage_frontier.py\` and \`.github/workflows/angiosperm-effective-visit-quality.yml\`. Analytic derivative checked vs finite differences, nondecreasing G for k≥0 verified numerically F=1..10000, declining examples for k<0, and separate finite reproductive-budget cap validated. **SUCCESS** https://github.com/zuizui0223/onishibari/actions/runs/43723236266 ; current extended workflow with independent Silene source arithmetic also **SUCCESS** https://github.com/zuizui0223/onishibari/actions/runs/43723425050 . This is model validation, not evidence about the mechanism in the field.

## 2. Real multi-year case with actual flower mass and lifetime carryover: Silene latifolia
**Primary publication:** Delph & Herlihy, online 2011 / *Evolution* 66(4) 2012, DOI https://doi.org/10.1111/j.1558-5646.2011.01510.x . Original source data (five files of individual IDs/flower censuses/fruit, seed, siring, viability; metadata inspected but individual files not downloaded): https://doi.org/10.5061/dryad.20k62r83 . The original study planted artificial-selected **small- and large-flower lines** and monitored performance for **3.5 years**, including PGI-based paternity assignment of 3,560 seeds *to lines*, 2004 and 2005 flower counts, physiological traits and survival to summer 2007.

**Source-faithful quantified differences:**
- Original Table 2 male flower mass: large-line **10.81 mg (SE .532)** vs small-line **7.18 mg (SE .297)**, ratio **1.5056**. The original table says "flower mass", but do not automatically label it calibrated **total flower dry-mass investment**.
- Earlier Herlihy & Delph (2009) line-comparable estimated male pollen grains per flower, as quoted in the 2011 paper's fecundity section: large-line **29,750** vs small-line **13,040**; ratio **2.2814**.
- Dividing the **different-study** pollen/flower estimates by line-mean flower masses produces a *constructed cross-publication index*, **2,752.08 vs 1,816.16 pollen grains per mg flower, ratio 1.5153**, in favor of larger-flowered line. THIS IS NOT a confidence interval, same-individual mechanistic efficiency estimate, or proof that an extra milligram causes more pollen.
- 2004/2005 male flower *snapshot sums* across twice-weekly surveys do not statistically differ between lines (p=.936/.592), although **early-versus-late flowering distribution** differs. The count may include repeatedly observed flowers and is **not unique flower production**. In fact authors reported 1.18 fruits / counted flowers in one female line-year, highlighting that a snapshot sum is not a fruit-set denominator/probability.
- Line-level siring share was *not detected to be different from the expectation based on estimated pollen supply* (2004 mean share difference -0.02, p=.65, 2005 +.027, p=.56). Absence of detected disproportion is not proof of exactly equal rates or exact genetic fathers.
- Large-flowered females produced greater **total seasonal seed mass** in both 2004 (t=2.96, p=.004) and 2005 (t=2.51, p=.015), despite having fewer observed flowers. More seeds per fruit was seen in both years. Female seed mass is not established seedling recruitment.
- Four sex×line groups had different survival over 3.5 years, overall log-rank χ²=10.73, df=3, p=.013. The large-flower line tended to survive better; the study implicated **genetically correlated leaf physiology and water use**. It did **not randomly manipulate floral tissue expense**, so does NOT show flower creation directly caused a death.

These are all original authors' findings and published means, **not a new causal reanalysis**. \`scripts/silene_latifolia_tradeoff_source_audit.py\` checks exactly specified source values and arithmetic; current successful CI https://github.com/zuizui0223/onishibari/actions/runs/43723425050 .

### Important identifiability improvement
The *Silene latifolia* experiment has something stronger than a one-year flower-vs-seed cross-section: **pollen supply, floral mass, female seed production and subsequent survival** were collected for the sex-selected line system. However the experimental manipulation is selection-line ancestry and not "add one flower to a fixed genotype while controlling physiology." It is therefore an excellent direct comparative case and a poor **isolated causal cost-of-floral-biomass** estimate. Its published paternal assignments are to **lines**, not individual genets, so don't silently transform line paternity into individuals' realized fitness.

## 3. Complementary intervention evidence, with no fictitious joint raw data
- **Silene dioica**, Barbot et al. 2022, doi 10.1098/rspb.2022.1987: artificial (nonrewarding/non-gamete-bearing) fake flowers increased total insect attraction without detectable increase in maternal seed or paternal siring success. This isolates *visual advertisement* better than adding real flower-producing organs, but null statistical significance is not proof of zero; same genus as S. latifolia, **not a phylogenetically independent clade replicate**.
- **Mimulus ringens**, Karron & Mitchell 2012 doi 10.1093/aob/mcr193: experimental flower display size alters male and female returns differently (published total male siring exponent .46 vs total mothered seeds exponent .98). Real flower number is not dry-mass manipulation.
- **Brassica rapa**, Harder & Hamnett-Mount 2026 doi 10.1002/ajb2.70204: controlled hand-pollen fractions yield post-pollination maturation limits; in sampled developing fruit, only 61% of ovules fertilized and 96.3% of fertilized ovules developed. The 2026 model includes development cost and resulting reproductive limits; floral construction-cost extension in this repository is expressly hypothetical.
- **UK PoMS** survey target-floral-unit elasticity beta=.1693 within controlled squares (not within-genet randomized). **Zurich** within-garden/date NB betas near proportional, showing scale dependence. No effective pollen-per-landing and no floral tissue cost measured.

One must not multiply or meta-analyze the pooled visitation beta, conditional fertilization fraction, and separate male siring elasticity: they have incompatible units, treatments, experimental plants, seasons and sampling grains.

## 4. Pre-registered biological tests that can distinguish three kinds of failure

**A. Delivery vs dilution (parameter gamma).** For each independently identified genet in a population and flowering day, manipulate *visible artificial floral display* at fixed genuine receptive flowers and sham controls, then record equal day/night observation hours, total landings V, distinct individual flower contacts, compatible conspecific pollen per stigma (pollen tubes or donor genotype where possible) and genetically assigned pollen export in the male morph. **Positive β_visits with γ>β_visits** would indicate deteriorating effective pollen delivery per visitor at higher display; positive β alone does NOT locate failure. Distinct "no arrivals at all", "visitors but no compatible grains" and "compatible grains but fertilization failure" must be coded.

**B. Development, chronology, resource vs hormones.** Independently randomize viable compatible pollen addition and carefully controlled early-fruit load/fruit removal with sham where appropriate, follow *unique* tagged flowers by date, ovules, embryo abortion, mature seed count and germination. A pollen-addition×fruit-load interaction can locate a potential maturation constraint, but removing fruits changes both resources and hormonal signals; it cannot isolate hormone mechanism without further treatments. In sensitive wild Daphne, retain the existing 50%/10 bud/genet limit and prioritize observational/pilot measurements; do not remove developing fruit without permission.

**C. Residual reproductive value.** Repeated annual genet tracking measures survival, next-year flowering and viable descendants. Don't infer lifetime buffering merely because individuals are perennial; estimate an integrated multiyear fitness or reproductive value with detection/uncertainty and genetically appropriate maternal/paternal assignment. Avoid converting repeated open-flower snapshots into unique lifetime flowers.

### Minimal field record schema for causal matching
\`population × genet × year × flower_id × date × sex_morph × assigned_treatment × flower_mass_subsample × flower_exposure_hours × observed_visits × compatible_deposition × embryo_status × mature_viable_seed × male_paternity × next_year_survival_and_flowers\`. Preserve zero-sampling effort as **missing**, not zero-visit observation; destructively sampled flower cohorts cannot subsequently generate fruits, so use comparable independent randomized flowers with a joint missingness model. Parentage genotype IDs/locations must stay private.

### Planned confirmatory contrasts
1. Added *dummy* display causes more visiting insects but no increase in compatible stigma deposits: the **advertisement-to-effective-delivery gate** is directly supported (conditional on proper power and positive controls).
2. Added live pollen causes embryo increases, while mature viable seeds do not increase; fruit-load removal rescues late fruit: a **post-fertilization developmental gate** is supported, though hormones vs resources remain unresolved.
3. Added real flowers increases line/genet siring even when dummy signal does not: pollen supply, mate access or display-independent mechanisms beyond image attraction contribute.
4. Current siring rises while adult survival/next-year offspring fall in the same lineage under randomized manipulation: real **interannual reproductive-value trade-off**; otherwise a selected-line comparison remains confounded by physiology.

## 5. Data access and integrity
Dryad metadata currently exposes the 2011 S. latifolia five original CSV filenames (Flower counts 2004, flower counts 2005, Fruit and seed data, Individual data, Siring success and male flowering); *Brassica* original seven CSVs, and 2022 Silene dioica three TXT files. The public file-stream web retrieval attempt returned **403**, and no authenticated download or individual-level re-estimation was performed. This report cites original publication details and checked arithmetic rather than inventing individual-level results. The code generates **only hypothetical numerical model curves** and source-value ratios.

**Net contribution v0.5:** The dividing line is not floral saturation alone. A positive total-visitation elasticity cannot cause declining total **effective** service by dilution alone under this model; declining quality (loss of compatible pollen), reproductive developmental expenditure and survival/phenological trade-offs are separately testable mechanisms. Need longitudinal and direct functional transfer, not just flower counts, photos or satellite biomass.
