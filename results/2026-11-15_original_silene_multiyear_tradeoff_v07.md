# Current female seed output versus later survival: original-data Silene latifolia multi-year audit (v0.7)

Date: 2026-11-15 JST. Status: **NEW exploratory reanalysis of publicly released individual CSVs, completed and sensitivity tested; no causal flower-expense intervention and no lifetime parental fitness estimate.**

## 1. Original data and verification
- Publication: Delph & Herlihy (online 2011, *Evolution* 2012) doi:10.1111/j.1558-5646.2011.01510.x; Dryad doi:10.5061/dryad.20k62r83, official **published-version public ZIP**: https://datadryad.org/api/v2/versions/11856/download .
- All **five original CSVs** independently validated by fixed uncompressed member SHA256 via \`scripts/silene_latifolia_original_multiyear_refit.py\`, no source IDs or geography exported. Source records: 208 phenotype–physiology–survival individual rows, 208 plant flower-census trajectories (2004), 128 in 2005, 104 fruit+seed rows, 275 male flowering/paternity **line-level** rows. 2004 census \`Total Flowers\` is a **sum of standing open flowers at repeated survey snapshots**, not distinct flower output, and must NOT be multiplied by mean flower tissue mass as a true allocation budget.
- Selection lines are *large vs small flowers* in each sex; other leaf physiological traits covary genetically. The lines are not randomized expense of an additional flower. Genet and selected line confounding remain. Paternity in this study is assigned to selection **lines** rather than individual fathers.
- Two successful CI original-record audits: initial raw join [#43726938227](https://github.com/zuizui0223/onishibari/actions/runs/43726938227) and expanded timepoint/sex/missingness reanalysis [#43727087333](https://github.com/zuizui0223/onishibari/actions/runs/43727087333). Expanded artifact ID 14166706826; ZIP SHA256 f19293c898eb85d26547d96d30c4c9d734bba7d02014dcfc0bb2199522ea33de.

## 2. Individual/year join audit: no pseudoreplication
- Base **208 independent genets**: males 105, females 103; selected line large 104, small 104.
- 2004 flower count table matches all 208 unique base genets; 2005 table matches 128 unique genets.
- Original 2004 mean sum of repeated open-flower census units by sex×line: female-large 29.29 (n52), female-small 43.73 (n51), male-large 330.40 (n52), male-small 324.91 (n53). **The male and female observations have different morphology/census scales**, and this is not unique flower biomass.
- Last \`Alive 7/2007\` outcome: **27 Y, 165 N, 16 missing** of 208. Explicitly observed final outcome n=192, alive rate 27/192=14.06%. Fifteen of the 16 missing at the last date were originally alive in April 2005; do not equate all 16 missing with prior mortality.
- Original fruit/seed table has **104 entries**, 103 matching female and **one male-coded** by source join. The male-coded row was **flagged and excluded**, not silently relabeled as female.

## 3. Female output: 103 female genets, observational seed-mass association
For 103 females with a matching original fruit/seed record and nonnegative seed mass, a robust OLS of \`log(1+seed_mass) ~ log(1+2004 summed standing-open-flower snapshots) + selected flower-size line\` yields:

**β_display = +0.448393**, HC3 95%CI **[+0.321470,+0.575315]**, p < .001. This is **within the female sex** controlling selected line, but not a randomized display change or direct floral dry-mass investment estimate. The seed-mass source year was **not independently verified** in the exported metadata; do not assert strict within-year chronological mediation.

Selected large-flower females (n52) vs small (n51) have original median seed mass **4.1715 mg** vs **1.956 mg**; median fruit counts **21.5** vs **30.0**. This descriptive contrast shows why a higher fruit number does not necessarily imply higher seed output, but selection-line ancestry, flower size, seed size, maternal resource differences and sampling matter.

## 4. Temporal survival: same 2004 display, five differently timed later outcomes

Explore \`logit(alive_at_time) ~ log(1+2004 census sum) + sex*selected line\`, independent genet and HC3 robust SE. Outcome coding respects original Y/N and does not recode missing as dead.

| Follow-up time | Original n with outcome | n alive | β_display | HC3 95%CI | unadjusted p |
|---|---:|---:|---:|---|---:|
| April 2005 | 208 | 192 | +.344094 | [-.138248,+.826437] | .162052 |
| May 2006 | 208 | 159 | -.029928 | [-.359626,+.299770] | .858790 |
| Nov 2006 | 192 | 83 | -.097622 | [-.445122,+.249878] | .581905 |
| May 2007 | 192 | 46 | -.112400 | [-.499198,+.274399] | .568986 |
| July 2007 | 192 | 27 | **-.408044** | **[-.795186,-.020902]** | **.038849** |

This is a *time-point specific exploratory association*, not a time-varying hazard model. The last association becomes an unadjusted p=.039 only when **27 survivor events remain**, and the earlier timepoints are compatible with null. A conservative five-outcome Bonferroni correction on the last p is **.1942** (BH across this simple family of five is also approximately .1942). The collection of significance tests is not evidence for uniquely timed causal survival suppression.

**Missing data bounds (purely hypothetical extremes, NOT imputation):** with all 16 missing final statuses set dead (208 total), β_display=−.376788 [−.701388,−.052188], p=.022901; with all 16 set alive, β=−.223012 [−.591604,+.145581], p=.235682. The latter changes sign of *statistical conclusion* (sign remains negative but uncertainty includes zero). The missing statuses are not ignorable without a tracking-loss model.

**Sex-stratified July 2007 exploratory model** \`alive ~ log1p(display)+selected line\`: females n95, survivors19, β=−.435968 [−.927760,+.055825], p=.0823; males n97 survivors8, β=−.364106 [−.998005,+.269793], p=.2603. The point direction is similar, but effects are individually uncertain, particularly with only eight male survivors. The pooled sex-adjusted β cannot show the strongest survival mechanism or selection for smaller flowers.

## 5. Interpretation: possibility of an investment trade-off, NOT established cost
There is a credible **candidate** pattern within *Silene latifolia*: the present seed-mass association with summed floral display is positive; the last recorded adult-survival association is negative and appears stronger late in follow-up. But maternal reproductive output and the particular 2004 flowering variable have different source year/measurement grain and are genetically/physiologically confounded. It would be **overclaiming** to say one extra flower causes mortality or to infer an evolutionary optimal allocation from a significance contrast.

This added stage also prevents the error "large flowers have more pollen so males gain" from being equated with total genetic fitness. Sex-specific long-term parentage would require **individual-father assignments**, absent here. The old v0.5 floral mass per flower and earlier-study pollen per flower remain cross-publication aggregates, not joint plant-level flower construction cost.

## 6. What will falsify rather than restate the trade-off?
- On the same genetically distinct plant, measure unique *new* flowers (not weekly repeated open-flower snapshots), **dry mass of floral construction/maintenance**, reproductive resource reserve and after-fertilization fruit/seed biomass. Control for original growth/size before new flowering.
- Randomly reduce true floral resource demand while keeping visible advertisement and available gametes unchanged to separate resource expense from display effect; or use common-garden controlled matched lines. Flower removal changes resource and gamete supply simultaneously, so requires sham/alternative interventions.
- Measure viable maternal recruits and genetic fathers plus survival over >1 year and include multi-state survival/detection; explicit stage-specific fit rather than multiple p tests. *Daphne pseudomezereum* must first receive verified functional sex-type and cross-pollen controls; existing conservation limits ≤50% buds/≤10 buds per genet remain binding.
- A true trade-off would require a within-genet experimentally induced rise in floral investment **raising present reproduction but reducing residual reproductive value** under stable source plant physiology; this dataset does **not** make that causal claim.

## 7. What is now established across independent empirical works?
- **Brassica original 86-plant/2,201-flower/9,421-ovule reanalysis**: artificially hand-crossed late flowers are less likely to develop fruit; fertilized ovules within sampled developing fruits overwhelmingly mature, but many are not fertilized. Original treatment effect and anthesis pattern are separate.
- **Veratrum original 97 pairs**: plant-wide seed gain from pollen addition is much smaller than a flower-only effect; comparable units need matched genet design and handling of zero seed ratios.
- **Silene latifolia original 208 genets**: floral count/survival links can be computed directly over 3.5 years, but final negative association depends on missingness and timepoint scanning. Female seed mass association is positive but not experimentally identified as benefit of floral *mass*.
- **Not established:** universal causal point of diminishing return measured in units of **floral dry mass**, or general sex-specific total lifetime fitness returns. Cross-study differences cannot be chained into one multiplicative percentage.

Original data, code and results are in **Draft PR #3**, main remains unchanged.
