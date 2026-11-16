# v0.13 | Floral display redirects pollen before it increases maternal seeds

Date: 2026-11-16 JST. Status: source-original five-experiment Clarkia concinna reanalysis plus a separate exploratory post-hoc interaction robustness check. All biological patterns and intervention designs are from published research, NOT new field experiments or independent species replications.

Primary source: Miller & Kay (2025), Madroño DOI https://doi.org/10.3120/0024-9637-72.1.a1 ; original published Dryad DOI https://doi.org/10.5061/dryad.8pk0p2nz2 , version 330787. Five individual original files retrieved via the official anonymous Dryad published-version ZIP endpoint, with source member SHA256 verified. No original plant, pollinator or location IDs exported.

## Original study grains and separate endpoints

| Experiment | Original records | Actual intervention/exposure | Endpoint |
|---|---:|---|---|
| 1 | 1,626 visitor-plant bouts; 270 nominal genus×insect-tag groups | Open flower number observational | Consecutive within-plant floral probes |
| 2 | 137 individually manipulated focal plants; 136 complete stigma observations (34 × 4 arms) | 1,2,4,6 natural pollen-bearing male-phase flowers | Fluorescent pigment from focal self vs neighbor plants, ordinal rank 0–4 |
| 3 | 130 plants assigned (26 × 5 arms); 106 fruits measured for seeds, 105 with both seeds and unfertilized ovules | 0,1,2,4,6 pollen-bearing male-phase flowers | Fertilized seed number and weight per ONE female-phase fruit |
| 4 | 121 sessions, each observing 20 flowers for 20 min | Patch flower number 20–10,000 (observational) | Insect visit frequency |
| 5 | 204 plants, 34 per 3 display treatments × 2 patch groups | Focal flower only; focal + six emasculated; focal + six intact pollen-bearing flowers, crossed with natural/partly thinned small vs large patches | Pollen grains on focal emasculated stigma, first pollen time, corolla wilt time |

These are DISTINCT original experiments, not 1,626→136→105→204 observations of a single tracked flower. Color pigment is not real genetically compatible pollen, and source grain counts in Experiment 5 lack pollen-donor genotypes.

## Consecutive probing and self versus neighbor pollen movement

Experiment 1: 673 of 1,626 pollinator–plant bouts include >1 consecutive flower, and 204 bout plants had only one open flower (no opportunity for a second probe). Exploratory negative-binomial model of the number of additional within-plant probes beyond the first, adjusting genus and date with cluster-robust uncertainty by nominal genus×insect tag: log open flowers beta +1.1314, 95%CI [+1.0429,+1.2199], predicting about 2.19× additional probes when display doubles. This is partly structural and NOT a causal increase of genetically proven self-pollination or arrival of more unique bees.

Experiment 2: four assigned treatments with 34 complete stigmas each. Own-plant fluorescent dye was detected in 9/34 (26.5%), 13/34 (38.2%), 19/34 (55.9%), 17/34 (50.0%) stigmas as the number of male-phase flowers left was 1,2,4,6. The corresponding neighbor-dye detected fractions were 73.5%, 73.5%, 61.8%, 58.8%. HC3 logistic model of own-dye detection on log male flowers: beta +.6445 [.1260,1.1629], p=.0148; neighbor-dye detection beta −.4221 [−.9674,+.1232], p=.1293. Dye deposition does not prove self-fertilization or successful fathering.

Experiment 3: 26 assigned plants for each treatment 0,1,2,4,6 true male-phase flowers. Among sampled mature focal fruits, seed number averaged 37–43 among the arms; adjusted binomial response of seeds/(seeds+unfertilized) to log(1+male flowers), basal stem diameter, floral position and treatment date beta −.2045 [−.5437,+.1346], p=.2372 (n105 complete). Not an equivalence test; 24/130 uncollected/missing fruit seed values are NOT zero offspring. These data do not measure total plant maternal or male genetic fitness.

## Neighborhood patch: visiting frequency and the remaining fertile opportunity window

Experiment 4: 121 observation sessions/20-min/20-flower, 34 zero-visitation windows. Negative-binomial model with calendar date: beta of log patch flower number +.2662 [.1199,.4124], p=.00036, predicting 1.20× recorded visits for each patch-size doubling. Patch exposure is observational and site dependent, not true plant-level floral mass.

Experiment 5: 204 independent focal emasculated stigmas, exactly 34 per arm (1-only, 6e extra emasculated, 6 extra intact male flowers) × large/small patch. Excluding the ONE stigma without grains as censored, 203 positive-pollen flowers have first-pollen timing. Source date and stem diameter adjusted simple OLS suggests small patch log time to first pollen +.3501 [.2036,.4966], p about 3×10^-6, or about 1.42× longer. Original paper instead fitted a Cox model, reaching a similar delayed-pollen conclusion; our OLS excludes one censor and is NOT a fully validated survival/hazard regression. New descriptive remaining lifespan fraction after first pollen = (wilt time − first pollen time)/wilt time was .7156 in large vs .6418 in small; adjusted small-patch beta −.0706 [−.1142,−.0270], p=.00149. This is fraction of observed corolla lifespan, NOT proven receptive time, compatible pollen, or fertilization.

## Targeted pollen-bearing vs emasculated display effect may reverse by neighborhood

Experiment 5 group mean counts of pollen grains deposited on the one focal stigma:

| Source surrounding patch | 1 only | 6e: six non-pollen-bearing flowers | 6: six pollen-bearing flowers |
|---|---:|---:|---:|
| Large patch | 119.85 | 135.85 | 100.41 |
| Small patch | 96.94 | 91.79 | 108.06 |

With the real display size held at six additional flowers, pollen-bearing 6 minus emasculated 6e on log(1+grains) changes from a negative association in large patch (beta −.5844, CI [−1.0034,−.1654], exploratory p=.0063) to positive but individually uncertain in small patch (beta +.3532, CI [−.1225,+.8289], p=.1456).

The unplanned, targeted ONE-degree-of-freedom interaction (6 vs 6e in small minus the same in large) is +.9376, HC3 95%CI [.2952,1.5800], p=.00423. Original authors' TWO-df omnibus interaction test was p≈.06; the narrower post-hoc contrast must NOT be promoted as pre-registered proof. Treat it as a source-sensitive exploratory hypothesis. Separate visual/reward-only comparison of 6e versus 1 was not clearly associated with greater grains in either patch (large +.251 [−.106,+.607], p=.169; small −.098 [−.641,+.446], p=.724).

**Robustness on exactly these original 204 plants**, all exploratory:
- Additional original location fixed effects: interaction +.9957 [.3544,1.6370], p=.00234;
- Negative-binomial raw pollen count GLM: +.4439 [.0044,.8834], p=.04774;
- Winsorized top 5% source grains: +.9325 [.2984,1.5665], p=.00395;
- Date-block bootstrap on unadjusted interaction: +.9249 [.3957,1.4424] using 13 dates;
- Treatment-label permutation within original date×patch cells (8,000): two-sided p=.00313.

Same data under alternate specifications do NOT represent independent replications. Patch type is partly natural / partly thinned and potential food/odor effects of emasculation are not measured. Most critically, a stigma's ALL pollen grains cannot be partitioned by donor with these study-5 records; the difference might involve visitation, self pollen, outcross pollen, or other altered rewards, not simply donor fertilization efficiency.

## Evolutionary interpretation and onishibari next discriminant

Across one species, source-original trials find: more open flowers correlate with more within-genet probes; increasing male-phase flowers causes more own-origin dye on stigmas; increasing male-phase flower number does not yield a clear benefit for the sampled mature seed-containing flower; surrounding flower density can delay pollen arrival relative to floral lifespan; and pollen-bearing versus emasculated additional flowers show an exploratory patch-dependent effect on ALL grains. No stage-to-stage individual causal conversion is observed, and none is true floral biomass vs lifetime genetic fitness.

**Directly comparable original biomass and male-function counterexample (Podolsky 1993, Ecology, DOI 10.2307/1939578):** The short-styled staminate and long-styled fertile flowers of the hummingbird-pollinated shrub Besleria triflora **did not differ in floral biomass investment accumulated until corolla abscission**, but staminate flowers exported *substantially less pollen* over a floral lifetime. Thus, "equal construction mass implies equal direct male reproductive contribution" is false in this experimental case. In Podolsky's preceding 1992 Science study (DOI 10.1126/science.258.5083.791), adding staminate flowers to display increased outcross pollen deposition and pollen tubes on fertile stigmas, via attraction and reduced interception of pollen passing through the display. Together these older papers explicitly test **indirect female-function benefits of male-like flowers despite poor direct pollen export**. This is a known andromonoecious within-genet mechanism; do not presume *Daphne pseudomezereum* has equivalent morphology, reproductive system or self-compatibility.

**Independent broader pollen-fate prior research (Johnson & Harder 2023, Proc R Soc B, DOI 10.1098/rspb.2023.1148):** Authors synthesized pollen-removal and transfer outcomes across **228 species**. Their abstract reports mean pollen transfer of removed grains to conspecific stigmas **2.4% for separate pollen monads** versus **27.0% for orchids with solid pollinia**. These are stage-specific cross-species pollen *delivery* measures, not compatible fertilization or genetic fathering, and the 228 spp are original authors' analysis, NOT this repository's newly refitted data. Authors have published the underlying data tables as supplementary Figshare collection 10.6084/m9.figshare.c.6837249. This evidence suggests an unusually useful next comparative axis: source pollen supply → anther removal → pollen delivery → verified fathers, rather than flower area alone.

Direct sexual-dimorphism precedent is Podolsky (1992), Science DOI 10.1126/science.258.5083.791, where andromonoecious Besleria triflora staminate flowers changed pollinator attraction AND transfer to fertile flowers; this mechanism was published long ago and is not evidence that Daphne has that same sex system or functional retained pistil. Morgan (1992), Evolution DOI 10.1111/j.1558-5646.1992.tb00628.x, theoretically analyzed sex allocation to attractive vs male-only and female-only floral structures; Cunha & Aizen (2023), AJB DOI 10.1002/ajb2.16180, already reported across 1,241 animal-pollinated hermaphroditic species that pollen per flower increases with display but ovules do not.

Falsifiable next field design: factorial display-only vs pollen-supplying flowers and controlled surrounding flower density, record individually identified fresh vs returning pollinators, genetic compatible stigma pollen, maternal viable seeds, genetically verified sire offspring and next-year survival on the SAME independently identified genets. Real dry floral construction cost must be measured separately. The focal Daphne pseudomezereum has NO own verified floral-function, SI, paternity or floral dry-mass causal evidence. Preserve existing ceiling ≤50% of confirmed buds and ≤10 manipulated/destructively sampled buds per wild genet; an abundant common-garden model may be needed for larger factorial manipulation.

## Reproducibility

Source table/schema verifier: scripts/clarkia_concinna_five_studies_audit_v13.py; GitHub Actions successful run https://github.com/zuizui0223/onishibari/actions/runs/43830848260.

Primary five-stage reanalysis: scripts/clarkia_concinna_five_stage_refit_v13.py; successful run https://github.com/zuizui0223/onishibari/actions/runs/43831634645; original results artifact 14209345458 with ZIP SHA256 ac0e1d33d41564ab3650a2637fc2ef7b05ad9fe176025e6d8e482d2d7bc77dbb.

Specification stress test: scripts/clarkia_concinna_study5_robustness_v13.py; successful run https://github.com/zuizui0223/onishibari/actions/runs/43832085072; artifact 14209156898 with ZIP SHA256 d526a4351513277ea3011f2167ba5712730e9bdfbc1949d6d726062af7c5869a.

Original archived studies and source evidence are published prior work. Novelty of this project would require an independent same-genet or comparative cross-clade causal test, not merely a new secondary contrast.
