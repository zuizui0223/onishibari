# Flower investment is not a single conversion: visual attraction, fertility, maturation, and male vs female gain (v0.4)

Date: 2026-11-15 JST. Status: **published randomized experiments validated; new explicit cost-aware theoretical sensitivity and original-article p-value consistency audit successfully executed; NO new individual-level Daphne, Brassica, or Silene raw fitness observations fitted.**

## 1. Three perturbations identify different links (not one meta-analysis)
| Experimental comparator | Controlled intervention | Empirical conclusion | Limit |
|---|---|---|---|
| Barbot et al. 2022, *Silene dioica*, doi:10.1098/rspb.2022.1987 | Random assignment to visually double floral display with crepe-paper dummy flowers **without new pollen/ovules**; 43 plants per sex×control/treatment before noted exclusions; 20 hand-cross supplemented females | No. visiting insects rises (reported χ²₁=6.04, p=.014) and total visits rises (χ²₁=6.91, p<.01), but no detected difference in paternal reproductive success (χ²₁=.87, p=.35) or maternal success (χ²₁=1.23, p=.27). Authors additionally report greater maternal genetic mate diversity. | Failure to detect an RS response ≠ causal zero; visitation not necessarily compatible pollen transfer. No floral dry-mass measurement. Authors already made this finding. |
| Karron & Mitchell 2012, *Mimulus ringens*, doi:10.1093/aob/mcr193 | Actual simultaneous flowers manipulated (2,4,8,16); paternity 1,935 offspring | Published total paternal siring scaling exponent **0.46** vs total maternal seeds **0.98**. An 8× display predicts ~2.60× paternal siring and ~7.67× maternal seeds from authors' fits. Extra real flowers are not identical to added *dummy* display. | Species- and conditions-specific exponents; no measured floral construction cost, lifetime future cost. Other taxa can have near-linear male gains. |
| Harder & Hamnett-Mount 2026, *Brassica rapa*, doi:10.1002/ajb2.70204 | Random 25/50/75/100% flowers receive abundant artificial cross pollen; 86 plants Exp1; 46 plants Exp2 | Lower treatment pollen fraction limits fruits, while high-pollen treatments exhibit maturity capacity limits, particularly for later opening flowers (mean fruiting window ~9 vs ~5 d); in **sampled developing fruits** 61% of ovules fertilized and among fertilized ovules 96.3% developed. | The article and original Eq 1 already demonstrated co-limitation. Floral dry mass and its effect were **not** randomized. Mechanism of older fruits inhibiting younger fruit (resource vs hormonal carpic dominance) unresolved. |

Source links:
- Silene randomized floral signal and selection: https://pmc.ncbi.nlm.nih.gov/articles/PMC9709571/
- Mimulus paternity: https://pmc.ncbi.nlm.nih.gov/articles/PMC3278290/
- Brassica fruit-fertilization-resource sequence: https://pmc.ncbi.nlm.nih.gov/articles/PMC13206284/
- Counterexample to universal diminishing male gain: *Sagittaria latifolia* male mating array experiment, linear siring gain; doi:10.1111/j.1558-5646.2011.01344.x
- Broader lifetime cost definition: Dorken et al. 2025, doi:10.1111/nph.70166. Current-year fruit output ≠ lifetime net fitness; the review did not find a simple universally greater reproductive cost in dioecious females.

### Original source statistical mismatch audit: a genuine publication-text caution
Directly recalculate survival probabilities for reported χ²₁ likelihood-ratio tests using `erfc(sqrt(χ²/2))`, without access to underlying statistical model output. In Barbot et al. (2022) online article prose, **three** displayed combinations cannot simultaneously be correct:
- Real flower count → any recorded visit: χ²₁=**3.25**, printed **p=.71**, but chi-squared(1) identity requires **p=.07142**.
- Floral-size PC1 → total visits: χ²₁=**6.44**, printed **p=.11**, but χ²₁ requires **p=.01116**.
- Added dummy visual display → female mating success: χ²₁=**4.97**, printed **p=.26**, but χ²₁ requires **p=.02579**. This aligns with the authors' prose describing a significant maternal mating-partner increase without increased maternal reproductive success.

No raw models were re-run and no editorial erratum was confirmed; **likely decimal-point transcription errors, not verified corrections**. The intended *core experiment* contrast uses the consistent reported p-values: increased pollinator counts but no detected reproductive-success benefit. Reproducible script: `scripts/audit_silene_dummy_flower_chisquare.py`.

## 2. Original 2026 Brassica Eq 1 plus an explicit costly-flower extension
Let F=total flowers, B=reproductive resource budget (not directly estimated), c_b=resources spent on constructing/maintaining each flower, p=**fraction of real flowers receiving effective female fertilization opportunity** (not observed insect landing rate), o=viable ovules per flower, z=effective ovule fertilization in a pollinated flower, c_f=fruit-development cost, c_s=cost per fertilized seed. The original 2026 paper holds available development resources D separate from flower costs and gives max fruit capacity D/(c_f+z o c_s). This branch's clearly labeled **scenario extension** subtracts floral expense first: D(F)=B−c_b F. Then, assuming fixed p, o, z and no extra trait variation,

`S_f(F) = z o min[p F, max(0,B−c_b F)/(c_f+z o c_s)]`

Threshold at which increased female floral investment changes from flower/arrival/fertilization limited to development-capacity limited:

`F* = B/[c_b + p(c_f + z o c_s)]`, for p>0.

Below F*: `dS_f/dF = z o p > 0`; above F*: `dS_f/dF = -(z o c_b)/(c_f+z o c_s) <0` when c_b>0. For c_b=0, success **plateaus** above F*, rather than falling. If p=0, female success is zero for all F, and an internal F* has no biological meaning. The kink assumes homogeneous flowers and resource sharing; **later-flower developmental order and hormones are not modeled**.

The separate *male* component is `M(F)=q_m a_m F^β_m`, where q_m is independent male pollen-export opportunity, a_m an arbitrary fitness scale, and β_m=0.46 is borrowed solely as a *Mimulus* sensitivity illustration. The unified comparison `S_f + w M` is in **arbitrary units** because offspring genetics, competition, and paternal/maternal value have not been calibrated. Do not claim an empirically inferred sex allocation or optimum.

### Completed **hypothetical** sensitivity, not empirical counts
Illustrative settings only: B=120 abstract units, c_b=.4, c_f=2, c_s=.5, o=10, z=.61 (borrowed conditional *Brassica* ovule estimate, not a Daphne whole-plant probability). No real population corresponds to these budget terms.

| Constant p assumed | Continuous F* | Best whole-number flower F for **current female output only** |
| ---: | ---: | ---: |
| 0.10 | 132.60 | 133 |
| 0.25 | 72.18 | 72 |
| 0.40 | 49.59 | 50 |
| 0.75 | 28.66 | 29 |
| 1.00 | 22.02 | 22 |

**Counterintuitive conditional prediction:** when the effective female pollination fraction rises, *fewer* flowers can be optimal under a fixed total budget, because each flower now requires a greater expected maturation expense. This is an algebraic model property, **not** evidence that improving pollination evolutionarily reduces flower numbers. More male benefit can shift the optimum upward (at p=.4, arbitrary combined weight w=3 shifts the optimum F from 50 to 128, with female output falling; relative weights wholly hypothetical).

### Equal average pollination is not equal reproductive opportunity
For whole-season/whole-year pollen fraction states with **the same mean p=.40**, the model's optimal integer flower F is:

| Year states and probabilities | Mean p | Model F optimum | Expected female seed-model units at optimum |
| --- | ---: | ---: | ---: |
| p=.40 every year | .40 | **50** | 120.79 |
| p=.05 or .75, each 50% | .40 | **29** | 69.89 |
| p=.05 (25%); .30 (50%); .95 (25%) | .40 | **63** | 90.69 |

**No universal claim that pollination unreliability must favor producing more flowers.** The optimal F can move in either direction depending on the *full distribution* of successful events, even when the mean is fixed. For a fixed F and fixed z, the `min()` capacity function is concave in p, so averaging p before imposing the cap systematically overstates or equals expected output: `E[S(F,p)] ≤ S(F,E[p])`. These are within-annual-bout budget scenarios, not a fitted Markov multi-year model; no automatic buffer by longevity or fitness of surviving seedlings is implied.

Reproduction and tests: `scripts/flower_investment_switch_frontier.py`, `.github/workflows/angiosperm-investment-switch-audit.yml`, source-check and sensitivity outputs in the GitHub Actions artifact. Float equality test corrected to `math.isclose`; subsequent run **success**:
https://github.com/zuizui0223/onishibari/actions/runs/43722276633
artifact ID **14164943006**, ZIP SHA256 `5b575867c0622e766ca077c5e9f6750d665f4d6dbaddf766b29dca4148e15e6a`.

## 3. The new biological distinction: advertisement vs reproductive substrate
The prior framework's `flower number → visitor number → reproductive success` chain was underidentified. The Silene experiment demonstrates that simply increasing visual advertisement *can change visit counts without a detected increase in sired or produced offspring*, while the same source found real flower number associated with mate acquisition. A strong focal generalizable question is therefore:

> **Do floral displays cease to deliver fitness because attractive signals no longer convert visitors to effective pollen, because sexual gamete numbers are limiting, or because developmental capacity prevents additional pollinated flowers from maturing?**

The answer need not be the same for male versus female fitness, for morphs, or for early versus late flowering. Stage gates are **not mutually exclusive**; mechanisms can co-limit a single plant and switch across environmental conditions. The ordinary flower-number regression, NDVI, or 10-min landing elasticity alone cannot locate these gates.

## 4. Falsifiable matched-treatment design; no unsourced onishibari reproductive claims
**Phase A, existing conservation-protected 2027 pilot:** Verify genetically distinct onishibari sexual morph function first. Genet/branch/bud IDs; 5 assigned pollination controls already defined, F-pollen donor positive control; stigma, pollen tube, seed formation and germination as different destructive/follow-up flower cohorts. Continue to obey existing cap: at most **50% of confirmed buds and 10 per genet** across all experimental/destructive allocations; skip underpowered genets rather than relaxing. Never infer physiological self-incompatibility from bagging alone.

**Phase B, separate low-damage advertisement mechanism:** Within confirmed morph and date, randomized small visually standardized fake-flower addition vs sham and known observation effort; measure total visits, **conspecific compatible pollen per stigma**, effective pollen export/paternity. Artificial flowers change signal, not genuine ovule or pollen availability. Record lure spectral appearance, reward absence, and pollinator response; no direct seed effects can be assumed. Unlike increasing real flowers, this is a clean intervention on advertisement.

**Phase C, resource/feedback mechanism:** Apply compatible pollen supplementation and **carefully controlled early-fruit thinning** where abundant common-garden material and permissions permit, not in small sensitive wild onishibari genets by default. Monitor tagged late flowers; rescue of later fruit set suggests prior-fruit mediated constraint but **does not distinguish** direct resource competition from hormone/carpic dominance. A further hormone manipulation plus sham would be required for that. Fertilized ovule failure and seed viability need independent sampling.

**Phase D, lifetime/net contribution:** Count offspring established next spring and future flowering/plant survival; assign paternal successes from donor population; assess benefit/cost to future reproductive value. Do not replace long-term progeny survival with seed set or adult persistence alone.

### Identifiability table for eventual experiment
- Dummy display ↑ visits, no compatible deposition change: break at *effective pollen transfer*, not necessarily seed resource.
- Compatible cross pollen ↑ embryos but not mature fruit: break *after fertilization* if embryo comparisons valid; confirm by tracking appropriate flowers.
- Fruit thinning ↑ late mature seeds: prior-fruit interference/resource feedback; cannot automatically attribute to hormones.
- Real flowers ↑ genetically sired offspring while dummy flowers do not: supply of pollen / donor fertility contributes, beyond appearance.
- Seed quantity ↑ but established offspring and next-year flowering ↓: a *lifetime* trade-off, not mere short-term reproductive improvement.

## 5. Remaining blockers and novelty ceiling
- The original Silene dataset is at Dryad doi:10.5061/dryad.0rxwdbs4d and Brassica at doi:10.5061/dryad.hx3ffbgtq. Current run **does not** download/re-fit their raw individual CSV data. Do not imply stage-specific numerical model fits on those individuals.
- The underlying Brassica cost and seed-cap model is **original authors' theory**. Our c_b-cost extension and stochastic stress are explicit exploratory bookkeeping, not evidence for a new evolutionary law.
- The 2012 Mimulus β_m=.46 does not generalize to all males: in 2011 *Sagittaria latifolia* mating-array work, male siring was instead reported approximately linear with flower production. More clades and genuinely measured male success must be compared.
- An independent 2025 review (Dorken et al. doi:10.1111/nph.70166) shows that current-year and residual reproductive value require independent measurement and that stereotypical female-vs-male reproduction costs in dioecious plants are not consistently supported.

Main claim for paper planning: **locate and compare the bottleneck distributions (advertisement, pollen transfer, fertilization, fruit development, recruitment) and sex-specific marginal gain curves**, rather than pooling incompatible visits/seed/biomass endpoints into one biased meta-analysis. The current branch now has a falsifiable design, source consistency audit and mathematical frontiers; it does not have empirical cross-clade causal estimates of the switching surfaces.
