# Flower investment conversion gates: original Silene and Clarkia raw-data refit (v0.8)

Date: 2026-11-16 JST. Status: **original public-data secondary analysis completed**, not new field data and not a valid cross-study causal meta-analysis.

## Research question
When do increasing (a) visual floral advertisement, (b) genuine pollen/ovule-producing flowers and (c) compatible pollen delivered to flowers fail to increase genetic mates, mature seeds, effective siring or future viable offspring? These are **distinct interventions and biological denominators**.

## Sources, verification and replication
- Barbot et al. 2022, *Silene dioica*, DOI 10.1098/rspb.2022.1987, Dryad 10.5061/dryad.0rxwdbs4d published version 251371. Original groups 43 each MC/MT/FC/FT, plus 20 additional female hand-pollen controls (FS); total 192 individual trait records. Genuine flower numbers are unaffected by the *dummy-flower visual signal manipulation*. Original 556 fruit rows and 1,721 confidently assigned offspring–candidate-father rows are **not** 556 independent plants or the full 2,030 genotyped offspring.
- James & Geber 2025 *Ecological Monographs*, four *Clarkia* species, DOI 10.1002/ecm.1649, Dryad 10.5061/dryad.k3j9kd56r version 332601. Original 430 FRUIT rows, 175 genuine within-focal-plant pairs of natural-control vs supplemental-hand-pollen flowers, plus 80 single-control rows. 19 paired cases missing seed/ovule outcomes, so n=156 matched plant fruit pairs in 38 blocks.
- Joffard et al. 2026-10-13, *American Journal of Botany* DOI 10.1002/ajb2.70276: separate, independently published *S. dioica* female display manipulation assigning paternity for 1,320 genotyped seeds, reporting no detected difference in within-fruit paternal diversity or among-fruit sire profiles. This is **not** original 2022 FLODIS or evidence of unchanged total sired offspring.
- All original observations retrieved via documented **anonymous official Dryad published-version ZIP API**. Uncompressed original members checked by SHA256 previously independently inventoried; no raw plant, plot, maternal or candidate father IDs saved to GitHub. Scripts: scripts/silene_clarkia_design_audit_v08.py, scripts/silene_clarkia_original_stage_refit_v08.py; workflows angiosperm-silene-clarkia-audit-v08.yml and angiosperm-silene-clarkia-stage-v08.yml. Design audit success run https://github.com/zuizui0223/onishibari/actions/runs/43788290584 ; full source refit success https://github.com/zuizui0223/onishibari/actions/runs/43789194671 ; artifact 14191921387 (time-limited). 3,500 bootstrap draws, 5,000 permutation draws; original source plants/blocks are the independent resampling grains.

## 1. Silene: dummy floral visual signal vs actual reproductive success

Random dummy flower addition increased *apparent* floral display without increasing real flowers, an intervention on pollinator attraction rather than gamete supply. Mean genuine flower count females FC 4.395 vs FT 4.504 (n43 each); males MC 19.291 vs MT 18.752 (n43 each), no identifiable difference by random label permutation.

| Endpoint / source-grain | Fake-added minus control | Uncertainty / caution |
| --- | ---: | --- |
| Female insects observed, n42 vs 42 | **+0.571 visitors** | individual-bootstrap 95%CI **[+0.071,+1.095]**, permutation p=.033 |
| Pooled insects observed, n168, adjusted sex and natural flower count | **IRR 1.272** | HC3 Poisson 95%CI **[1.004,1.610]**, p=.0458; no sex-treatment interaction in this sensitivity |
| Female total viable seed index, 43 vs 43 | **+0.286 on log(1+seed) scale** | bootstrap 95%CI **[−0.182,+0.796]**, p=.257 |
| Male MEMM genetically adjusted siring index, 43 vs 43 | **+0.051 on log(1+male fitness index) scale** | 95%CI **[−0.259,+0.362]**, p=.760 |
| Female genotyped mating partners standardized for eight offspring/mother, n41 vs43 | **+0.486 genetic mates** | 95%CI **[−0.003,+0.986]**, permutation p=.0546 |
| Male genotyped mating partners standardized, n42 vs43 | +0.226 | 95%CI [−1.094,+1.538], p=.740 |

Authors' 2022 more detailed model reports an increase in female mating partner number, despite no detected increase in female seed output; their prose has a χ-square/p-value decimal inconsistency recorded in v0.4. Our *different* raw simple unadjusted and covariate-adjusted tests place the effect near but not clearly below .05 (adjusted female mating treatment beta +0.418, CI [−0.016,+0.853], p=.059). Do not claim the original result is disproved or that a statistically null fitness response is *exactly zero*.

In exploratory regressions controlling genuine flower count and gametes per flower:
- Female total viable seeds: fake treatment beta +.081 [−.109,+.272], p=.403 (n85), while true flower-number predictor beta **+1.506 [1.285,1.728]** and ovules per flower beta **+1.877 [1.425,2.328]**.
- Male genetically estimated siring index: fake treatment +.104 [−.167,+.375], p=.453 (n86); true flower-number beta **+.690 [.462,.918]**, pollen per flower −.063 [−.324,+.197], uncertain.
- Removing the source-reported exceptional manipulated male with 70 average real flowers (next highest 48.67) changes unadjusted male log-siring contrast from +.051 to +.021 [−.286,+.333], p=.897. The source 1,721 assigned-offspring records have 85 candidate fathers; 900 assigned MT vs 821 MC, but such **raw counts are NOT equal unbiased paternal success** due maternal sampling and genetic assignment; use authors' MEMM median estimates for main results.

**Interpretation:** random fake visual display can increase insect attraction **without a detectable proportional improvement in viable maternal seeds or fathered offspring**, whereas observational actual flower number is associated with available reproductive substrate. This does **not locate the precise stigma-compatibility/pollen-tube bottleneck**, nor prove fake flowers are fitness-neutral.

### Hand-pollen control: was female reproduction already pollen-limited?

The original 2022 Silene data include **20 additional independent female plants receiving supplemental hand pollen**. AllFruits has **556 fruit records from 105 mothers**. Using the original README definition that nbOv is *unfertilized ovules*, filled seed fraction is nbSeeds/(nbSeeds+nbOv). Fruits were averaged *within each mother* before comparing group treatments, so fruits are not independent treatment units:

| Female treatment | Distinct female mothers with source fruit measurements | Mean mother-level seed set |
| --- | ---: | ---: |
| Natural-pollen control FC | 42 | **0.6116** |
| Fake-display FT | 43 | **0.6312** |
| Extra compatible hand pollen FS | 20 | **0.5821** |

FS minus FC: **−0.0295**, 95% independent-mother bootstrap **[−0.1113,+0.0507]**, permutation p=.459. FT minus FC: **+0.0196**, 95% bootstrap **[−0.0457,+0.0861]**, p=.564. The original authors likewise reported **no evidence for pollen limitation** in the experimental garden. This adds a crucial competing explanation: if female fruits already obtain sufficient effective pollen, attracting additional visitors need not increase seed filling. A nonsignificant result does not prove complete absence of pollen limitation or that attracted visits are inherently ineffective. FS versus FC compares **different mothers**, not a within-mother pollen intervention; sample size FS n20 limits power. Fruit seed set is different from total female viable-seed production including fruit number and germination.

## 2. Clarkia: compatible pollen addition converts to seeds in three of four species

Every accepted pair represents *the same plant* and two different fruits, one control and one hand-supplemented. Treatment was assigned haphazardly, not necessarily through formal randomization. Each seed-set fraction is filled seeds per source counted or assigned ovules, from individual fruit, not whole-plant total seeds. Source 430 fruit rows -> 175 source matched pairs -> **156 with complete seed and ovule outcomes**, 38 background-design blocks. No missing records treated as 0.

| Species | n matched plant pairs | Control fraction | Supplemental fraction | **Additional seed-set percentage points** | 3,500-block-bootstrap 95%CI | Exploratory BH q (4 species) |
|---|---:|---:|---:|---:|---|---:|
| C. cylindrica | 37 | .458 | .607 | **+14.9** | [+3.3,+26.6] | .0165 |
| C. speciosa | 31 | .696 | .643 | **−5.3** | [−13.4,+1.9] | .2000 |
| C. unguiculata | 52 | .343 | .617 | **+27.4** | [+17.1,+37.5] | .0008 |
| C. xantiana | 36 | .437 | .545 | **+10.8** | [+3.1,+19.4] | .0136 |

Removing all source-assigned (rather than directly counted) ovule totals gives +14.2, −5.3, +26.7 and +11.2 percentage points respectively, broadly preserving the pattern. *C. speciosa* has higher natural seed set (69.6%), but this is not proof of saturation or true compatibility physiology. Four closely related species are NOT four independent evolutionary origins.

### Is flower count itself the moderator of pollen limitation?

Exploratory OLS of **within-plant seed-set improvement** on mean-centered log(1 + real plant flower count), clustered by the source experimental block, finds no clear within-species slope:

| Species | beta per ln(1+flowers) | 95%CI | Habitat-adjusted beta | 95%CI |
|---|---:|---|---:|---|
| C. cylindrica | −.0075 | [−.1047,+.0896] | −.0143 | [−.1108,+.0822] |
| C. speciosa | +.0649 | [−.0013,+.1310] | +.0547 | [−.0194,+.1287] |
| C. unguiculata | −.0196 | [−.1252,+.0861] | −.0429 | [−.1456,+.0598] |
| C. xantiana | +.0349 | [−.0557,+.1255] | +.0208 | [−.0662,+.1077] |

Habitat-adjusted sensitivity includes source neighborhood background Clarkia dry biomass, grass and forb numbers. This **does not prove zero flower-count effect**; source pairs select plants with ≥2 reproductive flowers, and whole-genet total seed/paternity and *true floral dry mass* are unavailable. But simple "more flowers always eliminates pollen limitation" is **not supported** across these four comparable within-plant source designs.

## 3. How the cumulative empirical results modify the thesis

The previous v0.6/v0.7 source verified **Brassica** (late artificially cross-pollinated flowers fail fruit maturation), **Veratrum** (flower-local seed benefits of adding pollen do not equal whole-plant seed gains), and **Silene latifolia** (2004 flower-count/seed output positive association, exploratory final 2007 survival association negative but not robust to timepoint scan/missingness). This new v0.8 adds:
1. **Advertisement-to-functional reproduction decoupling** using randomized fake flowers on Silene.
2. **Pollen-to-seed response heterogeneity** across four Clarkia species with paired fruit-level supplementation.
3. **No simple effect of REAL flower count on the size of pollen limitation** at current Clarkia sample size; a species/partner-context gate must be considered instead of assuming flower-count saturation is universally deterministic.
4. **Paternal genetic diversity ≠ paternal total offspring**, supported by the separate Oct 2026 Silene experiment (DOI 10.1002/ajb2.70276).

**Not estimated:** a global common transition probability from 1 mg of constructed flower biomass to viable maternal or paternal offspring across species. Outcomes/experiments are not algebraically chainable or valid as one fitness meta-effect. The original authors already reported most of these mechanisms; our contribution is a directly source-checked comparable stage evidence map, new bounded descriptive analyses, and explicit inference boundaries.

## 4. Concrete next experiment and general test

For rare wild *Daphne pseudomezereum* maintain already prespecified conservation limits (≤50% of confirmed buds and ≤10/genet, with independent-genet identification and positive hand-cross controls). **Verify the actual function of the retained pistil first; do not assert physiological self-incompatibility from floral appearance or bags alone.** Separate (i) sham versus fake visual advertisement at identical genuine flowers and pollen supply, (ii) compatible cross pollen supplementation and genetic pollen-tube/embryo fate, (iii) flower and floral dry-mass allocation versus verified male paternity and female viable descendants, and (iv) current vs next-year survival and flowering. Only a matched multi-year genotype and multiple independently evolved clades can distinguish fitness transition thresholds from evolutionary legacy or phylogenetic constraints.

Source IDs and original geolocations are excluded from this public report; this is a reproducible proposal and secondary analysis, not a claimed new Daphne field study.
