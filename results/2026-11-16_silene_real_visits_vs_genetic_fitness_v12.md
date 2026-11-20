# What converts floral visits into genetic reproduction? Visitors, repeated contact, and pollen renewal (v0.12)

Date: 2026-11-16 JST. **NEW secondary analysis of original 2022 individual plant data completed**; separate tagged-bee 2007 literature identified as key prior work. No new Daphne field observations, and no original individual-tagged bees were measured in the 2022 source.

## Central estimands

The three quantitatively different metrics must never be interchanged:
1. Observed number of visiting pollinators ("nbVis" in original FLODIS README; individual identities are NOT tracked in this table);
2. Total observed floral visits ("nbVisTot"; multiple contacts can accumulate);
3. Effective movement of viable and compatible pollen between genetically independent plants, reflected in genetic mating and viable maternal/paternal descendants.

The first two are **highly correlated** in original FLODIS data and do **not distinguish individually tagged fresh bees from same-bee repeat arrivals**. The 2026 Joffard bee-sequence assay measures **same identified bee's sequential stigma deposition**, but is a **separate experiment** and cannot be multiplied into a mother-level genetic fitness pathway.

## Historical novelty check: tagged bees in Japanese Cirsium are already published

Makino, Ohashi & Sakai (2007, original online 2006), *Functional Ecology* DOI https://doi.org/10.1111/j.1365-2435.2006.01211.x, studied 27 individually observed Cirsium purpuratum plants for 2 h each, with individually tagged Bombus diversus bees. They directly split all recorded plant visits into **numbers of different visitors** versus **frequency of individual returns**, further categorizing occasional (1-3 times) and regular (4+ times) visitors. Larger floral displays correlated with greater density of both occasional and regular bees, whereas surrounding local flower density increased occasional visitors but was negatively associated with return rates of regular visitors. This DIRECTLY preempts any claim that separating visitor density vs revisit rate is our newly discovered principle. The original paper marked all bees possible but required regression-based imputation for unmarked visits; n=27 plant systems, not a multi-species universal result. It did not follow actual genetic mothered/sired offspring or compatibility and fruit resource capacity on the SAME focal individual. That remains the unresolved opportunity.

Ohashi & Yahara (1998), DOI listed in primary abstract https://pubmed.ncbi.nlm.nih.gov/21684905/ , showed bumblebee visits per Cirsium flower head could stay nearly constant as display changes, despite sublinear plant-level visitation. Mitchell et al. (2004), *Functional Ecology* DOI https://doi.org/10.1111/j.1365-2435.2004.00812.x, experimentally manipulated Mimulus floral displays of 2/4/8/16 true flowers and observed more within-plant sequential probes for larger displays, but still more than 30% of 16-flower-plant bouts ended after one probe. Karron et al. (2009), DOI https://doi.org/10.1093/aob/mcp030, combined visit order with paternity in Mimulus and found among-flower selfing variation attributable to geitonogamous movement. These cases show why visual display, new arrival, within-genet probing and realized mating are separable axes. No raw multi-species dataset linking all axes is claimed here.

## Original 2022 Silene dioica plant-level record reanalysis

Original paper Barbot et al. (2022) *Proceedings of the Royal Society B* DOI https://doi.org/10.1098/rspb.2022.1987. Original Dryad published dataset DOI https://doi.org/10.5061/dryad.0rxwdbs4d, published version 251371, member dataID_FloDis.txt with uncompressed SHA256 1d245bb2dd6222bdaa05b0e23521eafd286cfa533b3826f8a5d8a65f715d3f9a.

Original 192 plant trait rows include n=43 each FC, FT, MC, MT randomized fake-display vs unchanged-display control, plus n=20 separate FS hand-pollen females. This analysis compares **only 172 fake-display and its control original individuals**. Four of 172 lack both observed visitor count and visit count, leaving 168. nbVis <= nbVisTot for all 168; 41 observations have zero visitors and zero visits. Exactly 85/168 have total visits exceeding counted pollinators. However original README provides no tag-level identifier of individual bees or marked bout re-entry, so this arithmetic is **not directly observed independent new-bee arrivals**.

### Matching 4 genetic fitness outcomes in the SAME source plant

Four endpoints were fitted on exactly common plant-complete-case rows across each set of four models, using HC3 robust standard errors and exact leave-one-plant-out PRESS prediction (no seed-level pseudoreplication). All responses log(1 + source genetic or viable seed index). All models control fake-display treatment indicator, log(1 + real flowers), and log(1 + source ovules/flower for female or pollen grains/flower for male). These adjustment covariates are **observational and potentially post-treatment**; they do NOT create causal randomization of visit composition.

Candidate models per endpoint:
- baseline: real floral display + gamete supply + fake-display assignment;
- baseline + total observed visit count;
- baseline + number of observed visiting pollinators;
- baseline + both observed visitation measures.

| Outcome | independent original plants | Baseline LOO MSE | +total visits LOO MSE | +visitor number LOO MSE | +both LOO MSE |
|---|---:|---:|---:|---:|---:|
| Female viable maternal seed index | 83 | .20363 | .20726 | .20791 | .20922 |
| Female genetic mate number (8-offspring-standardized) | 82 | .03129 | .03186 | .03223 | .03217 |
| Male genetic total siring MEMM estimate | 84 | .40618 | .39285 | .39106 | .39924 |
| Male genetic mate number | 83 | .21267 | .20454 | .20447 | .20810 |

For **female total viable seed index**, neither visitors nor total visits adds a clear slope after source real flower and ovule abundance adjustment:
- ln(1 + total visits) beta **+.0420**, HC3 95%CI **[−.1144,+.1985]**, p=.599.
- ln(1 + observed visitor count) beta **+.0314 [−.1579,+.2207]**, p=.745.
- With both added, LOO MSE is WORSE than the simpler baseline model; baseline .20363 vs joint .20922. This does not prove zero response or pollen sufficiency by itself; original hand-pollen control in 2022 paper and v0.8 provides its own separate diagnostic.

For **male genetic total siring**, both separately associated with higher source genetic success:
- ln(1 + total visits) alone beta **+.19397 [+.00799,+.37994]**, p=.0409, conditional on treatment + real flower count + pollen grains/flower.
- ln(1 + observed visitor count) alone beta **+.25610 [+.04763,+.46456]**, p=.0160.
- BUT after including both source visitor metrics, their individual coefficients are uncertain: total visits beta **+.11767 [−.20763,+.44298]**, p=.478; observed pollinator number **+.13065 [−.24174,+.50305]**, p=.492.
- Female reproductive/male mating success is NOT a unit that can be substituted for deposited stigma pollen or lifetime siring.

For **male distinct genetic mates**, separate visits beta +.14977 [.00613,.29341] p=.0410, separate observed visitor number beta +.19232 [.03270,.35194] p=.0182, but both conditional coefficients uncertain. Female genetic mate predictor visitor count alone beta +.01263 [−.07407,+.09933], p=.775.

### Why the key individual visitor effect remains unidentified

Observed visitor number and total visits are extremely collinear for **females** (Spearman rho **.983**) and materially for **males** (rho **.835**); neither measure contains tag IDs for fresh unique individuals. In male siring, separate visitor and visit models have similar LOO errors (.39106 vs .39285), whereas including both worsens to .39924. A paired plant bootstrap of LOO MSE difference "visits-only minus joint" yields **−.00640, 95%CI [−.01892,+.00438]**, so extra joint complexity offers no clear predictive benefit. This is **failure to statistically distinguish contact composition**, not empirical evidence that fresh arrivals have zero additional benefit.

Similarly for females, source genetic mate outcome has n=82 plants and independent visitor metrics contribute negligible prediction after treatment, floral number and ovules/flower. The predictor group includes genuine natural flower size/number, not the experimental fake-only flower manipulation.

### Nonindependence and modelling precautions

- n=172 plant records is not n=1,721 offspring assignments; use female mother plant / male source MEMM genetically corrected per-plant endpoint, not raw offspring counts.
- Source observable distinct "visitors" may count repeated observations of a bee: nbVis is not a bee identifier. Do not label nbVisTot - nbVis a count of verified revisits.
- Read 2022 original v0.8 authors' mixed model alongside our simplified exploratory HC3 OLS: magnitudes/p-values are not exactly identical.
- The four endpoints are correlated observations in a single species with same randomized display treatment; the conditional associations are NOT separate species or phylogenetic replication.
- LOO CV uses *the same source plants*, not a fully independent garden or a different year; split-by-plant avoids pseudoreplication but does not test other locations.
- Do not infer all visible landing contacts deliver genetically compatible pollen, seed fertilization or successful future recruitment.
- The source "nbVis" field and the source "nbVisTot" field are not independent tagged-bee arrival data; our expected new-bout contrast is not available in this archival source.

## Connection to the independent 2026 pollen-carryover experiment

Joffard et al. 2026 AJB https://doi.org/10.1002/ajb2.70276 and original https://doi.org/10.57745/CF1VNW directly identify 25 bumblebees for 186 stigma visits, finding progressive loss of deposited grains and increased risk of zero grains on successive probes, replicated in 25-held-out-bee fits and after controlling 51 recipients in repo v0.11. It **does not measure** the 2022 identified bee-arrivals or link these exact sequenced bees to genetic mothered/sired offspring.

Together, the studies yield a falsifiable but **not yet experimentally established** prediction: net reproductive gain depends on independent pollen-loaded visitor **recruitment**, not merely repeated contacts to the same flowering genet. Because Makino et al. already studied tagged revisit behavior, novelty must be the **downstream genetic fitness conversion and capacity/temporal filtering** in the same marked individuals, not visitation decomposition itself.

### Minimal biological next assay in focal Daphne pseudomezereum

Mark unique genets and individual visitor identities without assuming morphological sex function or physiological SI. For each **2-hour plant observation bout**, record:
1. Unique visiting bee identity (tag/video tracking where possible) and number of plant returns by that bee;
2. Per-plant flowers open and viable pollen output (dry mass separately); within-genet flower order per bout;
3. Conspecific compatible pollen deposited on genotyped female stigmas, later viable seed and paternal genotype;
4. Fruit resource-development history and subsequent survival;
5. Surrounding flower density, which Makino et al. identified as an important predictor of visitor turnover/revisitation.

With wildplant intervention cap no more than 50% of verified buds and no more than 10 buds/genet, and observational tagging as feasible, avoid assuming a large population is available. Distinct arrival-versus-repeat effects require actual unique pollinator tracking in exactly the same target plant and time window, not a proxy.

## Reproduction status

Code: scripts/silene2022_visitor_vs_visit_mating_conversion_v12.py
Workflow: .github/workflows/angiosperm-silene-visitor-revisitation-v12.yml
Success: https://github.com/zuizui0223/onishibari/actions/runs/43828178689
Source artifact ID 14206869498, artifact ZIP SHA256 dd458f26b03cd051bc0a61e00b978a0833ae9f0c879a93f483defdad6f539519. Source records excluded from outputs, only aggregate estimates. v0.12 will have its own stage-evidence ledger source-value verifier.

**Net status:** strong meaningful negative result for identifying fresh visitor contributions from standard visitation observations, alongside an observational male-genetic success signal. The hypothesis is plausible, not causally proven, and its novelty is narrow because prior individual-tagged studies exist.
