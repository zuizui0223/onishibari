# Onishibari public records: actual acquisition and limits (2026-11-14)

**Status: public-data snapshot obtained successfully from GBIF and iNaturalist APIs. These are real record counts, NOT floral sex/fruit-set/pollinator observations.**

## Verifiable GBIF occurrence counts
Run [#43528264121](https://github.com/zuizui0223/onishibari/actions/runs/43528264121), API query 2026-11-14 08:35 JST. Species identifiers resolved by GBIF names matching with exact species-level accepted names. Source script `analysis/gbif_occurrence_audit.py`; full snapshot (and original media links) retained as [artifact](https://github.com/zuizui0223/onishibari/actions/runs/43528264121/artifacts/14085194012) temporarily.

| Species | GBIF key | Japan all present | Japan georeferenced | Japan 2016–2026 georeferenced photo records | Human observation subset |
|---|---:|---:|---:|---:|---:|
| *D. pseudomezereum* | 5524168 | 649 | 362 | **26** | 25 |
| *D. jezoensis* | 5523246 | 428 | 246 | **32** | 32 |
| *D. kiusiana* | 5523189 | 650 | 324 | **16** | 15 |

Important: 26 photographs **are GBIF occurrence records with media**, not 26 independent plants; 25/26 are iNaturalist-syndicated. The GBIF API includes other non-GBIF providers; data coverage and taxonomy are imperfect.

### Observed month in the restricted photo subset (counts of records, not flowering)
| Month | *D. pseudomezereum* | *D. jezoensis* | *D. kiusiana* |
|---|---:|---:|---:|
| Jan | 4 | 1 | 2 |
| Feb | 8 | 0 | 4 |
| Mar | 5 | 2 | 5 |
| Apr | 2 | 24 | 2 |
| May | 2 | 3 | 0 |
| Jun | 2 | 0 | 1 |
| Jul | 0 | 2 | 0 |
| Aug | 0 | 0 | 0 |
| Sep | 0 | 0 | 0 |
| Oct | 1 | 0 | 0 |
| Nov | 0 | 0 | 1 |
| Dec | 2 | 0 | 1 |

Within this biased sample, January–April accounted for 19/26 (73.1%) *D. pseudomezereum* photo records; for *D. jezoensis* April accounted for 24/32 (75%). This may reflect species range/climate, different observation effort, uploaded photo selection, phenology, or all of them. **No flowering-date shift or ecological association has been established.**

## iNaturalist target species audit
[Run #43528570859](https://github.com/zuizui0223/onishibari/actions/runs/43528570859), source `analysis/inat_onishibari_audit.py`. Japan place ID 6737, *D. pseudomezereum* taxon 447647, observed 2016-01-01 through 2026-11-14.

| Count category | Observations |
|---|---:|
| All observation events (any grade) | **61** |
| With photo(s), any grade | **61** |
| Research-grade photos | **28** |
| Research-grade photos explicitly annotated flowering | **4** |
| Research-grade photos explicitly annotated fruit/seed | **4** |
| Research-grade photos with male sex annotation | **1** |
| Research-grade photos with female sex annotation | **0** |

Annotated sex is community metadata, **not proof of functioning female/male flowers**. Zero female annotations means **missing annotation**, not missing females. Because 25 GBIF target photo records are sourced from iNaturalist, **do not add GBIF and iNaturalist record counts**. Restriction to research-grade and exact observation dates explains different totals.

Research-grade iNaturalist photo months were Jan 4, Feb 9, Mar 6, Apr 2, May 2, Jun 2, Oct 1, Dec 2. The observation date is not a validated date of flowering.

## New decision (based on retrieved data)
- **GBIF plus iNaturalist observations alone cannot currently support a robust macro-level model of the floral sex ratio.** There are too few sex-annotated records and uncertain reuse of the same individual plant.
- Photo triage is still feasible as a pilot: 28 research-grade observation events, up to several photos per event. Classify flowering/fruiting, visible anthers/pistils, and reliability separately. Do not infer female function solely from the presence of a pistil.
- A large-scale study should separate macro-scale **range, climatic niche, observation-season availability** (potentially hundreds of occurrence records) from population-level **sex morph and compatibility** (requires targeted field measurements).
- Priority sequence: (a) screen photo visibility/sex evidence, (b) field survey local sex ratios/flower number, (c) randomized matched pollen-mix experiments, (d) pollen tube & embryo assays, (e) seasonality/pollinator comparison.
- Can compare the published mating-system dichotomy: *D. jezoensis* exhibits low hermaphrodite seed fertility and inbreeding depression; in *D. kiusiana* an eastern lineage displays smaller flowers with elevated selfing/homozygosity (https://pmc.ncbi.nlm.nih.gov/articles/PMC10184445/). **Never extrapolate these mechanisms to onishibari without reproductive experiments.**

## Sources
- GBIF Occurrence API https://techdocs.gbif.org/en/openapi/v1/occurrence
- iNaturalist observed records interface https://www.inaturalist.org/observations?taxon_id=447647&place_id=6737
- Kew POWO accepted taxon *Daphne pseudomezereum* https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:831337-1
