# v0.16 | Source-original Daphne jezoensis day/noon/night pollination: unit and identifiability gates

Date: 2026-11-16 JST. **Prospective source audit / reanalysis specification only.** No new experimental data, onishibari measurements or verified effect estimate are claimed here.

## Reason to extend the broad reproductive-stage analysis
The 228-species pollen-fate data measure anther removal and conspecific-stigma delivery, not fruit set, compatible pollen, outcross siring, or a diel mechanism. A separate, source-original same-species study could test whether *access window*, flower exposure and mature reproductive output can be aligned. It does **not** automatically provide the missing same-flower full conversion chain.

## Source to inspect
- Public Zenodo record: https://zenodo.org/records/8253921 (verify title/authors/version/license before citing as a formal original-data outcome).
- Candidate supplied tables: Nopporo_daynoon.csv, Nopporo_control.csv, Nopporo_outcross.csv; paternal assignment table(s) if present.
- Inspect archive version, every original SHA256, table columns, type, missingness, recorded design, intervention timing, and sampling denominators. Do not assume a filename is a biologically defined treatment group.
- The 216-row candidate day/noon table count is a **raw-file inventory count**, not 216 independent individuals or experimental replications; it must be rechecked before downstream analysis.

## Stop/go hierarchy
1. **Identifier integrity**: distinguish site, year, reproductive genet, ramet, inflorescence, flower, pollinator observation, and offspring; audit uniqueness and whether multiple rows of a unit are repeated or independent. No invented pairing.
2. **Treatment semantics**: confirm exact timing of bagging, day/noon/night exposure, whether all treatments ran simultaneously on the same genet, and whether control/outcross files have matching flower cohort, date and plant IDs. If no reliable join, fit independent panels and decline causal mediation.
3. **Valid denominator**: identify female/receptive flowers versus all male+female blossoms, monitoring period and any sampled flowers; missing outcome != observed zero. For binomial fruit set use fruits / verified receptive flowers (not visits or inflorescences).
4. **Within-genet comparison only if identifiable**: paired contrasts use treatment-defined exposures in the same genet×year, each with appropriate observed denominators. Resample/permutate matched genets, not records. If treatments are at distinct genets/site, use stratified independent comparisons instead.
5. **Independent outcomes**: compare the exposure/time contrast with natural-open versus hand-outcross to bound pollen limitation where positive controls are valid. Source assigned paternity (if verified) is a different outcome; do not multiply correlations across unlinked flowers to claim a mechanistic pipeline.

## Estimands, conditional on source coding
- Primary maternal outcome: paired difference in mature fruit probability, p(day-access)-p(night-access), using the actual experiment's operational exposure definitions and no extrapolation to true 24-h day and night.
- Secondary maternal outcome: 1 - p(open)/p(verified hand-outcross) only when cohorts, morph and dates are compatible; note resource constraints mean this is an index, not a direct estimate of all pollen limitation.
- If parental genotypes are present: outcross-siring fraction by access window and genotype identity, with maternal families as clusters. No inference from visible pollen alone.
- Compare effect **heterogeneity** by verified sex morph, local flower density, bloom stage and season only when the original data carry these variables with adequate replication.
- Mixed or beta-binomial model with genet and year if source units support it. Include observed opportunity/effort, source site and treatment order. Prefer exact matched contrast and cluster bootstrap for primary effect; no p-values if only a small number of independent experimental blocks.

## General claim to test, not assume
A conspicuous or larger flower can increase encounter opportunity, but a novel visitor can fail to deliver compatible pollen, and a pollen-sufficient fruit can still fail to mature. Distinct sex morphs might depend on different stages; large flowers and retained pistils do not prove any of these payoffs in *Daphne pseudomezereum*. The independently published *D. jezoensis* data are an evolutionary comparator, never a measurement of the focal species.

## Non-negotiable quality criteria
- No release of precise locations or individual genet identifiers.
- Published source results attributed to their original authors; a secondary fit is explicitly labeled secondary.
- Source hashes and group counts verified in CI before any new numerical claim is added to the cumulative evidence ledger.
- Keep separate draft PR; do not merge to main.
