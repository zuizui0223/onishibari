# Conditional maintenance of residual female function: toy threshold model

**THEORETICAL SCENARIOS ONLY. No *D. pseudomezereum* mating-system data, no estimated colonization frequency, no true sex-system fitness or equilibrium.**

**REVISION 2026-11-14:** The older 10.9% threshold assumes mature selfed seeds are produced whenever compatible mates are absent; it does NOT model whether effective pollinators can transfer self pollen. This can fail biologically even in self-compatible plants. See [the two-channel model](../analysis/two_channel_pollen_mate_sensitivity.py) for a mate × pollinator sensitivity test; neither model is fitted to Onishibari.

### Only empirical input: viability after a seed already exists, from a different species

- *Daphne jezoensis* Shibata (2020): outcrossed offspring retained 31/40 = 0.775 and selfed 6/19 = 0.316 through the three-year stage.
- Selfed survival relative to outcrossed survival **r = 0.40747**, implying inbreeding depression **δ = 0.59253** at this measured horizon.
- A diploid selfed child receives two parental nuclear gene copies; a child produced via one outcross parental route receives one. The crude relative return of a selfed versus outcrossed seed is **2r = 0.81494** (below 1 when δ > 0.5).
- These are **post-seed** values, not rates of autonomous selfing, fruiting, pollen transfer, complete lifetime fitness or the target species' parameters.

### Two-state conditional budget (simplified)

- p: fraction of reproductive opportunities under severe absence of compatible mates (not measured and NOT equivalent to geographic isolation distance alone).
- S: expected mature selfed seeds produced **per mate-limited opportunity** by a residual-female-function morph (including whether self pollen can be delivered).
- c: loss of viable outcross offspring sired in a mate-rich opportunity, comparing a residual-female morph with an otherwise identical fully male morph.
- k: extra fixed maintenance/future-fitness cost, in the same normalized gene-copy-equivalent units per reproductive opportunity.
- r: post-seed self viability relative to outcross viability, borrowed from *D. jezoensis* solely for illustration.

For nuclear genetic transmission accounting (gene-copy proxy), an assumed difference between retaining a pistil (R) and full male loss (M) is:

    ΔG = p × (2 r S) − (1 − p) × c − k

Female function is predicted to have a positive *direct return in this toy budget* for p above:

    p* = (c + k)/(2 r S + c), if p* < 1 and S > 0.

This is NOT a gynodioecy-equilibrium equation and ignores selfing and outcross female output when mates are available, pollen/seed discounting, frequency dependence, nuclear–cytoplasmic sex determination, and variation in paternal mating success. It must not be treated as an adaptation test.

### Illustrative outcomes, fixed maintenance cost k = 0 (entirely assumed)

| S: mature self seeds when mates scarce | c: male return lost when mates present | Gene-copy break-even p* | If counting surviving offspring only, p* |
|---:|---:|---:|---:|
| 0.25 | 0.02 | 8.9% | 16.4% |
| 0.25 | 0.05 | 19.7% | 32.9% |
| 0.25 | 0.10 | 32.9% | 49.5% |
| 0.50 | 0.02 | 4.7% | 8.9% |
| 0.50 | 0.05 | 10.9% | 19.7% |
| 0.50 | 0.10 | 19.7% | 32.9% |
| 1.00 | 0.02 | 2.4% | 4.7% |
| 1.00 | 0.05 | 5.8% | 10.9% |
| 1.00 | 0.10 | 10.9% | 19.7% |

**Example**: if S=0.5, c=0.05 and k=0, a residual-female morph has positive genetic-transmission return only when severe mate scarcity occurs in more than **10.9%** of opportunities. If counting only surviving descendants, the analogous threshold is **19.7%**. These are distinct units and NOT empirical rates.
If fixed cost k=0.02 in that same gene-copy scenario, p* becomes **15.3%**.
- If S=0 (no viable selfed mature seeds when mates absent), retaining this function never yields a strict positive direct benefit in this minimal model.
- If c=k=0 and S>0, any p>0 yields a positive direct return but that does not demonstrate an evolutionary constraint or persistence.

### What to measure in *D. pseudomezereum*

1. Autonomy: bagged flower self-seed formation and viability in the absence of visitors; if there is no autonomous selfing, mate absence is not necessarily a selfing opportunity.
2. Viability: genetic self versus outcross assignment, mature seeds per manipulated flower, three-year and preferably lifetime survival; estimate r with maternal/genet-level uncertainty.
3. Costs: repeated whole-plant flower production, pollen viability, actual genetic siring and next-year performance, with reproductive manipulation and a sham; histological pistil formation costs are not the same as fruit production costs.
4. Frequency of mate scarcity: mapped, synchronous flowering neighbors, actual compatible pollen deposition, pollinator activity, and colonization history; ordinary geographic distance is not a direct p measurement.
5. Sex-system genetics and frequency-dependent selection require a separate, full population-genetic model after phenotypic and genetic data exist.

### References
- Shibata (2020) PhD thesis, data source: https://doi.org/10.14943/doctoral.k14152
- Shibata et al. (2018) near-functional dioecy and no assured population selfing: https://doi.org/10.1007/s10265-017-0978-5
- Busch and Delph (2012) reproductive assurance vs automatic transmission advantage: https://pmc.ncbi.nlm.nih.gov/articles/PMC3278291/
- Broader theoretical review of mating-system transitions: https://pmc.ncbi.nlm.nih.gov/articles/PMC3278288/
