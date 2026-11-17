# Pre-registered next phase: Eastern Asian Daphne cytonuclear discordance and sexual-system divergence (draft v0.1)

**Date:** 2026-11-15. **Status:** prospective design; no new specimens or independent nuclear loci sequenced.
**Primary question:** Is Korean C more closely related genomically to Japanese P or Kamchatka/Hokkaido K? What process explains its nuclear-rDNA ITS genealogy opposing the plastome genealogy?
**Secondary question:** Is sex expression/flower size dimorphism conserved or variable after lineages diverge?

## Existing evidence and limits
- Source-marker analyses validated under expanded outgroups, taxon sampling and AU alternative-tree tests; results in `results/2026-11-15_topology_robustness_extended.md`.
- Morphology/rhizomes shared by C and K, chloroplast C-K; ITS P-C; see `results/2026-11-15_morphology_cytonuclear_interpretation.md`.
- Korean 2023 chloroplast ON244034 and 2025 nuclear ITS PX056848 use *different* individuals.
- Korean public SRA SRR19117239 is 99.96% plastome-mapped and 0 ITS-mapped; cannot independently check nuclear genealogy.
- False sample labeling risk: GenBank taxon IDs of D. jezoensis often return under D. pseudomezereum.
- Some Honshu/Shikoku mountain plants previously named var. koreana likely do **not** correspond to Korean type var. koreana.

## First low-cost action
Request from the corresponding author of Han et al. (2025), Jung-Hyun Lee (contact in article: quercus@jnu.ac.kr), original **unfiltered** 2×150bp Illumina total-DNA FASTQs and specimen IDs for
- Japanese P `Lee-Dpp230303` (2023 Mie; assembled plastome PV752465),
- Kamchatka K `K005971094` (2003 herbarium; assembled plastome PV752464),
- Korean C `Byeon-Dpk250701` (2025 Pyeongchang; ITS PX056848), including any unfiltered WGS reads if available.
Additionally ask whether unfiltered original Korean 2023 `Lee 7891` Gangneung reads underlying ON244034 are retained, since public SRR19117239 lacks nuclear sequence.
**No correspondence is sent by this research branch.** These identifiers originate from the published article and public repositories.
If permission to reuse raw read data is granted, pilot recovery of Angiosperms353 loci and rDNA copies before any field sequencing.

## Provisional prospective sampling (adjust according to available access)
- P `D. pseudomezereum` (Japan warm-temperate): 3 populations × 6 genets, if multiple morphs available.
- C `D. kamtschatica var. koreana` (Korea): 3 populations × 6 genets, include geographical contact gradient if supported by collections.
- K `D. kamtschatica var. kamtschatica` (Hokkaido/Kamchatka): 3 populations × 6 genets (at least Hokkaido + one Kamchatka reference if feasible).
- Unresolved Honshu/Shikoku montane OTU (not automatically labeled Korean C): 2 populations × 6 genets.
- `D. jezoensis` reference: 2 populations × 6 genets; distinguish females and male-biased hermaphrodites.
Total planned = 13 population units × 6 genets = **78 genets**, core mechanism P/C/K = **54**. These counts are a design target, not observations or guaranteed sampling availability. Sample different *genets* where rhizomatous clones occur (avoid adjacent ramets as independent samples); verify genotype and de-duplicate if necessary.

## Molecular data for every voucher
- Paired plastid and independent low-copy nuclear sequences **from the same individual**. Use Angiosperms353 target capture (not requiring a reference genome), backed by a pilot on each focal lineage. The 2026 European Daphne Angiosperms353 study doi:10.1093/aob/mcag358 demonstrates feasibility in genus but not coverage for P/C/K.
- Aim for ≥100 recovered orthologous, phased where possible, low-copy independent nuclear loci; the exact acceptable number needs simulation/pilot. Adjacent SNPs from one locus are NOT independent genealogies.
- Genome skimming chloroplast from off-target reads; track organelle haplotype per specimen, sequence length and IR handling.
- Flow cytometry chromosome ploidy/aneuploidy controls; an ancient allopolyploid history cannot be inferred from ITS and cp alone.
- Voucher, sampling locality, phenophase, DNA source, sequencing batch, collector, taxonomic determination; deposit FASTQs with sample-to-voucher maps.
- Pilot checks: allelic dropout and capture bias, reciprocal divergence, paralogy and ITS copy variation, PCR artifacts, clonality, and chloroplast heteroplasmy where appropriate.

## Inference gates set before viewing multilocus results
1. Estimate **nuclear species tree** under MSC (ASTRAL or equivalent, with gene support/concordance, loci treated as independent where appropriate), include uncertainty and alternate population partitions.
2. Check whether discordant nuclear gene-tree alternatives occur at comparable frequency after accounting for coalescent branch length and tree-estimation error. Symmetric discordance is **compatible**, not proof, with ILS; ancestral structure and gene flow complicate simple expectations.
3. Only after fixing rooted nuclear species-tree hypotheses use Dsuite ABBA-BABA/f-branch and migration/network MSC comparisons (SNaQ, PhyloNet, or equivalent) with suitable outgroups. Assess uncertainty, ghost/unsampled donors and admixture direction.
4. Compare cp-vs-nuclear haplotypes among same individuals and geography. Genuine organelle capture can leave little current nuclear admixture; non-significant D alone cannot reject it.
5. **Conditional conclusion discipline:** if majority low-copy nuclear loci favor C–K, investigate rDNA nuclear introgression/concerted evolution vs ILS; if P–C, investigate chloroplast capture vs ILS; if unresolved, avoid assigning mechanism.
6. Do NOT estimate speciation rates from 3 focal lineages. Distinguish population genetic differentiation, phylogenetic divergence, reproductive isolation, and hybrid speciation.

## Ecological phenotype on the same voucher/population
- Flowers: calyx tube/lobe size and display per flower and per individual, pollen fertility, stigma function, flower number, sex phenotypes, sex expression shifts over years, sex-specific plant size; quantify phenological timing and lifespan context.
- Pollinators: day vs night exposure, visits per flower-hour, conspecific pollen deposition and pollen tube growth, controlled self/cross and bagging treatments; visiting rarity does not automatically equal low pollination effectiveness.
- Fitness: outcross/self seed production and viability, adult survival, seedling recruitment; test whether male-function-biased hermaphrodites have residual functional female fitness.
- Analyze sex expression as continuous **functional** measures instead of treating historically morphological 'dioecy' as a confirmed reproductive state.

## Practical stop rules
- If raw P/K/C nuclear reads remain unavailable, **report marker-level discordance** and stop before mechanism ranking. Draft a genome-enabled study rather than publishing unsupported introgression.
- If fewer than ~100 truly independent recoverable nuclear loci or phylogenetic uncertainty remains high after pilot, revise inference claims and prioritize targeted resequencing or RAD with locus-dropout diagnostics.
- Deposit code/data with immutable accession IDs and source SHA checksums, maintain no false claim of currently available complete nuclear data.


## Addendum — archive feasibility + locus target adaptivity (2026-11-15)

The initial fixed target of ≥100 nuclear loci should be treated as a **pilot benchmark, not a calibrated full-analysis target**. Under a simple 3-taxon MSC, the expected fraction of discordant gene trees is `2exp(-t)/3` for internal branch t. A two-sided binomial test of an **illustrative 2:1 imbalance among the two discordant alternatives** has ~82% power conditional on 80 observed discordant independent loci; a mean of 80 such discordants requires ~198 total loci if t=0.5 or ~327 loci if t=1.0, even before missingness or tree error. These are *assumption-based sensitivities*, not an introgression test or guarantee.

The source-verified 2026 Bonnaudet aligned nuclear dataset offers a **50-locus proxy reference FASTA** across D. sureil, D. tangutica, D. mucronata and D. laureola (4×50=200 reference sequences, median 4-way fully observed length 973.5 bp). These species do NOT include P/C/K. Use panel only for a first 20–50-locus recoverability test when original genomic reads are obtained, with paralogy/coverage QC and only then plan the full panel.

Korean Yoo et al. 2023 SRA SRR19117239 is the **sole run** indexed under BioProject PRJNA835716 in both ENA/NCBI (359.8 Mb vs 7.34 Gb initially reported). The public reads have effectively no recoverable nrITS; prioritize contacting authors for unfiltered raw 7.34 Gb and for Han et al. 2025 original P/K total-DNA FASTQs. Prepared **unsent** requests are in `docs/author_unfiltered_read_requests_unsent.md`.

Reproducible inputs/results: `results/2026-11-15_archives_proxy_panel_power.md`, and successful runs 43682778058 (archive), 43682931727 (50-locus panel), 43683047543 (conditional power).

**Go/no-go:** If no original or new independently derived P/C/K nuclear reads become accessible, stop at robust cytonuclear marker discordance. Do **not** identify ILS vs introgression from nonfocal nuclear loci or from two marker compartments.
