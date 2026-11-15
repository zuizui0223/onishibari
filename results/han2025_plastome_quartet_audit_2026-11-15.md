# Han et al. (2025) whole-plastome quartet source audit (observed, not species-tree inference)

2026-11-15. Study: https://doi.org/10.11110/kjpt.2025.55.3.157.
Workflow run: https://github.com/zuizui0223/onishibari/actions/runs/43677393796
Artifact ID: 14147185854, name `daphne-plastome-public-sequence-audit`, ZIP sha256 `7cca838598edfdf7c1e642e0b4f70c32627c07ae6c05f60e13eac097907a9151`.
Executable source: `scripts/daphne_plastome_audit.py` (NCBI EFetch + MAFFT).
This document records a previous successful source-sequence and site-pattern run. ML analyses were initiated in a subsequent workflow revision; do not impute ML support numbers from these site patterns.

## Public genomic records fetched successfully
P: PV752465.1 D. pseudomezereum var. pseudomezereum (171019 bp)
C: ON244034.1 D. pseudomezereum var. koreana / working D. kamtschatica var. koreana (171152 bp)
K: PV752464.1 D. kamtschatica (171100 bp)
O1: MN201546.1 D. laureola (171613 bp)
O2: KY991380.1 D. kiusiana (171491 bp)

The C plastome (Yoo et al. 2023) was obtained from Gangneung and the C ITS (Han et al. 2025) was obtained from Pyeongchang. They are not paired at individual level.

## MAFFT alignment and strictly ACGT 2:2 site patterns
| outgroup | aligned columns | complete 4-taxon ACGT | variable complete sites | CK 2:2 sites | PC 2:2 sites | PK 2:2 sites |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| D. laureola | 173356 | 169311 | 1585 | 27 | 4 | 0 |
| D. kiusiana | 173582 | 168974 | 1605 | 27 | 6 | 0 |

Ungapped pairwise identity (O1 alignment; denominator is pairwise comparable sites):
C-K 0.999708 (171068 columns), P-C 0.999327 (170845), P-K 0.999327 (170839).
These are **not directly equivalent to the identities cited in the original paper**, which may differ in site filters or alignment.

## Inverted repeat duplication check
Yoo et al. 2023 reported C: LSC 84963, IRb 41725, SSC 2739, IRa 41725. A direct reverse-complement identity check on the downloaded ON244034 FASTA confirmed the IRb and IRa intervals are perfectly complementary.
Assuming the published reference orientation (and validated IR coordinates), the 27 CK two-by-two informative positions partition:
- 24 in LSC
- 1 in IRb
- 1 in SSC
- 1 in IRa
This leaves 26 CK informative columns after exclusion of one duplicate IR copy, versus 4 or 6 PC columns depending on outgroup. The 3 adjacent PC columns at C positions 74892–74894 may reflect a local multi-base event or alignment context; avoid treating every nucleotide as an independent mutational event.

## Biological inference ceiling
The chloroplast quartet strongly favors C-K at the **genealogical marker level**, contrasting with 6–8 informative sites all supporting P-C in the nrITS quartet run 43677310692.
This is reproducible **gene-tree discordance**, not demonstrated chloroplast capture. Both chloroplast capture/introgression and ILS can create it; cpDNA is one linked inheritance unit, and nrITS is a multicopy repeat.
The next inferential unit must be multiple independent low-copy nuclear loci across multiple individuals, with coupled organelle/nuclear samples and explicit MSC / network comparison.
