# Daphne cytonuclear project: public raw-read archive and nuclear-locus readiness audit
Date: 2026-11-15, reproducible evidence. This is a **data availability / method** advance, not detection of historical introgression.

## Source-verified 2023 Korean archive crosswalk
Original published whole-leaf-DNA sequencing: Yoo, Park & Oh (2023), DOI https://doi.org/10.1080/23802359.2023.2179356, sample voucher **Lee 7891**, Gangneung. Paper states **7.34 Gbp** original raw reads.
NCBI BioProject https://www.ncbi.nlm.nih.gov/bioproject/PRJNA835716 identifies 1 experiment (public **359 Mbases**). ENA portal API cross-checked **three** exact accession levels:
- BioProject `PRJNA835716` => 1 public run `SRR19117239`
- BioSample `SAMN28108557` => same 1 public run
- run `SRR19117239` => selfsame single public run
ENA run metadata: sequencing `WGS`, source `GENOMIC`, selection `RANDOM`, 1,191,530 **spots/pairs**, 359,842,060 bases. Thus publicly accessible sequence volume is **4.9025%** of the 7.34 Gbp described in the paper. Do not confuse ENA's `read_count` for spots with the 2,383,060 individual reads produced by `fasterq-dump`.

In our prior direct analysis of public SRA `SRR19117239`, the same 2,383,060 individual reads aligned **99.96% to the ON244034 same-voucher plastome** and **0 to the 2025 PX056848 ITS**. The publicly archived reads therefore cannot recover matched nuclear ITS; even though ENA calls the library GENOMIC/RANDOM, that is insufficient evidence of useful nuclear sequence content.

NCBI SRA name/voucher search: zero indexed SRA hits for 2025 specimen labels `Lee-Dpp230303` (Japanese P), `K005971094` (Russian K), `Byeon-Dpk250701` (Korean C), and for 2025 plastome accession IDs `PV752465`, `PV752464` or `ON244034` searched as strings. Searching `PRJNA835716` retrieved 1 SRA record. These are **negative indexed-search results**, not proof no alternative archives, privately retained sequences, or author-accessible raw reads exist.

Source GitHub run: https://github.com/zuizui0223/onishibari/actions/runs/43682778058
Artifact ID **14148703894**, ZIP SHA256 `c2d59da1ecc50d145ca42596e21e6d1f2928f301940bea57d5da646554e0d2ef`. Script `scripts/daphne_ena_archive_audit.py`.

## A practical nuclear-locus proxy resource created from 2026 family data
Bonnaudet et al. 2026 (`10.1093/aob/mcag206`), Zenodo `10.5281/zenodo.17517526`, source `3_AlignedLoci.tar.gz` (1,486,795 bytes, sha256 `bf35cb4f123db761a1c56d202b5adaa99c1180c6a55f21b1fb580ea1f4ffea05`) was parsed directly:
- **349 actual aligned nuclear locus FASTAs**, all parsed.
- **162** include all four published *Daphne* proxy samples `D. sureil` 27665_L001, `D. tangutica` 27673_L001, `D. mucronata` 27691_L001, `D. laureola` 07178_L001.
- **150** pass a deliberately simple screening rule: ≥200 4-way ACGT-complete aligned columns, ≥250 ACGT bases in each species, variable-site ratio 0.001–0.15.
- Top **50** selected by common 4-way complete column count and variation: **200 proxy sequences total** (four species × 50 loci). Minimum 4-way complete length **737 bp**, median **973.5 bp**. This can become a pilot alignment-reference panel for read recruitment from future original **P/C/K** unfiltered libraries.
- Output: `all_four_daphne_alignment_coverage.csv`, `candidate_proxy_loci_top50.csv`, `candidate_proxy_fourDaphne_sequences.fa`, and `proxy_panel_report.json`.
- Source run: https://github.com/zuizui0223/onishibari/actions/runs/43682931727
- Artifact ID **14149301848**, ZIP SHA256 `95c6986bf903215258c0fbb44ee6a8041753b2ec804703fd455dd9f0dd01fbab`.
- Script `scripts/daphne_353_proxy_loci.py`.
- **Interpretation ceiling:** proxy species **DO NOT include** P (=Japanese D. pseudomezereum), C (=Korean var. koreana), K (=D. kamtschatica). No gene tree for P/C/K recovered. The panel is not a verified single-copy ortholog set in the focal lineages. Recoverability and duplications must be checked against FASTQ and pilot assemblies.

## Why a strict 100 independent nuclear loci target is not sufficient in every scenario
Assuming neutral, panmictic three-species MSC, with internal branch t (in coalescent units):
- each of the two discordant gene-tree topologies has probability `exp(-t)/3`;
- total discordant gene-tree frequency is `2exp(-t)/3`.
Conditional on **n observed discordant independent gene trees**, two-sided exact binomial asymmetry test H0 proportion 50:50, illustrative alternative proportion 2:1 (not a biologically estimated effect).
- n=40 observed discordants => **53.0%** detection power (alpha=.05)
- n=60 => **66.4%**
- n=80 => **81.9%**
- n=100 => **90.3%**
To have *80 discordant gene trees on average*, simple t=0.25 needs ~155 total independent loci, t=0.5 needs ~198, t=1.0 needs ~327, t=1.5 needs ~538.
These calculations assume **true, independent gene trees without error** and fixed 2:1 asymmetry. They are a **method-design sensitivity illustration**, not a complete proof, an achieved sequencing target, nor a guarantee of test power at the expected locus count. An asymmetry can arise from population structure, gene-tree errors, unsampled donors, copy artifacts etc., not solely introgression. No nuclear asymmetry can rule out purely organelle capture.
- Source run: https://github.com/zuizui0223/onishibari/actions/runs/43683047543
- Artifact ID **14149111550**; SHA256 `8ad030c4b61916c8d8412a4a2084cc6f65e278abbc465f5d4b0f0b9ecfef0af2`
- Script `scripts/daphne_discordance_power.py`.

## Revised next step: source acquisition then adaptive inference
1. Ask Yoo et al. for the **unfiltered 7.34 Gbp** original read set for voucher Lee 7891, its sample sheet, and whether archived public SRR19117239 was restricted to organellar reads. Ask Han et al. for original P/K reads underlying the two new 2025 plastomes; verify matched specimen IDs and whether Korean C 2025 ITS specimen has raw genomic reads. Emails are drafts only, not sent.
2. After receiving access, use the verified 50-locus, 4-proxy-species panel to assess ≥20–50 recoverable low-copy nuclear targets, cross-mapping specificity, contamination, coverage, paralogy and genome-size implications. Do not treat low-coverage nonrecovery as proof of absence.
3. Adapt full capture design to **estimated** coalescent branch length and gene tree concordance precision. A nominal 100 loci is only an initial benchmark; may need ~200–500 independent loci for some gene-tree imbalance questions. Different models/network designs require their own power tests.
4. Pair organelle and nuclear loci from the **same vouchers**, replicate P/C/K/Hokkaido/Korean populations, separate *D. jezoensis* misindexed taxon records, and phenotype functional sex system independently of classification labels.
5. Until that information is available, the strong conclusions remain **ITS P–C, plastome C–K**, causal history unidentified.
