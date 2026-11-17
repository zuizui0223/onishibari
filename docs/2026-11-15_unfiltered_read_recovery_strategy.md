# Minimal-cost next action: recover original unfiltered nuclear reads from 2023/2025 Daphne libraries

Date: 2026-11-15. **This is an evidence-based acquisition protocol, not an accomplished data recovery. No email has been sent.**

## Why this is now highest priority
GenBank/NCBI public records provide conflicting **plastome** and **nuclear nrITS** genealogies: C–K plastome, P–C ITS. AU tests, alternate outgroups and IR de-duplication confirm both marker signals, but they cannot identify ILS vs gene flow.

A crucial new quantitative provenance discrepancy is confirmed:
- Yoo, Park & Oh (2023), DOI https://doi.org/10.1080/23802359.2023.2179356, explicitly describe **7.34 Gbp raw DNA sequencing** for a Korean var. koreana leaf sample (voucher **Lee 7891**; Gangneung; cp ON244034). Input was total DNA, sequenced Illumina NovaSeq.
- SRA `SRR19117239` linked to that study contains only **359,842,060 bases**, equal to **4.90%** of the 7.34-Gbp reported original raw base volume. It may be a *filtered/selected subset*; the raw amount and deposited amount differ.
- We directly downloaded and aligned the 2,383,060 SRA reads: 99.96% mapped to the same specimen's plastome ON244034, while **zero mapped to nuclear ITS PX056848**. This suggests SRA sequence is overwhelmingly plastid. **Do not conclude original unfiltered raw genomic reads lacked nuclear data.**
- The Korean nrITS used in Han et al. (2025), `PX056848`, is from **Byeon-Dpk250701, Pyeongchang**, not Yoo et al. 2023 Lee 7891. Thus P–C ITS and C–K cp analyses use different Korean **individuals**.

## Decision tree
1. Query the authors or associated archives for **unfiltered original 7.34-Gbp total-DNA FASTQs** for Korean **Lee 7891**. Ask for sequencing sample sheet / all runs, if present, and whether plastid-only reads were deposited in SRR19117239.
2. Query Han et al. (2025) for **unfiltered paired-end Illumina total-DNA FASTQs** used for:
  - Japan P = **Lee-Dpp230303** (Mie, March 2023; cp PV752465 and ITS PX056852).
  - Russian K = **K005971094** (Kamchatka, 2003 herbarium; cp PV752464 and ITS PX056858).
  - Korean C = **Byeon-Dpk250701** (Pyeongchang, July 2025; nuclear ITS PX056848); assess whether this individual was sampled for total-DNA sequencing (published ITS was obtained by PCR/Sanger).
3. Once obtained, perform strict read provenance, FastQC, deduplication and k-mer/genome-size estimate. Test read coverage against a panel of *Daphne* Angiosperms353 alignments from Bonnaudet et al. (2026, DOI https://doi.org/10.5281/zenodo.17517526); only 4 nonfocal Daphne taxa are in that public 89-tip nuclear reference, but 349 separate locus trees/alignments are provided.
4. Nuclear ortholog recovery is conditional: a genome skim at ~7.34Gb may or may not yield usable coverage of ≥100 loci, depending on genome size, library preparation, contamination, organellar filtering and phylogenetic divergence. Start with a **pilot** of 20–50 independently mapped low-copy loci and estimate missingness/paralogy before extending.
5. If adequate unfiltered P/K/C reads are not accessible, sequence NEW same-voucher P/C/K/Hokkaido populations with Angiosperms353/target capture rather than infer gene flow from plastid+nrITS alone.

## Technical preregistration
- Separate *gene-tree inference* from *species-tree inference*: one linked cp plastome and one multicopy ITS ribosomal array provide at most **two compartment-specific genealogies**, not a nuclear multi-locus sample.
- Only estimate MSC branching + independent gene concordance and uncertainty from many **unlinked low-copy nuclear loci**. Avoid concatenated SNP pseudo-replication.
- Infer hybridization only if multiple independent nuclear loci / haplotype patterns favor migration or explicit MSC+introgression fits over no-migration alternatives; check ghost lineages, ancestor structure and locus-estimation errors.
- After genealogy, measure sex expression and flower size as functional, individual-level traits; do not treat historic 'dioecious' labels as demonstrations of strict functional dioecy.
- `D. pseudomezereum` is a woody **summer-deciduous shrub (落葉小低木)** rather than an herbaceous perennial; survival and clonality must be measured as woodiness/life-history traits.

## Read request outline (not sent)
**Han et al. (2025) investigator**: ask whether raw unfiltered Illumina total-DNA FASTQs for specimens Lee-Dpp230303 (P) and K005971094 (K), and any analogous data for Korean Byeon-Dpk250701 (C), have existing public SRA IDs or can be shared for a comparative ortholog/cytonuclear analysis. Ask permission conditions/credit and sequencing library design.

**Yoo et al. (2023) investigator**: point out that paper states 7.34Gb initial total-DNA sequencing, whereas SRR19117239 has 359.8Mb of public reads, 99.96% of which align to plastome. Ask whether the *original unfiltered genomic* FASTQs from voucher Lee 7891 are archived and can be reused for nuclear low-copy locus extraction, and whether other data records exist.

No sending or representation that author permission has been obtained.
