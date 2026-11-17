# Failed recovery of paired nuclear ITS from Korean 2023 plastome source: positive-control validated

Status: completed 2026-11-15; **negative result for nuclear ITS**.
Source specimen: Korean historical var. koreana voucher Lee 7891 (Gangneung). Plastome ON244034; NCBI SRA SRR19117239.
Original paper: Yoo et al. (2023), https://doi.org/10.1080/23802359.2023.2179356
Comparison: Han et al. (2025), nuclear ITS PX056848 is a **different** plant Byeon-Dpk250701 from Pyeongchang.
Executable source: `scripts/daphne_same_individual_its.py` and `.github/workflows/daphne-same-individual-its.yml`.

## Reproducible result
- [Initial run (ITS)](https://github.com/zuizui0223/onishibari/actions/runs/43678673355)
- [Successful positive-control run](https://github.com/zuizui0223/onishibari/actions/runs/43678829837)
- Positive-control artifact ID `14146689880`, ZIP SHA256 `b1bdb96a7649589521c6bac77e026414e8fa4cb4584ded7dda701bc4edf16337` (contains pipeline logs, nrITS position-by-position counts and masked consensus)
- Input: `SRR19117239`, 1,191,530 read pairs = 2,383,060 individual reads, 359,842,060 total bases.
- Bowtie2 --very-sensitive-local vs 2025 Korean PX056848 nrITS 579 bp: **0/2,383,060 reads mapped** (0.00%). 0/579 positions with any read coverage, all six P-C diagnostic sites have depth zero. This is **not** evidence that 2023 and 2025 nuclear genotypes differ.
- Bowtie2 --very-sensitive-local positive control vs same-voucher plastome ON244034: **99.96%** overall alignment (only ~1,050 of 2,383,060 reads unmapped, from bowtie2 read-level log). Samtools flagstat on mapped-only BAM confirms 2,382,010 mapped reads.
- Note on paired alignment: only 0.34% plastome pairs mapped concordantly; many reads mapped in discordant or multimapping configurations, consistent with plastome repeats and circular-genome/read presentation. Therefore use **read-level** mapping rate, not properly-paired percentage, as positive control.
- Data thus appears effectively chloroplast-only or overwhelmingly chloroplast-enriched in its *public SRA form*. It cannot provide the same-individual nrITS haplotype needed to resolve the 2025 specimen mismatch.
- Do **not** conclude nuclear DNA was never present in the original total-leaf-DNA library: cp read extraction/filtering before SRA submission is possible.

## Resulting decision
Exclude public `SRR19117239` from the focal nuclear ITS and independent low-copy nuclear gene analyses.
The plastome `ON244034` remains useful as chloroplast-marker data. If investigators have archived **original unfiltered total DNA reads** from the Gangneung voucher Lee 7891, request these specifically rather than re-downloading the SRA.
For Han et al. (2025) P and K specimens, request Illumina raw total DNA reads associated with vouchers Lee-Dpp230303 and K005971094 and their provenance; published plastomes alone cannot replace 100+ unlinked nuclear loci.
For future sampling, use same-voucher paired plastome/nuclear target enrichment, multiple individuals from each taxon/region, including Hokkaido K and South Korean C, and assay morphology/function/ploidy separately.

## Inferential boundary
This experiment is a successful **quality-control demonstration of source limitations**, not a successful nuclear marker reconstruction or discrimination of chloroplast capture vs ILS.
