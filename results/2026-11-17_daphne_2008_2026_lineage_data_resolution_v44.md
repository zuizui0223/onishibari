# v0.44｜2008–2011年ジンチョウゲ属の「常緑→夏緑→冬緑」系統仮説を原資料別に再監査

Updated 2026-11-17 JST. Source-first audit of **published statements**, not original 2010 specimen sequences and not a new molecular species tree. Target *Daphne pseudomezereum* (Onishibari) and relatives.

## Main result: historical species branching is based on what kind of sequence evidence?

| Source | Original wording/data scope | Which independently evolving genealogies? | What can be defended |
|---|---|---|---|
| KAKEN 20380093 2008 annual (https://kaken.nii.ac.jp/report/KAKENHI-PROJECT-20380093/203800932008jisseki/) | Evergreen *D. kiusiana/miyabeana* diverged far from wintergreen *D. pseudomezereum/jezoensis* and summergreen *D. koreana*; inferred from genetic analysis but no accession or nucleotide-site summary | Unspecified; **not independently auditable** | Original **reported topology/phenology hypothesis** |
| **KAKEN 2009 annual** (https://kaken.nii.ac.jp/report/KAKENHI-PROJECT-20380093/203800932009jisseki/) | Explicitly lists chloroplast rbcL, matK, trnL–trnF, psaB–rpl4, rpl14–rps16, trnH–psbA as analysed for Daphne phylogeny; interprets evergreen early, then summergreen, then wintergreen Onishibari | **Six regions of one linked plastid lineage**; no low-copy nuclear locus described in this annual phylogeny methods | One reported plastid-genealogical hypothesis; not a multilocus nuclear species tree |
| Japanese Forest Society 2010 conference B08 (https://doi.org/10.11519/jfsc.121.0.332.0) | Seven barcoding targets, **six plastid regions plus multicopy nrITS2**; reports comparative barcode screening | At most plastid genealogy + one nuclear-rDNA repeat history **if the same vouchers were actually sequenced and both markers analysed together**; no linked original accession/voucher table retrieved | DNA barcode study is documented; **not seven independent nuclear loci** |
| KAKEN 2010 final report issued 2011 (https://kaken.nii.ac.jp/file/KAKENHI-PROJECT-20380093/20380093seika.pdf) | Evergreen basal, summergreen Korean-Naniwazu then wintergreen Onishibari/Jezoensis; four focal Daphne as gynodioecious with more limited self-compatibility in wintergreen P/J | No novel individual matched nuclear-locus details, cross treatment N, surviving seeds or accession links provided | **Historically claimed phenotype/phylogeny correlations**, insufficient to independently rank character transitions and causal mating mechanisms |
| Lei et al 2019 online/2020 issue (https://doi.org/10.1111/ppl.12972) | P and J treated as southern/northern races of one species **based on an unpublished phylogeny credited to T Kawahara**; Japanese montane *D. koreana* from Nara/Tokushima | An underlying analysis was cited, but original trees, taxon-voucher crosswalk and loci unavailable. Paper authorship says Kawahara inferred a chloroplast/nuclear phylogeny; **cannot recover independent low-copy nuclear gene trees from that statement** | Proposed conspecific P/J not a proven delimitation result |
| Han et al 2025 (https://doi.org/10.11110/kjpt.2025.55.3.157) and previously run original source audit (PR#2) | Focal P Japanese ITS PX056852, C Korean ITS PX056848, K ITS PX056858 yields P-C; P plastome PV752465, C plastome ON244034, K plastome PV752464 yields C-K | plastome = **one linked inheritance history**; nrITS = **one multicopy rDNA history**; Korean C plastome vs nuclear ITS sourced from different individuals and localities | Focal **marker-level cytonuclear discordance** robust; not introgression/ILS ranking |
| Bonnaudet et al 2026 family Angiosperms353 (https://doi.org/10.1093/aob/mcag206; https://doi.org/10.5281/zenodo.17517526) | 349 nuclear locus gene trees; 162 include all FOUR sampled Daphne; 151 support main proxy quartet | 162 proxy gene genealogies but **none of focal P/C/K** | Demonstrates broad Daphne target-capture feasibility, not focal three-species history |

## The new 2010 raw-accession attempt in v0.42–43

Original team author/period search of NCBI nuccore: **zero** entries indexed to exact Kawahara (2008–13) or Kawahara/Yamashita/Yoshimura (2008–15). A broad ITS2 query returned 60 records, plastid spacer query two. **All 62 indexed candidate EFetch GenBank flatfiles have now been checked**: two GU736359/360 map to Zhang/Volis/Sun Tibet-wide *Stellera* phylogeny, the other **60 indexed ITS2 records** (multiple accession-series candidates) to Alonso/Herrera *Daphne laureola* mating-system study. **0/62** matched the Japanese 2010 project's first-reference authorship. These negative **indexed searches and nonmatching records within this specific query** do not prove the original FASTA does not exist, was never submitted, or cannot be retrieved under alternative author names / years / taxa.

Original-source machine results:
- v0.42 https://github.com/zuizui0223/onishibari/actions/runs/44011060026
- v0.43 first 14-record validator PASS https://github.com/zuizui0223/onishibari/actions/runs/44011439558; expanded full-62 source originally refitted https://github.com/zuizui0223/onishibari/actions/runs/44011840289; see current v43 workflow for revised full-62 validator
- Source-coded full report results/2026-11-17_daphne_2010_barcode_sequence_recovery_audit_v43.md.

## Mechanism implications

1. **Do not infer a unidirectional sex system series (hermaphrodite → gynodioecy → dioecy) merely from the 2011 plastid-only hypothesis.** In particular the focal P male-flower functional pistil status is unresolved quantitatively despite its reportedly low self-compatibility.
2. **Do not infer leaf habit drives sexual differentiation** from evergreen vs wintergreen categorical associations on a short, uncertain reticulate tree. Canopy, day/night pollination, mate availability, growth, resource state and ancestry are confounded.
3. **Do not equate all instances of "Daphne koreana":** Japanese 2020 montane vouchers (Nara/Tokushima; summergreen) are not proven identical to Korean *D. kamtschatica var. koreana*. Historical *D. kamtchatica var. jezoensis* refers to current *D. jezoensis*, not narrow *D. kamtschatica*.
4. **Do not multiply plastid marker count into independent species genealogy count:** rbcL/trnL–trnF etc contribute many nucleotides of one organelle, nrITS2 adds only an rDNA repeat history.
5. A defensible phylogenetic and eco-evolutionary test needs **same-genet plastome + independent low-copy nuclear loci** for focal P/K/C, geographically replicated, correct voucher/taxon concepts; coupled **hand-self/outcross functional seed outcomes** for P and independently verified morphological sexual-system states for each population.

## How to recover the historical primary data

- Request 2010 Kawahara et al's original matrix, 7-marker FASTA and chromatograms, cited 2010 phylogeny poster, old herbarium voucher codes and modern taxonomy crosswalk. Explicitly identify whether their phylogenetic tree used any nuclear evidence beyond ITS2 and whether sampled P, J, and Japanese mountain "koreana" are paired across loci.
- Request original 2008–10 mating experiments: Onishibari male-biased H and female F per-plant counts for hand-self, hand-outcross, bag only, any embryos, germination. Exact original evidence remains accessible only as an unquantified project report; **the literature is not wholly silent**, but its data provenance is unresolved.
- Request unfiltered Lee/Oh/Han P/C/K FASTQ source reads; original 2023 Korean study reported ~7.34Gbp, whereas public SRA SRR19117239 is ~360Mb and essentially plastid-only. Existing 2026 familywide proxy nuclear genes cannot fill focal gene-tree gaps.
- Contact drafts prepared but **not sent** in docs/author_unfiltered_read_requests_unsent.md; do not claim authors granted access.

## Conclusion

The 2010 Japanese Daphne project is an important historical **hypothesis-generation source** for ecological phenotype covariation, and its raw sequences are a valid recovery target. It is not yet an independent nuclear phylogenomic replication of the 2025 plastid-ITS discordance. Conversely, the 2025 discordance alone cannot rank chloroplast capture, nrITS introgression, incomplete lineage sorting, retained ancestral structure or hybrid history. For Onishibari, treat **morphological sex classification, physiological self-compatibility, viable maternal offspring and actual genetic paternal success as distinct variables**. No new focal cross or independent nuclear sequence obtained in this update.
