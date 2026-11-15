# Han et al. (2025) nrITS quartet reanalysis: observed results, not ILS/hybridization proof
Date: 2026-11-15 (JST)
Original: https://doi.org/10.11110/kjpt.2025.55.3.157
Source run: https://github.com/zuizui0223/onishibari/actions/runs/43677310692
Source artifact: daphne-nrITS-public-sequence-audit (artifact ID 14146782560; ZIP sha256 2ba6bfa659a4a04cb92b4715fbb5191e016b1eda4a6a3f7222c6fae0c8d4af09)
Executable: scripts/daphne_its_audit.py and .github/workflows/daphne-its-audit.yml

## Data, directly retrieved by NCBI eFetch in CI
- P = PX056852.1, D. pseudomezereum var. pseudomezereum, nrITS 586 bp, Japan
- C = PX056848.1, formerly D. pseudomezereum var. koreana (currently D. kamtschatica var. koreana), nrITS 579 bp, Korea
- K = PX056858.1, D. kamtschatica, nrITS 575 bp, Kamchatka
- O1 = MT623684.1, D. kiusiana, nrITS 780 bp, Kumamoto
- O2 = MT623688.1, D. kiusiana, nrITS 780 bp, Jeju
- O3 = MT623695.1, D. odora, nrITS 770 bp, Guangdong
- O4 = MT623699.1, D. papyracea, nrITS 780 bp, Guangdong

## Analysis
MAFFT --auto; align P/C/K + one outgroup at a time; 2:2 biallelic A/C/G/T parsimony-informative columns; IQ-TREE 2.0.7 -m GTR+F+I+G4 -B 1000 -alrt 1000.
Only a **single nrITS repeat region** was tested, so the four choices of outgroup do not count as four independent nuclear loci.

| outgroup | aligned columns | 4-taxon ACGT complete sites | variable complete sites | 2:2 PC | 2:2 CK | 2:2 PK | SH-aLRT PC/KO | UFBoot PC/KO |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| O1 D. kiusiana Japan | 797 | 571 | 62 | 6 | 0 | 0 | 97.2 | 99 |
| O2 D. kiusiana Korea | 797 | 571 | 63 | 6 | 0 | 0 | 98.0 | 98 |
| O3 D. odora | 787 | 551 | 53 | 6 | 0 | 0 | 97.2 | 99 |
| O4 D. papyracea | 797 | 567 | 56 | 8 | 0 | 0 | 99.4 | 100 |

All four yield unrooted split PC | KO; **with O as outgroup, P and C are sister in this four-taxon subset**. This supports robustness of the ITS topology against the selected outgroups, *not* nuclear-species-tree branching.

MAFFT alignment with O1 shows the six ACGT sites supporting PC at alignment 1-based columns **127,245,271,369,581,630**, distributed across the alignment. The original paper reported the raw pairwise identity as highest in C-K, which is not a contradiction: pairwise identity, site-pattern polarity, and full-taxon likelihood topology are different statistics, and alignment/overlap criteria differ.

## Interpretation
- **Robust within nrITS:** in this subset 6–8 informative sites all support PC, with high ML branch support.
- **Weak as a genome-wide inference:** nrITS is multicopy and subject to concerted evolution. One supported ITS genealogy cannot distinguish organelle capture, nuclear introgression, or ILS.
- **Missing paired plastid+nuclear observations for C:** nuclear C is from Pyeongchang (2025), plastid C ON244034 from Gangneung (2023).
- **Original fuller tree is not reproduced:** original 2025 analyses used other Daphne species and three more distant outgroups; our 4-tip analysis is a deliberately minimal, inspectable sensitivity experiment.
- **Inference next:** plastome site-pattern analysis with duplicated IR caveats, followed by multiple low-copy nuclear loci per individual and population. Explicit tests under multi-species coalescent with vs without introgression; compare tree concordance and network models. Genomic capture direction is not identifiable from the current two compartments alone.
