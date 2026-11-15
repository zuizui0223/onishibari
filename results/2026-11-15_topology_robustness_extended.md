# Daphne cytonuclear discordance: topology robustness audit (2026-11-15)

## Status
**Newly completed direct likelihood comparison.** All three quartet gene trees were optimized using the original 2025 focal GenBank accessions, two unlinked *compartments* but not two sets of independent nuclear loci. P=Japanese D. pseudomezereum; C=Korean D. kamtschatica var. koreana; K=Kamchatka D. kamtschatica.

Source run: https://github.com/zuizui0223/onishibari/actions/runs/43680208874
Artifact `daphne-three-topology-tests` ID 14147788090, SHA256 `924d27b69238aa544f5bfff303879fa99dd1e3d6ae727a4e981aeb7be29708b5`.
Script `scripts/daphne_topology_robustness.py` / workflow `.github/workflows/daphne-topology-tests.yml`.

MAFFT quartet alignment per marker: ITS (P PX056852, C PX056848, K PX056858, O MT623684 D. kiusiana); cp (P PV752465, C ON244034, K PV752464, O MN201546 D. laureola). Model: ITS GTR+F+I+G4, plastome K3Pu+F+I+G4. IQ-TREE 2.0.7 optimize 3 user trees `-z -n 0 -zb 10000 -au`.

| Unrooted split | ITS logL | ITS ΔlogL | ITS p-AU | plastome logL | plastome ΔlogL | plastome p-AU |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P–C (versus K–O) | -1414.777079 | 0.000 | 0.978 | -245280.0748 | 55.384 | 0.000122 |
| C–K (versus P–O) | -1423.254168 | 8.477 | 0.0228 | -245224.6912 | 0.000 | 1.000 |
| P–K (versus C–O) | -1423.254042 | 8.477 | 0.0222 | -245286.4151 | 61.724 | 2.14×10^-12 |

**Conclusions limited to these alignments**: ITS strongly supports PC and rejects CK and PK under the chosen nucleotide model (AU<0.05). Plastome strongly supports CK and rejects PC and PK. The tests do not compare demographic mechanisms (ILS versus introgression); two different genomic compartments can reflect different genealogies even without hybridization. The cp genome is one linked, typically uniparentally inherited unit; IR repeats duplicate some sites and per-site resampling does not grant independent genomic genealogies.

## Diagnostic ITS positions
Audit of the original 4-outgroup MAFFT alignments saved in artifact from run https://github.com/zuizui0223/onishibari/actions/runs/43677310692:

| C (PX056848) 1-based unaligned position | P and C nucleotide | K and O1/O2/O3/O4 nucleotide | within 5 aligned columns of any gap? |
| ---: | --- | --- | --- |
| 46 | C | T | no |
| 160 | G | A | no |
| 184 | C | T | no |
| 282 | T | C | no |
| 491 | G | A | yes |
| 540 | C | T | no |

These six positions occur across all 4 outgroup choices. With D. papyracea (O4), two **additional adjacent positions C134 and C135** appear informative for PC but are near alignment gaps and should not be overcounted as independent mutation events. Thus 5/6 core positions are far from gaps and 1 near a gap. They are 6 nucleotides within ONE linked multicopy nrITS region, not independent nuclear loci.

## Expanded taxon sampling
Completed expanded GenBank accession series MAFFT/IQ-TREE analysis:
https://github.com/zuizui0223/onishibari/actions/runs/43680405552

The `MT623670–MT623710` accession series contains 5 named `Daphne` ITS species with valid near-complete sequences: D. championii MT623676, D. genkwa MT623680, D. kiusiana MT623684, D. odora MT623695, D. papyracea MT623699. With focal P, C, K = **8 named OTUs**, not an exhaustive genus tree. Artifact ID 14147383995; zip SHA256 `115ee1c0b60d568ebb143a5fef99cffda84256ec5e4716f428cd1418c1b120d5`.

- MAFFT full ITS alignment length 826 columns.
- P–C formed an exclusive sister pair in the inferred unrooted tree; supporting internal edge SH-aLRT 97.1 / UFBoot 99.
- Filtering columns to >=80% taxa with ACGT calls retained 570 columns. P–C remained a sister pair, edge SH-aLRT 97.1 / UFBoot 99.
- The set is a limited taxon-sampling sensitivity check and includes divergent members (e.g. D. genkwa): independent nuclear species tree has NOT been resolved.
- Taxon names from GenBank must later be mapped to POWO with synonyms. In particular GenBank labels and accepted genus assignments can differ across older studies.

## Causal inference still missing
1. Exact same-individual plastid and nuclear sampling for C: C plastome ON244034 (Gangneung; 2023) vs C ITS PX056848 (Pyeongchang; 2025) are from different individuals.
2. Independent nuclear loci across multiple specimens from each P/C/K geography. The public 2023 Korean SRA SRR19117239 is effectively plastid-only: 0/2.38 million reads mapped to ITS, 99.96% to ON244034 plastome under positive control (run 43678829837).
3. MSC species tree, gene concordance and symmetrical-discordance tests, network / admixture alternative models, ploidy and chromosome markers, functional sex phenotype verification.
4. Nothing in these two independently analyzed **genomic compartments** is definitive evidence for ancient introgression, chloroplast capture, incomplete lineage sorting, or reproductive-system-mediated speciation.
