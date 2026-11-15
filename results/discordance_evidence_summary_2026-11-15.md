# Focal Daphne discordance audit — evidence and inferential limits (2026-11-15)

**Question:** Is the observed plastome vs nuclear ITS discordance in P = Japanese D. pseudomezereum, C = Korean D. kamtschatica var. koreana (originally D. pseudomezereum var. koreana), and K = D. kamtschatica caused by organelle capture/introgression, incomplete lineage sorting (ILS), or locus-specific nrITS history?

**Overall conclusion as of this analysis:** Discordance is robust **between the two linked genomic compartments**; the causal process is **not identifiable with a single multi-copy nrITS region and a single linked plastome, particularly since C was not sampled as paired individual across compartments**. This evidence does not permit ranking introgression above ILS.

## Public data and fully successful independent reanalysis
- Han et al. 2025: https://doi.org/10.11110/kjpt.2025.55.3.157
- Yoo et al. 2023: https://pmc.ncbi.nlm.nih.gov/articles/PMC9970196/
- nrITS workflow: https://github.com/zuizui0223/onishibari/actions/runs/43677310692 — artifact ID 14146782560 sha256 `2ba6bfa659a4a04cb92b4715fbb5191e016b1eda4a6a3f7222c6fae0c8d4af09`
- Plastome workflow with ML: https://github.com/zuizui0223/onishibari/actions/runs/43677644538 — artifact ID 14146976993 sha256 `71bb1a4d39fef834bbcb1fc7d68d35a18e10d8b1a0c5071fb6188a1223366968`
- GenBank flatfile metadata workflow: https://github.com/zuizui0223/onishibari/actions/runs/43677808968 — artifact ID 14146709725 sha256 `3444bad0a0f29d2d0abf69d43fc12b5c097aec78dfcb5573368c06cf41bf27c6`

### nrITS (P = PX056852.1; C = PX056848.1; K = PX056858.1)
Four-tip ML, one outgroup at a time, MAFFT and IQ-TREE 2.0.7, GTR+F+I+G4 with SH-aLRT 1000/UFBoot 1000:
| outgroup | alignment 2:2 sites PC | CK | PK | SH-aLRT of PC split | UFBoot |
|---|---:|---:|---:|---:|---:|
| D. kiusiana Kumamoto (MT623684) | 6 | 0 | 0 | 97.2 | 99 |
| D. kiusiana Jeju (MT623688) | 6 | 0 | 0 | 98.0 | 98 |
| D. odora (MT623695) | 6 | 0 | 0 | 97.2 | 99 |
| D. papyracea (MT623699) | 8 | 0 | 0 | 99.4 | 100 |
Only **6–8 informative nucleotide positions** determine the P-C split in the focal 4-taxon subanalyses. Altering outgroups does not change this within-marker result, but none of these are independent low-copy nuclear markers. Original paper used a broader taxon set; our minimal quartet is a sensitivity check, not a formal replication of that tree.

### Whole plastome (P = PV752465.1; C = ON244034.1; K = PV752464.1)
MAFFT and IQ-TREE 2.0.7, K3Pu+F+I+G4, with SH-aLRT 1000/UFBoot 1000:
| outgroup | CK 2:2 sites | PC | PK | SH-aLRT CK | UFBoot CK |
|---|---:|---:|---:|---:|---:|
| D. laureola (MN201546) | 27 | 4 | 0 | 100 | 100 |
| D. kiusiana (KY991380) | 27 | 6 | 0 | 99.6 | 100 |

C plastome (ON244034) has duplicated, reverse-complement IRb and IRa each 41,725 bp: when one IR copy is excluded, 26 CK-supporting 2:2 alignment positions persist. However, cpDNA is still **ONE linked gene genealogy** regardless of the number of aligned nucleotides or high bootstrap. Some nearby substitutions may be part of one mutational event.

### GenBank accession crosslink audit
Pulled original GenBank nucleotide flat-files for PV752465, PV752464, ON244034 via NCBI eFetch. **None of those three flat-files contains a direct DBLINK/BioProject/BioSample/SRA field.** This is only absence of an identifier in those GenBank **flat-files**, NOT proof that WGS raw reads or SRA submissions do not exist. For example the Korean 2023 paper independently reports SRA accession SRR19117239. To find P/K WGS fastq, conduct BioProject/NCBI SRA accession search using collectors/institution/title/vouchers or contact authors; match specimen IDs, not taxon name alone.

## How to distinguish mechanisms
1. Re-extract chloroplast and >=50 independent **low-copy nuclear loci** from the **same individual** in each lineage, across P/C/K multiple geographic populations, with 1–2 other relevant Daphne outgroups and contact-zone samples.
2. Test nuclear gene tree frequency spectrum and concordance under the multispecies coalescent; account for estimation error. Under simple 3-taxon panmictic MSC without migration, the **two discordant alternatives** have equal expected probabilities `e^-t/3`; repeated nuclear genealogies are necessary to estimate this.
3. Introgression signatures: compare D/ABBA-BABA where rooted unlinked SNPs and an appropriate outgroup exist; fit migration-network models (e.g., species networks, MSC+introgression) with different donor/direction candidates; check localized shared haplotypes and coalescent times. Nuclear D=0 does NOT rule out historical organelle-only capture.
4. Test secondary contact vs retained ancestral structure, ancestral population subdivision, nrITS repeat concerted evolution, cp uniparental inheritance, and unsampled ghost donor lineages.
5. Integrate morphological sexual dimorphism only after individual-level reproductive phenotyping and ensure taxonomic synonym harmonization (D. rebunensis not independent from D. jezoensis; Korean var. koreana not automatically all Japanese montane plants).

## What is **not** demonstrated
- No likelihood/Bayes-factor comparison of introgression vs ILS has been run.
- No high-resolution species tree has been estimated from independent nuclear loci.
- No timing/direction of gene flow has been established.
- No causation between flower sexual dimorphism, pollinator scarcity, and speciation has been demonstrated.

## 2026-11-15 update: same-specimen C nuclear ITS salvage failed, with control
Yoo et al. (2023) public raw reads SRR19117239 (same 2023 voucher as ON244034 plastome) were aligned to Han et al. (2025) PX056848 ITS: 0/2,383,060 individual reads mapped. Plastome positive control ON244034 produced 99.96% read-level mapping. See `results/korean_2023_sra_its_zero_validated_2026-11-15.md` and https://github.com/zuizui0223/onishibari/actions/runs/43678829837 . Therefore, the 2023 Korean specimen still has no observed matched nrITS; **do not imply the discordant nuclear/plastid trees come from the same individual**. A reference-free WGS rescue cannot solve this from public SRR19117239 because it is effectively chloroplast reads only. New 2026 genus-level HybSeq work confirms a different European Daphne allopolyploid but not focal Korean introgression: https://doi.org/10.1093/aob/mcag358 . The reported BioProject PRJNA1493253 did not appear in the NCBI SRA phrase-index query on this date.
