# New nuclear reference: Bonnaudet et al. (2026) Thymelaeaceae dataset, audited 2026-11-15

## Provenance
- Paper: Bonnaudet, Maurin, Forest, Baker, Gaudry, Rogers, Roy, Schoeplein & Anest (2026), *Annals of Botany* 138:1956–1975, DOI: https://doi.org/10.1093/aob/mcag206
- Data: Olivier Maurin, *Go and grow...*, Zenodo DOI: https://doi.org/10.5281/zenodo.17517526
- Reproducible source/taxon audit successful GitHub Actions run: https://github.com/zuizui0223/onishibari/actions/runs/43681600027
  - Artifact `daphne-2026-family-hybseq-inventory`, ID 14148079580, ZIP SHA256 `6599b89fb30faea893a3d82d693138dd17c24ec4745448c44e718478e0504bb1`
- Reproducible **gene-tree quartet audit** successful run: https://github.com/zuizui0223/onishibari/actions/runs/43681845833
  - Artifact `daphne-2026-nuclear-gene-tree-context`, ID 14148686605, ZIP SHA256 `82a6c6f6d9be651aed9781875b059104690ec0b4dd9115eacbcbd9bc65a957e1`
- Executable: `scripts/zenodo_phylogeny_inventory.py`, `scripts/daphne_2026_nuclear_context.py` and corresponding workflows.
- Source integrity: Zenodo 89-tip ASTRAL tree SHA256 `04ea7ca8ff434edfc6e73fda6abd00d97f2dbf73244474bc08ff2d9c4c358f7b`; dated NEXUS tree SHA256 `6cb9d42f1c12d38161c01f181e98f7ec866f809e8d5a314943af40da80da15df`; 349 IQ-TREE gene-tree archive SHA256 `afdc4869a532f3dbbdc26e8056516c83ce84a1fe2919e3916bcf929ef5fb593f`.

## What is available
Zenodo record `17517526` directly contains (sizes verified against Zenodo API):
- 1_Readme.txt: 286 bytes
- 2_Loci.tar.gz: 4,991,640 bytes
- 3_AlignedLoci.tar.gz: 1,486,795 bytes
- 4_IQTree.tar.gz: 2,527,850 bytes
- 5_iqtree-BS10_sp_V2.tre: 20,371 bytes, 89-tip ASTRAL species tree using specimen IDs
- 6_ThymelaeaceaeDated.tre: 60,841 bytes, 89-tip NEXUS chronogram containing taxon labels and a TRANSLATE table
The genus **cannot be identified from the 89-tip ASTRAL tree labels alone** (e.g. `27665_L001`); crosswalk with the chronogram's specimen-level taxlabels is required.

### Exactly four `Daphne` tips in the published 89-tip Angiosperms353 species tree
| Taxon | ASTRAL terminal ID | NEXUS label |
|---|---|---|
| Daphne sureil | `27665_L001` | `Thymelaeoideae_Daphne_sureil_27665` |
| Daphne tangutica var. tangutica | `27673_L001` | `Thymelaeoideae_Daphne_tangutica_var_tangutica_27673` |
| Daphne mucronata | `27691_L001` | `Thymelaeoideae_Daphne_mucronata_27691` |
| Daphne laureola | `07178_L001` | `Thymelaeoideae_Daphne_laureola_7178` |

**No P = D. pseudomezereum, C = Korean D. kamtschatica var. koreana, or K = D. kamtschatica is sampled.** Other relevant *Daphne* absent include D. jezoensis, D. kiusiana and D. miyabeana.
All four sampled Daphne form an exclusive clade in the ASTRAL species tree; induced unrooted quartet split: (D. sureil, D. tangutica) | (D. mucronata, D. laureola).

## Direct 349-locus gene-tree audit
Read **349 distinct IQ-TREE per-locus gene trees** from the deposited `4_IQTree.tar.gz` (1047 archive members including .log and .iqtree ancillary files). 162 gene trees included all four Daphne, 187 lacked >=1 Daphne; no parsing failures.

| Induced unrooted quartet split | Count of 162 | Fraction |
|---|---:|---:|
| sureil + tangutica vs mucronata + laureola | 151 | 93.21% |
| sureil + mucronata vs tangutica + laureola | 10 | 6.17% |
| sureil + laureola vs tangutica + mucronata | 1 | 0.62% |

The four-Daphne quartet is robust in this family-wide published dataset, with most locus trees agreeing on its unrooted split. **These are different, more deeply sampled *Daphne* lineages than the East Asian P/C/K trio.** It does **not** resolve whether P/C/K cytonuclear discordance reflects introgression or ILS. Also, gene tree estimation error, missingness (187/349 without all four), and model violations may distort minor-quartet frequencies. In the 11 discordants, frequencies 10 versus 1 are asymmetric (conditional exact binomial p = 0.0117 under an ideal equal-alternatives model), but **do not interpret this alone as evidence of hybridization** in these four species or in P/C/K.

The family-wide dataset is useful as (i) phylogenomic context for broad Daphne evolution, (ii) positive control for Angiosperms353 target capture across the family, (iii) reference alignment sources for exploratory ortholog identification if P/C/K total DNA reads become available, and (iv) demonstration of specimen crosswalk/coverage QC.

## Major new limitation and data recovery opportunity
The Yoo et al. 2023 Korean plastome paper states **7.34 Gbp of raw reads were originally generated** from total leaf DNA (https://doi.org/10.1080/23802359.2023.2179356). The linked public SRA `SRR19117239` contains **359,842,060 bases**, only 4.90% of that stated raw data volume. Our previous same-voucher QC mapped **99.96%** of these 2,383,060 available reads to ON244034 plastome, but 0 to ITS. This suggests the public deposited run is enriched/selected for plastid sequences, NOT that unfiltered original total DNA data were never produced. **Request the original unfiltered 7.34-Gbp FASTQs from the authors** before commissioning a new Korean nuclear sequencing experiment. Do not conflate reported raw read volume with currently public SRA volume; no original raw file has been obtained.
Likewise request original Japanese P and Kamchatkan K Illumina total-DNA FASTQs underlying Han et al. 2025 PV752465/PV752464, which were not identifiable as direct public SRA matches in our earlier audit.

## Life-form reminder
`D. pseudomezereum` is a **summer-deciduous woody shrub (落葉小低木)**, not a herbaceous perennial. This is biologically important for survival, clonality comparisons and lifetime reproductive selection. The Tsukuba Botanical Garden explicitly describes it as a deciduous subshrub and notes that 'dioecy' is not sharply demarcated: https://tbg.kahaku.go.jp/recommend/illustrated/result.php?mode=easy&name=%E3%82%AA%E3%83%8B%E3%82%B7%E3%83%90%E3%83%AA%EF%BC%88%E3%83%8A%E3%83%84%E3%83%9C%E3%82%A6%E3%82%BA%EF%BC%89&order=ruby&p=1&pn=4&season=3b .

## Inference gate
Until independently sampled low-copy nuclear loci exist **for P, C, and K from matched individuals**, historical chloroplast capture, ITS-only introgression/concerted evolution, and ILS remain observationally compatible. Do not combine these 162 nonfocal quartet loci with the 2025 focal ITS/Cp as if they were the same taxa.

## Independent row-level QC supplement
All 349 gene-tree assignments are now inspectable in `per_gene_coverage_and_quartet.tsv` (one row per gene tree, 349 rows and 4 columns; 162 covered, 187 missing >=1 Daphne). Latest successfully executed workflow run: https://github.com/zuizui0223/onishibari/actions/runs/43682044800 ; artifact ID **14148374481** (`daphne-2026-nuclear-gene-tree-context`), ZIP SHA256 `0024a44f3a6bc7459b067e25f522d77b8fdab0d1127d1fcddf4e8ee88aa8052c`. Counts independently rechecked by opening source ZIP: 151 + 10 + 1 = 162, total 349. This detailed table is an audit trail, not new data from the focal East Asian trio.
