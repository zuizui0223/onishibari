# SRA availability and evidence-grade audit, 2026-11-15

## Verified source
Workflow: https://github.com/zuizui0223/onishibari/actions/runs/43678442679
Artifact 14147790376, ZIP sha256 `363feeaf85c860f94af4d62154c26caed81070bd77fce01f43e9be3be99f8aed`.
Actual source: NCBI SRA ESearch/ESummary (`db=sra`, quoted exact species names, `[All Fields]`), captured as XML with `scripts/daphne_sra_index_audit.py`.

| Query | Indexed SRA experiments | Important qualification |
|---|---:|---|
| "Daphne pseudomezereum" | 3 | **Not 3 nominate Japanese D. pseudomezereum**; they are D. jezoensis study (2) and var. koreana (1) |
| "Daphne pseudomezereum var. pseudomezereum" | 0 | No directly indexed nominate SRA under this phrase |
| "Daphne kamtschatica" | 0 | No directly indexed K SRA under this phrase |
| "Daphne pseudomezereum var. koreana" | 1 | `SRR19117239` from Korean 2023 plastome study |
| "Daphne koreana" | 0 | Synonym search did not retrieve C directly |
| "Daphne jezoensis" | 2 | `DRR240745`, `DRR015021`; indexed organism name is misleadingly D. pseudomezereum |
| "Daphne rebunensis" | 0 | Synonym, no hit |
| "Daphne kiusiana" | 237 | Example titles 'CP-seq of Daphne kiusiana', library type AMPLICON; **not a 237-sample low-copy-nuclear phylogenomic set** |
| "Daphne reichsteinii" | 0 | New paper published Nov 10, 2026 reports sequencing deposit |
| "PRJNA1493253" | 0 | Reported new 2026 Daphne hybrid nuclear study project, not indexed by query yet |

## Named reusable runs
- Korean C plastome-source SRA `SRR19117239`, 359,842,060 bp, **effectively chloroplast-only in public reads (99.96% aligned to ON244034, zero nrITS reads)**, Illumina NovaSeq 6000, title 'The first complete chloroplast genome of Daphne pseudomezereum var. koreana'; BioProject PRJNA835716, BioSample SAMN28108557, voucher Lee 7891 (Gangneung); source: Yoo et al. 2023 https://doi.org/10.1080/23802359.2023.2179356
- `DRR240745` 1,761,565,001 bp, Illumina MiSeq, 'Genome Partial Sequencing of Daphne jezoensis' DRA010662 (specimen SAMD00003142; sequencing run).
- `DRR015021` 83,141,212 bp, Ion Torrent, same D. jezoensis reference specimen/project DRP001343.
- Neither `DRR240745` nor `DRR015021` is safe to label **D. pseudomezereum var. pseudomezereum** solely from their NCBI Organism name; the SRA *study and experiment title* label D. jezoensis.
- New paper Krtilová et al. (2026) `10.1093/aob/mcag358` reports `PRJNA1493253` with Angiosperms353 Hyb-Seq and sdRAD, including D. reichsteinii, D. alpina, D. petraea, other European Daphne; focal eastern P/C/K absent. NCBI index query did not confirm publicly accessible runs on 2026-11-15.

## New primary taxonomy/sexual-system evidence
Nitta & Uchida (2020), J. Jpn. Bot. 95:341–347, https://www.jstage.jst.go.jp/article/jjapbot/95/6/95_95_6_11054/_pdf/-char/en :
- D. kamtschatica in Hokkaido as well as Kamchatka, described with *hermaphrodite flowers*, long clonal rhizomes, 2n=18 in both sampled regions, flowers early spring, leaves present year-round.
- Morphologically coded D. pseudomezereum as *dioecious*; D. jezoensis as *gynodioecious*. Functional sex must be experimentally remeasured in target populations.
- Do not repeat old historical statement that all Mezereum species are dioecious: later taxonomic and floral analyses conflict.
- D. kamtschatica var koreana in Korea may differ from Japanese montane samples historically referred to same name; keep population provenance.

## Feasible next analyses
1. **Completed and failed nuclear salvage**: SRR19117239 nuclear nrITS recovery yielded 0/2,383,060 mapped individual reads (all six P-C diagnostic sites depth zero). Positive-control remapping of the same reads to the ON244034 chloroplast yielded 99.96% read-level alignment. The public SRA is effectively chloroplast-only/enriched and **cannot resolve the 2023-versus-2025 Korean specimen mismatch**. This does not imply the original *unfiltered* leaf DNA raw library lacks nuclear reads. See `results/korean_2023_sra_its_zero_validated_2026-11-15.md`, successful CI run 43678829837.
2. For true mechanism testing, collect >=3 geographic populations per operational taxon and ~5–10 individuals per population (allow fewer for rare sites), including P Honshu, C Korea, K Hokkaido and Kamchatka and D. jezoensis as adjacent compare; sequence paired plastomes and 100+ independently targeted nuclear loci from each same individual using Angiosperms353 Hyb-Seq (reference genome not required), plus ploidy.
3. Test P/C/K independent nuclear gene trees with MSC and introgression network models; cp capture can leave little nuclear trace and D-statistic may be nonsignificant even if capture occurred. Include evidence for and against ancient ancestral structure and ghost donor lineages.
4. Functional sex should be coded by self/outcross controlled crosses, seed viability, pollen export, observed flower dimorphism, and flowering timing: 'hermaphrodite' morphology and male/female performance are different axes.

## Non-inferences
- No independently sampled multi-locus nuclear dataset for all P/C/K was identified in this exact-name SRA search.
- No conclusion that all public sequencing for any taxon is absent.
- No direct probability ratio of ILS vs introgression is available.
