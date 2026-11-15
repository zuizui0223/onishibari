# オニシバリ・Daphne属 横断比較 v0.1（2026-11-15）

## 目的
オニシバリ Daphne pseudomezereum を中心に、花サイズの性的二型・性表現・自家和合性・送粉制限・生活史・遺伝的分岐の共進化（または非同期性）を検証する。
**今回完了した範囲:** 文献・分類データの初期監査、解析単位の整理、公開アクセッションの抽出。配列の実ダウンロード・再アラインメント・新たな系統推定は**未実施**。

## 重要な事実と注意
1. D. pseudomezereum の機能的性システムと自家不和合性は、今回確認した一次交配実験では確定できていない。古典的な「雌雄異株」記載と機能的雌雄異株を区別する。
2. D. jezoensis は形態的に雌性両全性異株（雌株＋両性株）、後者は雄機能優勢。14集団の研究で両性株は通常**自家和合性**だが結実能力が低い。雌性両全性異株が日本の地理的範囲で比較的保守的。性表現は完全固定ではなく6年間の追跡で約15%が変化。
3. D. kiusiana は両性花を持ち、自殖可能。同一分類群内で韓国南西部（西）と韓国南東部・九州（東）の遺伝系統と花形質が異なる。東側で花が小型化し自殖関連の形質・遺伝的指標が顕著。推定された自殖シンドロームの成立は少なくとも約10万年前。
4. D. miyabeana は自家和合性で自動自殖も可能。開放地と閉鎖林で自然結実率は約31.22%と9.01%。訪花不足への対処が単純な性分化に限定されない比較例。
5. D. laureola はイベリア半島の集団間で雌性両全性異株／両性型が分化し、13集団のITS・葉緑体解析は複数の移行を示唆。**同一種内の系統地理学**であり、種分化イベントと混同しない。
6. D. kamtschatica var. koreana（旧 D. pseudomezereum var. koreana）についてHan et al. 2025は、葉緑体で D. kamtschatica に、核ITSで D. pseudomezereum に近縁と推定。交雑は一仮説であり、ILS・ITS固有の問題も考慮する。
7. Han et al. 2025は本州（主に関東山地）・四国の「var. koreana」相当集団について、朝鮮産型とは異なり未記載分類群の可能性を指摘。**日本の山地型を朝鮮産 var. koreana と自動統合しない**。
8. D. rebunensis はKew POWOで D. jezoensis の異名、D. koreana は D. kamtschatica var. koreana の異名。旧名のまま独立種数を水増ししない。
9. Daphneの広義の属境界に関して核ITS/葉緑体による系統解析の違いがあり、主題のMezereum周辺を核多遺伝子系統で解決することを優先。
10. 古典的な花サイズ差の記述は**花冠でなく萼筒（calyx tube）**を指す。単一研究の代表値（雄 8–10 mm、雌 4–6 mm）を全種・全個体群の効果量として扱わない。

## 研究上の3つの反証可能な対立仮説
- H1 性システムの保守性: 遺伝・生態的分化が大きくても性表現の実測指標は変わらない。D. jezoensis の集団間比較が参照例。
- H2 繁殖保証への転換: 送粉機会の不足下で自家和合性、花の縮小、雌雄器官近接が共変化する。D. kiusiana の東西系統比較が参照例。
- H3 過去の交雑と形質モザイク: 葉緑体と複数核遺伝子の系統不一致、性機能・花サイズの祖先形質の不一致が地理的な接触帯に集中する。Han et al. 2025 の ITS vs plastome は出発点だが、これだけでは交雑を確証できない。

## 最低限必要なデータ設計
1. Taxon backbone: Kew POWO accepted_name, taxon rank, synonym, geographic population unit, voucher, year。
2. Sex axis: morphological sex (F / H / M), functional sex (母性適応度 vs 父性適応度), individual sex switching, within-population sex frequency。
3. Flower axis: calyx-tube length, lobe length/area, flower number, total floral display, herkogamy; F/H と M/H を混同しない。
4. Mating axis: controlled outcross vs self vs bagged autonomous self; fruit set, viable seed set, offspring survival, inbreeding depression。
5. Pollinator axis: observed visits per flower-hour, day/night, pollen deposition, effective pollen transfer; 0 visits is pollination absence proofではない。
6. Genomics: plastome and nrITS first, then >=数十独立の核座位を各集団複数個体で取得し、MSC and network / introgression testsを比較。
7. Divergence: geographic isolation, niche, phenology, interpopulation cross compatibility; differentiation FST≠speciation。種分化率の推定はサンプルサイズ不足のため当面行わない。

## 分類学上の解析単位（当面）
- 東アジア焦点: D. pseudomezereum (Japan nominate), D. jezoensis, D. kamtschatica var. kamtschatica, D. kamtschatica var. koreana (Korea), Honshu/Shikoku mountain plants (working OTU), D. miyabeana, D. kiusiana (western/eastern lineages).
- 広域対照: D. laureola (Iberian populations), D. mezereum, D. gnidium, D. rodriguezii.
- D. rebunensis/D. koreana は過去の文献名として synonym に保持するが、受容種として数えない。

## 監査対象資料（一次研究優先）
- Hamaya 1955, Notes on Thymelaeaceae (historical calyx tube and sexual morphology): https://www.jstage.jst.go.jp/article/jjapbot/30/2/30_30_2_3828/_pdf/-char/en
- Sinclair et al. 2016 D. jezoensis male-biased hermaphrodites: https://doi.org/10.1111/plb.12463
- Shibata et al. 2021 D. jezoensis 14 populations: https://doi.org/10.1007/s00606-021-01780-1
- Shibata & Kudo 2022 sex-differential plasticity negative test: https://doi.org/10.1111/plb.13458
- Han et al. 2023 (online Dec 2022), D. kiusiana selfing syndrome: https://doi.org/10.1093/aob/mcac142
- Sakata & Nakahama 2018 D. miyabeana: https://doi.org/10.1111/1442-1984.12212
- Alonso & Herrera 2011 D. laureola: https://doi.org/10.1111/j.1558-5646.2011.01246.x
- Han et al. 2025 D. kamtschatica var. koreana reclassification: https://doi.org/10.11110/kjpt.2025.55.3.157
- Kew POWO D. jezoensis: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:831255-1
- Kew POWO D. kamtschatica var koreana: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77371659-1

## 次の実解析 (着手前)
1. GenBank 6 accession の核 ITS 3本と plastome 3本を実取得し、配列長・分類名・voucher 監査。
2. 追加Daphne 近縁分類群・複数個体の独立核遺伝子候補を収集、欠損率と混入を確認。
3. 系統不一致のsupport検証（partition-specific bootstrapや核遺伝子 gene concordance）、交雑の対立説明を整理。
4. 文献証拠を実測/二次記載/未検証に分離、花形・性表現を population×year 単位で抽出。
