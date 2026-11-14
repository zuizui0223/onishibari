# オニシバリの公開写真による形態的性型の評価基準 v0.2

Created 2026-11-14. **写真から自家不和合性・花粉稔性・雌性適応度を推定しない。** オニシバリ直接の実験的繁殖成績は未取得。

## 最初に解決する分類上の問題

- 「雌雄異株」：名古屋大学博物館および複数の植物図鑑の表記。観察者向けには結実する雌株と雄的な株の区別として実用的。
- 「雌性両全性異株 (gynodioecy)」：観察図鑑 jplants はオニシバリに**雌株＋両性株（機能的にはほぼ雄株）**を記載。雄的な両性株に雌しべが形態的に存在していても、生き残る種子を残すかは別問題。
- 形態的雌しべ残存と、花粉管・胚・成熟種子という機能を区別しなければ、雌雄異株・雌性両全性異株の分類は一意にならない。
- 近縁種 *D. jezoensis* の人工授粉・性型追跡結果を、*D. pseudomezereum* についての直接実証として転用しない。

原典：
1. https://www.jplants.jp/onsibari.html （雌性両全性異株、雄しべ8個・4個×2段、雌花の退化葯、両性株の低結実という記載。非査読の観察図鑑）
2. https://www.num.nagoya-u.ac.jp/media/plants/mar.html （「雌雄異株」と記載。弱い香りあり）
3. https://hanatachi.sakura.ne.jp/pl1c/jincyouge/sub5/onishibari_b.html （雌雄別株、雄花・雌花の比較写真）
4. https://pmc.ncbi.nlm.nih.gov/articles/PMC5497021/ （近縁ナニワズで4/8葯が花筒外へ突出するという学術論文）

## 形態的な最低限の観察指標

| field | admissible states | rule |
|---|---|---|
| `flower_stage` | bud/open/senescent/no_flower/unclear | 萼の裂片が開いているか、ピンぼけで不明か |
| `anthers_visible` | mature_large / reduced / not_visible / uncertain | 花筒内の2段の雄しべ、葯の色とサイズが明瞭な場合のみ |
| `pollen_grains_visible` | visible / not_visible / uncertain | 粒子を確認できなければ「ない」とせずuncertain |
| `stigma_visible` | visible / not_visible / uncertain | 正面・断面写真でのみ評価、見えない=存在しないではない |
| `ovary_or_ovule_visible` | visible / not_visible / uncertain | 通常の公開外観写真ではほぼ評価できない |
| `fruit_visible` | ripe / unripe / absent / uncertain | 赤い果実の有無。成熟果実はその時点の雌機能の実現を示すが、生涯一定の性型は証明しない |
| `sex_morph_photo_candidate` | male_biased / female_like / indeterminate | 大型かつ成熟した葯が明瞭なときmale_biased候補、退化葯・花の葯付近が明瞭なときfemale_like候補、その他indeterminate |
| `functional_sex` | unknown / field_anchored_herm / field_anchored_female | 写真単独なら**unknown**。花粉稔性・人工授粉/果実・反復年の個体IDに基づく場合のみfield_anchored* |
| `confidence` | high / medium / low | 撮影方向、フォーカス、複数写真間の一貫性を基準に |

画像が複数ある場合は、写真レベルの評価を保持し、同一観察ID内で統合する。観察が同一株かどうかは独立に未確定。開花写真のサイズ差については**尺度がない画像からmm差を直接推定しない**。

## 手動判定の工程
1. 公開iNaturalistの同定を一次確認する。オニシバリ、ナニワズ、コショウノキが混入していないか、花・葉・萼の形態を確認する（識別不能な場合はuncertain）。
2. ライセンスのある写真のみ閲覧用コンタクトシートに掲載。各画像の出典/ライセンス/著作者情報は添付TSV。無ライセンスや転載不可の写真を再配布しない。
3. 画像を**開花/果実/葉/不明**と、**葯可視/柱頭可視/不可判定**の二段階で手動判定する。性型が不明のときは空欄にせず`indeterminate`。
4. 各観察について確認できる全写真を見て、二名による独立判定を推奨。食い違うものは保留して原寸で解決する。
5. 性型頻度を公表する場合、写真の選別確率が性型で異なる可能性を評価する。写真の雌雄候補比を自然集団の性比と呼ばない。
6. 完全な自家不和合性を主張するには手動自家授粉・他家授粉・花粉管/胚発達・発芽データを分離する必要がある。袋掛けゼロ結実は証明ではない。

## 追加すべき実証データの優先順位
- 個体IDの両性/雌型、萼幅・萼筒長、葯のサイズ・花粉の有無・雌しべの縦断面。
- 同じ株における自殖、他殖、混合花粉、除雄対照での他家花粉管の到達率、果実・成熟種子・発芽率。
- 雌株1–3–10m以内の同期的に開花する両性株密度と花粉付着/自然結実（花粉追加で回復するか）。
- 昼夜の訪花者、実際の柱頭接触、花粉運搬、GC-MSの香気。

## 研究上の重要な限界
**一枚の花の写真から「自家不和合性」「雄機能優位」「雌しべが退化しない進化的理由」を確定することはできない。** 本プロトコルは資料画像を研究に活用できるかを監査するためのもの。

## 2026-11-14：1955年原著による実測形態・同定の再校正

浜谷稔夫 (1955), [Some taxonomical notes on Thymelaeaceae from Japan and adjacent regions](https://www.jstage.jst.go.jp/article/jjapbot/30/2/30_30_2_3828/_pdf/-char/en) は**オニシバリ自身を含む近縁5分類群**の花と葉を観察し、雄花と雌花の識別に花サイズだけでなく**柱頭・下段雄しべ・子房の位置関係**が有用と記した。

- **節を通じた概略**：雄花の萼筒8–10 mm、雌花4–6 mm、雄しべの長さは雄1.5–2 mm、雌1–1.5 mm。これは**オニシバリだけの平均・SD・検定結果ではない**。種ごとに標本数は記されていない。
- **雌花**：萼筒に対して子房が大きく、柱頭は下段雄しべに届く高さになる例が多い。柱頭の位置を外から見えないときは判定保留。これだけでは雄花の雌しべの種子形成能力は分からない。
- **種同定**：萼裂片長/萼筒長は**オニシバリもナニワズも0.5〜1.0の範囲で重なる**。典型的にはオニシバリは0.5寄り、ナニワズは1.0寄り。絶対閾値で同定しない。
- **葉脈**：オニシバリでは基部の側脈が不規則で早く枝分かれし、ナニワズではより規則的。例外・中間型があるので花と葉の複合証拠を用いる。
- **開花月**：1955年検索表はオニシバリ12–3月、ナニワズ4–5月とするが、地域・標高・年次による重なりを許容し、開花月一つで分類を変更しない。
- **最新の分類基盤**：Kew POWOはオニシバリとナニワズを別の受容種として記載。2020年の未発表分子系統「一種の南北型」を確定データとして採用しない。

### 画像判定の必須フラグ

追加すべきメタデータ：
`taxon_determination_basis`（花/葉/標本/推定のみ）、
`lobe_tube_ratio_measurable`（yes/no）、
`lobe_tube_ratio`（単位なし、画像比からの目安で絶対mmと混同しない）、
`basal_lateral_venation`（irregular/regular/unclear）、
`ovary_stigma_lower_stamen_relation`（reaches/below/unclear）、
`independent_confirmed_taxon`（yes/no）、
`sex_morph_functional`（必ずunknownから始める）。

iNaturalist 61観察のうち33件は `needs_id` であるため、**形態未確認を保留したまま、性型や花期の種間差を算出しないこと**。 [詳細な原著監査](../docs/2026-11-14_hamaya_1955_direct_evidence.md) を参照。
