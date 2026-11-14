# onishibari

**オニシバリ (*Daphne pseudomezereum*) の性表現・送粉・繁殖生態研究**

> Research draft. 現時点で**オニシバリ固有の人工授粉・送粉者・遺伝・花サイズに関する実測データは未収集**。以下の定量比較はすべて近縁の別種について公表された数値に基づく。

## 核心の問い
雌株より花が大きいとされる雄的な両性株で、雌しべが形態的に維持されながら、なぜ雌機能が低下するのか？ 特に**自家花粉の干渉 (pollen clogging)**、花粉不足、自家不和合性、近交弱勢、胚珠・胚の機能低下を区別する。

## 文献と解析
- [文献レビュー](docs/2026-11-14_evidence_review.md)
- [1989–2023年の自家花粉干渉仮説の証拠監査](docs/pollen_interference_critical_review.md)
- [4種の公表平均値から計算した13組の比較表](docs/comparative_results.md)
- [種名・文献・測定条件つきの公開数値CSV](data/derived/published_aggregates.csv)
- [再計算する Python スクリプト](analysis/compare_published.py)
- [自家花粉対他家花粉の混合実験・仮想感度解析コード](analysis/synthetic_design_sensitivity.py)

実行（Python標準ライブラリのみ）：
```bash
python3 analysis/compare_published.py --check
python3 analysis/compare_published.py --output docs/comparative_results.md
python3 analysis/synthetic_design_sensitivity.py --check --simulations 100
python3 analysis/synthetic_design_sensitivity.py --output docs/synthetic_design_sensitivity.md
```

各比率は**同じ論文・同じ種**の公表平均値間の記述的比較のみ。個体データ・標準誤差の共分散がないため統合メタ解析や新しい有意差検定を行っていない。

## 2027年野外実験計画
- [全体プロトコル](protocols/2027_field_protocol.md)
- [最優先：除雄花への他家＋他家対他家＋自家（等量）＋花粉付着順序の要因実験](protocols/2027_pollen_interference_factorial.md)
- [空の野外観測CSVテンプレート](data/templates/onishibari_observation_template.csv)

## 更新された直接根拠
- オニシバリ花粉のSEM径21–26µmは既存研究を2023年の花粉比較論文が引用するもの：<https://www.jstage.jst.go.jp/article/jjpal/68/2/68_68_2_KJ004/_article/-char/ja/>。本文2023年の新規測定はナニワズなので種を混同しない。
- 神奈川県自然環境保全センターによるオニシバリの**弱い香り**の観察記載：<https://www.pref.kanagawa.jp/docs/f4y/01about/ikimono/winter_trees.html>。香気成分や昆虫誘引は未測定。
- ナニワズ2017年は茎体積に性差を検出せず（F23.3、H24.3 cm³; p=0.68）、花数（F24.4、H39.3）、果実数（F2.2、H0.6）に差。オニシバリの値ではない：<https://doi.org/10.1093/aobpla/plw089>。
- 2026年*Silene acaulis*の研究では局所的な雄機能個体密度が雌機能を左右した：<https://doi.org/10.1002/ajb2.70282>。これをオニシバリの局所性比を測る動機とする。

## 2026-11-14 公開データの実取得（今回）
- [GBIF / iNaturalist の記録数・写真機会の監査](docs/2026-11-14_public_occurrence_audit.md): オニシバリのGBIF日本国内649記録（座標あり362）、2016–2026年の座標付き写真記録26。iNaturalistでは同期間に写真61観察、研究用28観察。**GBIFとiNaturalistは重複するので足さない。**
- [写真28観察の一次目視スクリーニング](docs/visual_photo_screen.md) と [ラベルCSV](data/derived/inat_visual_screen_v01.csv): 最初の写真で開花を確認14、果実4、花/蕾不明瞭1、葉/蕾のみ7、ライセンスによる未評価2。花の機能的雌雄は**全28件未判定**。
- [Daphne属の競合する繁殖戦略と反証可能な仮説](docs/comparative_sex_system_synthesis.md): 低い両性株の雌機能（ナニワズ）と、小型花・高自殖（コショウノキ）の違いを比較。
- [GBIF API取得コード](analysis/gbif_occurrence_audit.py)、[iNaturalist取得コード](analysis/inat_onishibari_audit.py)、[写真コンタクトシート生成](analysis/build_inat_screening_sheet.py)、[一次判定検証・月別集計](analysis/summarize_visual_screen.py)。GitHub Actionsで実取得と写真スクリーニング資料の生成に成功。
- **これは観測機会と既存画像の監査であって、オニシバリの自家不和合性・送粉者・雌雄の適応度の新規実験結果ではない。**

## 2026-11-14：全104写真の監査へ拡張
- **実取得**：[iNaturalist全写真監査](docs/2026-11-14_all_photo_review.md)、[61観察の手動ラベル](data/derived/inat_visual_screen_v02.csv) を追加。日本のオニシバリ写真付き観察61件（研究用28、同定中33）、画像参照104枚、CCライセンス等で閲覧用100枚を確認。
- **写真で直接確認**：開いた花23観察、成熟果実8観察、未熟果実2観察（4月・6月）。花/蕾不明2、葉・枝等のみ22、未評価4。性型は**全61件未確定**。
- **注釈の欠測を発見**：実際に花が見える23観察のうち、iNaturalistの開花注釈が付いているのは6件のみ（17件には注釈なし）。果実確認10観察のうち注釈あり5件。注釈ベースの花期研究には欠測への対処が必要。
- **再現性**：[全写真取得](analysis/inat_all_photo_audit.py)、[ライセンス付き画像コンタクトシート](analysis/inat_all_photo_gallery.py)、[61件のID・性型unknownの検証](analysis/validate_all_photo_review.py)をGitHub Actionsで実行。画像本体はGitにコミットせず、出典・ライセンス付きの期限付きartifactへ保存。
- **性システム記載の差を明示**：[撮影・解剖の判定基準](protocols/photo_sex_morph_scoring.md)。植物資料には「雌雄異株」と「雌性両全性異株」の両表記があり、雌しべの外形と実際の種子親機能は区別する。
- 画像で見つかる橙色の葯様構造は形態候補にすぎない。**雄的な両性株の自家和合性や雌機能を今回の写真から証明したわけではない。**

## 2026-11-14：直接の原著と広域生態へ拡張
- **1955年の一次資料を発掘**：[オニシバリ自体を含む花の性差・花期・種同定の原著レビュー](docs/2026-11-14_hamaya_1955_direct_evidence.md)。雄花/雌花の萼筒8–10/4–6mm、雄しべ1.5–2/1–1.5mmは**節内の複数近縁分類群に共通する概括**であり、オニシバリ単独の平均・サンプル数ではない。雌花の子房と柱頭位置の違いも記載。
- **原著の分類に注意**：オニシバリとナニワズの萼裂片／萼筒比は**ともに0.5〜1.0と重なる**が典型値が違う。葉の基部側脈など複数の形態証拠を統合する。Kew POWOでは別の受容種。
- **秋から翌春の花芽形成**：浜谷(1955)には前年果実・旧葉・新花序が同時に存在する例あり。Lei et al.(2020)は冬緑・夏落葉の炭素収支を実験的に検証。両者を結ぶ性型別の季節資源配分は新しい**未検証仮説**。
- **広域のGBIF実解析**：[3種の分布格子重複の実測とスケール感度](docs/2026-11-14_gbif_spatial_overlap_sensitivity.md)。精度・日本範囲をQC後、オニシバリ347記録/65格子、ナニワズ240/99、コショウノキ293/99。オニシバリとナニワズの記録Jaccardは0.1°=0.0163、0.25°=0.0446、0.5°=0.1275。**真のニッチ重複・同所性・種分化の推定ではない**。
- **再現スクリプト**：[1955原著の証拠範囲と現代写真の月の照合](analysis/audit_hamaya_photo_phenology.py)、[GBIF座標を公開せず格子数だけで分布重複を評価するコード](analysis/gbif_daphne_spatial_overlap.py)。GitHub Actionsで確認する。

## 研究結果の運用ルール
1. *D. jezoensis*（ナニワズ）と *D. pseudomezereum*（オニシバリ）の実験値を混ぜない。
2. 袋掛けゼロ結実だけで自家不和合性と結論しない。
3. 同じ植物の**花の大きさ**と**株全体のサイズ**を区別する。
4. 除雄による改善・自家花粉の付着量・他家花粉管の到達・胚発達・発芽を分けて判定する。
5. コードが計算を終えても、オニシバリの新規生態学的実証が得られたと記載しない。
