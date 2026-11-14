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
python3 analysis/reproduce_shibata2020_seed_survival.py --check
python3 analysis/reproduce_shibata2020_paternity.py
python3 analysis/reproduce_shibata2020_siring_elasticities.py
python3 analysis/retained_female_function_scenarios.py --check
python3 analysis/synthetic_design_sensitivity.py --output docs/synthetic_design_sensitivity.md
```

各比率は**同じ論文・同じ種**の公表平均値間の記述的比較のみ。個体データ・標準誤差の共分散がないため統合メタ解析や新しい有意差検定を行っていない。

## 2027年野外実験計画
- [全体プロトコル](protocols/2027_field_protocol.md)
- [最新版v2：他家花粉量・総量・自家花粉生存性を分離する授粉実験](protocols/2027_pollen_identifiability_v2.md)（従来の[等量混合案](protocols/2027_pollen_interference_factorial.md)は識別上の交絡があるため参考資料へ変更）
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


## 2026-11-14 第2次追試：分布データの観測過程と分類改訂
- **2025年の新しい系統分類**：Han et al. (2025) は旧 *D. pseudomezereum* var. *koreana* を *D. kamtschatica* var. *koreana* とする分類を提案し、Kew POWOも受容。葉緑体とITSの系統不一致がある。オニシバリ・ナニワズとの比較で旧学名の集計を混同しない。[文献・送粉記録の批判的監査](docs/2026-11-14_taxonomy_pollination_critical_update.md)。
- **匂い・送粉者の前提を更新**：オニシバリの弱い香りの公的記述と、ハナアブ多数の訪花という非査読の野外観察例がある。「全く無臭・送粉者不在」と前提化しない。花粉輸送/結実寄与は未測定。
- **GBIF観測過程の感度を実計算**：[博物館標本・HUMAN_OBSERVATION・2000年以降・位置不確実性明記の層別と400回レコード数均等化](docs/2026-11-14_gbif_record_effort_robustness.md)。オニシバリ×ナニワズの0.25°占有格子Jaccardは全記録0.0446、標本のみ0.0438、HUMAN_OBSERVATIONのみ0.0250。各100レコードに揃えると中央値0.0241（5–95%反復幅0.0108–0.0482）。
- **位置品質上の制約**：報告済み位置不確実性が存在するのは、QC後オニシバリ347記録のうち21件だけ。現時点の格子差から真の接触頻度、繁殖様式・送粉者・性分化を推論しない。
- スクリプト：[記録層と件数を揃える再解析](analysis/gbif_daphne_effort_sensitivity.py)。Python標準ライブラリのみ、固定シード。実データはGBIF API、出力は**集約値のみ**。

## 2026-11-14：ナニワズ父性解析・長期適応度の一次資料を追加
- [Shibata（2020）博士論文の父性・長期繁殖コスト・近交弱勢の原著監査](docs/2026-11-14_shibata2020_paternity_lifetime_fitness.md)。**別種ナニワズの実測**でありオニシバリの新規観測ではない。
- 同論文の遺伝的父性では814子中482子（59.2%）が父性割当。雄的両性株の花数と父性成功が正に関連。平均花粉散布距離6.9–11.4m。**花の幅の因果効果ではない**。
- 他殖40種子→1年後38発芽→3年間31生存。自殖19種子→9発芽→6生存。累積近交弱勢 **δ≈0.593**。自殖は完全に不可能ではないが、子の適応度が低い。
- 7–11年追跡で果実数と翌年花数の負の関連を検出。ただし別の2022年の資源操作実験では明確な果実形成コストが見えなかった。**年数・因果性・集団差のある別設計として両論文を併記し、コストがゼロとも確定とも断定しない**。
- [再計算用種子データ](data/derived/shibata2020_seed_pipeline.csv)／[生存・近交弱勢検証](analysis/reproduce_shibata2020_seed_survival.py)、[父性割当元集計](data/derived/shibata2020_paternity_assignment.csv)／[父性件数検証](analysis/reproduce_shibata2020_paternity.py)、[雄成功度モデルの元係数](data/derived/shibata2020_sire_success_coefficients.csv)／[係数の効果倍率検証](analysis/reproduce_shibata2020_siring_elasticities.py)。
- **自家花粉干渉実験の旧比較を修正**：総花粉粒数を揃えるだけだと、活性のある他家花粉の量が違うので、単純なX+X2とX+Sの差は機構を識別できない。2027年試験は最新版v2の多処理と花粉管・胚・父性解析へ。従来の[仮想感度シミュレーション](docs/2026-11-14_design_sensitivity.md)のORも自家花粉の機構特異的な検出率ではない。

## 2026-11-14：残存雌機能が維持され得る条件を数値化
- [残存雌機能の仮想閾値モデル](docs/2026-11-14_retained_pistil_thresholds.md) と [再計算コード](analysis/retained_female_function_scenarios.py)。近縁ナニワズの**種子形成後の**生存率（他殖31/40、自殖6/19）から相対生存率r=0.40747、δ=0.59253、自己由来の単純な遺伝子コピーの帰結2r=0.81494を計算した。
- **あくまでも思考実験**：配偶者不足下で自殖種子S=0.5個/機会、通常の雄機能損失c=0.05、固定維持コストk=0という**架空の条件**なら、配偶者不足機会pが**10.9%超**で残存雌機能の仮想的な直接の遺伝的利益が正。自殖種子が一つも形成されなければこの利益は存在しない。真の性システムの進化的安定性は計算していない。
- [2027年の最小限の実証計画](protocols/2027_minimum_identifiability_workplan.md)。まず雄的両性花の他殖雌機能、袋掛け自律自殖と受精・発芽、翌年の投資コストを分けて確認する。**完全な雄花の対照が存在しなければ、「雌しべそのものの形成コスト」は直接測れない**。
- **空の記録フォーマット**：[花ごとの授粉操作](data/templates/2027_flower_reproduction_events.csv) と [種子ごとの発芽・生存追跡](data/templates/2027_offspring_survival_followup.csv)。実測値は一行も入れていない。
- CIは`python3 analysis/retained_female_function_scenarios.py --check`、自動生成Markdownとコミット済み出力の一致を検証する。2020年の種子・父性原著の検証と合わせて再現性を確保。

## 2026-11-14：送粉者不足と配偶者不足を分離（重要なモデル修正）
- [雌株F・雄的両性株H・完全雄型Mの比較を分ける研究解釈](docs/2026-11-14_two_channel_evolutionary_interpretation.md)。**F vs Hの雌株の維持**と、**H vs 仮想Mの雌しべ残存**は同じ問題ではない。
- 近縁ナニワズの2018年原著では、自殖由来の種子は形成されても**送粉者のサービスが必要**で、既存個体群で自殖による十分な繁殖保証はなかった。**他殖相手不在を送粉者不在と同一視できない**。オニシバリで自律自殖があるかは未測定。
- [新モデルの再現コード](analysis/two_channel_pollen_mate_sensitivity.py)：配偶者不在p、有効な花粉輸送v、訪花経由自殖S_visit、自律自殖S_auto、通常の雌機能F、雄機能損失c、維持コストkを区別。ナニワズの種子形成後の相対生存率r=0.40747だけを参考入力し、残りは全て仮定。
- 仮想条件（F=0.01、c=0.08、k=0.01、S_visit=0.5、S_auto=0）では、v=0.8なら残存雌機能が正となる配偶者不足pは**17.3%超**、v=0.2なら**25.1%超**、v=0なら利益なし。これは**オニシバリの環境依存性の実測でも進化的維持の推定でもない**。
- 同じナニワズ59種子の標本に基づく独立種子の仮定でもrのばらつきが大きい。母系間の非独立性は未補正。CIは[過程モデルの自動検査](.github/workflows/literature_check.yml)を実施し、全表を`daphne-pollinator-mate-two-channel-theory`の期間限定アーティファクトへ出力。
- **最優先の実証は、袋掛けだけの自律自殖と、人工自家授粉の自家和合性を区別すること**。次に他殖授粉後の雌機能、受精前干渉・自殖胚の死亡、雌機能保持に伴う雄成功度のコストを独立に測る。
- [旧10.9%閾値](docs/2026-11-14_retained_pistil_thresholds.md)は送粉経路を省略した思考実験として残し、適応証拠・検出力・性型平衡に使わない。

## 研究結果の運用ルール
1. *D. jezoensis*（ナニワズ）と *D. pseudomezereum*（オニシバリ）の実験値を混ぜない。
2. 袋掛けゼロ結実だけで自家不和合性と結論しない。
3. 同じ植物の**花の大きさ**と**株全体のサイズ**を区別する。
4. 除雄による改善・自家花粉の付着量・他家花粉管の到達・胚発達・発芽を分けて判定する。
5. コードが計算を終えても、オニシバリの新規生態学的実証が得られたと記載しない。
