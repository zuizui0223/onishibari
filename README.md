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

## 研究結果の運用ルール
1. *D. jezoensis*（ナニワズ）と *D. pseudomezereum*（オニシバリ）の実験値を混ぜない。
2. 袋掛けゼロ結実だけで自家不和合性と結論しない。
3. 同じ植物の**花の大きさ**と**株全体のサイズ**を区別する。
4. 除雄による改善・自家花粉の付着量・他家花粉管の到達・胚発達・発芽を分けて判定する。
5. コードが計算を終えても、オニシバリの新規生態学的実証が得られたと記載しない。
