# onishibari

**オニシバリ (*Daphne pseudomezereum*) の性表現・送粉・繁殖生態研究**

> Research draft. 現時点で**オニシバリ固有の人工授粉・送粉者・遺伝・花サイズに関する実測データは未収集**。以下の定量比較はすべて近縁の別種について公表された数値に基づく。

## 核心の問い
雌株より花が大きいとされる雄的な両性株で、雌しべが形態的に維持されながら、なぜ雌機能が低下するのか？ 特に**自家花粉の干渉 (pollen clogging)**、花粉不足、自家不和合性、近交弱勢、胚珠・胚の機能低下を区別する。

## 文献と解析
- [文献レビュー](docs/2026-11-14_evidence_review.md)
- [1989–2023年の自家花粉干渉仮説の証拠監査](docs/pollen_interference_critical_review.md)
- [4種の公表平均値から計算した比較表](docs/comparative_results.md)
- [種名・文献・測定条件つきの公開数値CSV](data/derived/published_aggregates.csv)
- [再計算する Python スクリプト](analysis/compare_published.py)

実行（Python標準ライブラリのみ）：
```bash
python3 analysis/compare_published.py --check
python3 analysis/compare_published.py --output docs/comparative_results.md
```

各比率は**同じ論文・同じ種**の公表平均値間の記述的比較のみ。個体データ・標準誤差の共分散がないため統合メタ解析や新しい有意差検定を行っていない。

## 2027年野外実験計画
- [全体プロトコル](protocols/2027_field_protocol.md)
- [最優先：除雄×人工他殖＋花粉付着順序の要因実験](protocols/2027_pollen_interference_factorial.md)
- [空の野外観測CSVテンプレート](data/templates/onishibari_observation_template.csv)

## 研究結果の運用ルール
1. *D. jezoensis*（ナニワズ）と *D. pseudomezereum*（オニシバリ）の実験値を混ぜない。
2. 袋掛けゼロ結実だけで自家不和合性と結論しない。
3. 同じ植物の**花の大きさ**と**株全体のサイズ**を区別する。
4. 除雄による改善・自家花粉の付着量・他家花粉管の到達・胚発達・発芽を分けて判定する。
5. コードが計算を終えても、オニシバリの新規生態学的実証が得られたと記載しない。
