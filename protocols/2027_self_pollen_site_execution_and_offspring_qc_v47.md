# v0.47｜2027年オニシバリの同花・別花自家／他家受粉：花粉源から生存子までの実行監査

2026-11-17 JST。**設計・空の入力表・合成テストのみ**。オニシバリ Daphne pseudomezereum の実際の人工授粉結果、種子、生存子、交配親の遺伝子型はまだ得られていない。2027年の現地操作には許可・安全管理・同定済みの花と株が必要。旧5腕v1のデータをv0.47に入力してはいけない。v0.45という別の7役割実験案専用である。

## なぜ必要か

近縁のナニワズ D. jezoensis では、2008年のFFPRI年報No.49 p.40に同花花粉による人工自家では結実せず、別花・同株花粉による授粉では結実したとある。しかし処理別の母株数・種子の生存率は未公開。Kikuzawa 1989も自花粉による柱頭閉塞と除雄効果を実験済み。2001年の D. laureola の別花自家 × 除雄の原Table3はHわずか5独立株で、自己花粉処理の除雄効果の原著p=.058。これらは**焦点オニシバリの完全SIや雄的花の雌性稔性を直接示すものではない**。

区別すべき独立の要因：受容花の雌性能力、花粉の同株/別株由来、同花/隣花の花粉源、供与花の開葯・花粉稔性、袋の破損・花への傷害、花粉管、初期胚、熟した種子、発芽・生存、遺伝子型上の実際の父母。供与花を採粉に使うことは本人の別の花への自然花粉輸送の機会も奪い得るので、保全上の操作に数える。

## v0.45の割付けと、v0.47の実行監査は別工程

H（雄的両性形態）の7役割：open、bag_no_hand、bag_sham、bag_self_same_flower、bag_self_other_flower、bag_hand_cross、pollen_donor_same_genet。6つは母性結果を追う受容花で、最後の専用供与花は**母性果実・種子の結果を採用しない**。Fは旧5腕。H 1genetで確認蕾14個以上、Fで10個以上という最小条件は、累積50%・10蕾操作上限からの帰結。花齢や花粉稔性を確認できない場合に無理に花粉を使わない。独立単位はgenet、花と父株は巣状構造の内部単位。

3つの独立private CSVを使う：

1. v0.45 事前割付け（既存 generatorで作成）— 花の4ID、事前割付け arm、開放/袋掛け/供与専用、予定の同花/別花供与ID。実際の結果とは切り離す。
2. v0.47 実施ログ（空ヘッダー data/templates/2027_v45_pollen_source_execution_header.csv）— 花の4ID、実際のarm、実施/逸脱/未実施、操作日時、袋の破損状況、花粉生存QC、実際の花粉源genet/branch/flowerの3ID、供与者の遺伝的同一性確認、花粉粒数。
3. v0.47 後続果実・発芽・親子確認（空ヘッダー data/templates/2027_v45_parentage_germination_followup_header.csv）— 受容花の4ID、成熟果実、成熟種子数、実測発芽種子、想定した父母関係が遺伝子型上確認された発芽種子の数、親子同定QC、観察日。

これらの原本は公開GitHubではなく **data/private** 内で保存し、解析結果にも個体IDや精密位置を公開しない。v0.47は保全上の割付け数・armの整合と、実施された花粉源が事前指定と一致したかを照合する。別花自家では、指定された専用供与花の実施ログがあること、その花の利用が受粉以前であることも検証する。未施用・花粉源違い・未確認の同株遺伝型・袋破損・花粉活性不明は**予定どおりの処理として採用しない**。

## 成果の厳密な階層

- 柱頭花粉付着だけ、花粉管だけ、早期胚だけ、果実だけ、成熟種子だけ、発芽だけではそれぞれ別の機能段階。
- strict として数えるのは **遺伝子型上の想定母性と父性経路が確認された発芽済み種子**だけ。後続生存・繁殖まで達した生涯適応度とは異なる。
- 遺伝子型照合の「passed」を記録できない、実施ログが未確認、子の発芽が未評価の場合、strict endpointは**欠測／未確認**であって必ず0ではない。
- 手動他家受粉がHで成功0でも、Fの同時期同じ供与者条件での陽性対照が成立し、十分な独立H genets、過少成功の検出力を確認するまで「H完全不稔」「雌雄異株化完了」と結論しない。
- 機能的他家雌性子が生じても、花粉親としての成功と、雌しべの炭素維持費・発生制約・生涯繁殖を測らなければ雌しべの適応的維持は同定できない。

## 防御的な合成テスト結果

合成24蕾・2genetの試験で、v0.45の割付け H7役割、F5役割の12役割が正常に追跡できた。以下はすべて偽データのエラー検出テストであり、花粉実験を実施した意味ではない。

- 別花自家の花粉源花IDが指定供与花と違う → strictな per-protocol 比較から除外
- 後続観察が欠測 → 種子0に代入せず未確認として集計
- 親子同定で確認した発芽子が実測発芽数を超える → ハードエラー
- 花粉供与専用花に成熟果実の受容花結果を結合 → ハードエラー

**成功したCI** https://github.com/zuizui0223/onishibari/actions/runs/44014541314

## 再現と利用制限

Code: analysis/audit_2027_self_pollen_site_v47.py
Empty execution header: data/templates/2027_v45_pollen_source_execution_header.csv
Empty fruit/genotype header: data/templates/2027_v45_parentage_germination_followup_header.csv
Synthetic-only invocation: python3 analysis/audit_2027_self_pollen_site_v47.py --self-test

将来の非公開実施後だけ使用可能な例：
    python3 analysis/audit_2027_self_pollen_site_v47.py --schedule data/private/2027_self_pollen_source_assignments.csv --execution data/private/2027_v45_actual_execution.csv --followup data/private/2027_v45_actual_followup.csv --summary-private data/private/2027_v45_audit_summary.json

実際の入力や成果を含むファイルはGitHubにコミットしない。プログラムの成功ステータスは**入力形式と演算ルールの検証成功**であり、授粉・結実・遺伝的雌性機能の生物学的証明ではない。HとFの形態ラベル、実際の母性機能、実際の自家和合性、遺伝的父性と生涯適応度を区別したまま記録する。
