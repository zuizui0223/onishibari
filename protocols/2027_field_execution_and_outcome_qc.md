# オニシバリ 2027年：割付け・実施・成熟種子を分離するデータQC

最終更新 2026-11-14。**観測・操作は未実施**。このページは、将来入力する野外データに対して、同定、無作為化、処理実施と雌機能の判定を分離するための再現手順。

## なぜ必要か
*Daphne jezoensis* では人工的に他家授粉を行っても59–91%の雄的両性株が果実を作らないという研究がある（Shibata et al. 2018; https://doi.org/10.1007/s10265-017-0978-5）。しかしこれをオニシバリでも実証するには、**観測なし・操作失敗・実際の結実ゼロ・未熟果実・発芽未実施を明確に区別**しなければならない。近縁 *D. rodriguezii*、*D. gnidium* の比較では、柱頭花粉の発芽率や果実形成、株内形質変動も異なり、花粉数・花粉管・種子成功は同一指標でない（Rodríguez-Pérez & Traveset 2011, https://pubmed.ncbi.nlm.nih.gov/20820845/）。

## ステップ0：現地蕾の同定と割付け
- 採集/操作許可、花・葉の証拠標本、各株の性型ラベルを現地で確認する。`sex_morph` は `male_biased_morph`, `female_morph`, `undetermined`、後者は割付け対象外。機能的性型はこの段階で判定しない。
- `data/templates/2027_bud_inventory_header.csv`の列で1行1蕾の実数を記録し、`population_id,genet_id,branch_id,flower_id`が重複しないことを監査する。
- `analysis/randomize_2027_pilot.py`で各genetの飛散前蕾の**最大50%（かつ最大10蕾）**を、各型の5処理へランダム割付けする。5処理を満たせない個体は無理に割り付けず `SKIPPED_...` として保存する。閾値の変更は事前登録し、結果を見てから変更しない。
- どの処理も `ASSIGNED_NOT_APPLIED` で始まる。**割付け済み≠授粉実施済み**。未操作の蕾を野外の自然繁殖のために残す。

## ステップ1：操作ログ（必須）
[`data/templates/2027_pollination_execution_header.csv`](../data/templates/2027_pollination_execution_header.csv) を同じ4IDで作る。

| 項目 | 仕様・解釈 |
|---|---|
| `actual_treatment` | 実際の処理名。未実施時は空欄 |
| `execution_status` | `performed` / `missed` / `deviation` |
| `execution_date` | 実施日。実施なしは空欄 |
| `bag_integrity` | `intact` / `breached` / `unknown` / `not_applicable` |
| `pollen_viability_qc` | `passed` / `failed` / `not_checked` / `not_applicable`。人工授粉の操作では `passed` がないとper-protocol品質未確認 |
| `cross_donor_genet_id` | `bag_hand_cross` と `open_plus_cross` で必須。母株と同じIDなら他殖と認めない |
| `self_donor_genet_id` | `bag_hand_self`で必須。母株IDと同一でなければ自殖と認めない |

母株と別IDでも**同一クローン、実際の遺伝的近縁・交配不和合**はあり得る。これらはID照合だけでは排除できないため、**他家供与者を複数独立genetで確保し、花粉の活性と供与者×母株の交配相性をパイロットで確認**する。使った花粉ロットと顕微鏡での品質判定は非公開実施ログへ保持する。バッグ破れや操作変更は「失敗例を除外して終了」ではなく、ランダム割付けに属するデータとして記録し続ける。

## ステップ2：成熟・発芽追跡
[`data/templates/2027_mature_fruit_followup_header.csv`](../data/templates/2027_mature_fruit_followup_header.csv) に花ごとに追跡。別花での花粉管破壊採取は**ここへ成熟果実データとして入れない**。

| 状況 | 記録 | 解析 |
|---|---|---|
| 評価完了、果実なし | `followup_status=completed`, `fruit_mature=0`, `mature_seed_count=0` | 観測された0として集計 |
| 評価完了、果実/種子あり | `completed`, `fruit_mature=1`, `mature_seed_count`に観測数 | 雌機能の候補。受精起源は遺伝子型と汚染対照で確認 |
| 観察継続中 | `pending`, 果実・種子欄は空 | **未観測**。0として埋めない |
| 果実/タグ紛失 | `lost`, 結果欄は空 | 欠測。果実・種子成功率の分母と欠測機構を別途監査 |
| 花粉管や胚の破壊採取 | `destructive_sampling`、果実結果なし | 別段階のoutcomeとして解析 |
| 成熟種子あり、発芽試験まだ | `seed_germinated_count`は空 | 発芽数0と区別。種子ごとの発芽日・生存年数は独立の家系追跡表へ |

**研究設計上の核心**：人工他家授粉をした雌型Fが成功しなければ、雄的両性株Hで0が観察されても陽性対照が成立しない。逆にHで種子が得られれば、F側の成否にかかわらずHに**少なくとも実験条件下での雌機能が存在する可能性**を支持する。ただし同定、遺伝的父性、種子の生存確認が前提。

## ステップ3：入力と監査（実測データ取得後）
```bash
# 架空データによる内部テスト。野外結果を生成しない
python3 analysis/randomize_2027_pilot.py --self-test
python3 analysis/audit_2027_pilot_outcomes.py --self-test

# これ以降のファイルは全て非公開の data/private/ 内
python3 analysis/randomize_2027_pilot.py \
  --inventory data/private/2027_bud_inventory.csv \
  --output data/private/2027_assignments.csv --seed 20271114

python3 analysis/audit_2027_pilot_outcomes.py \
  --assignment data/private/2027_assignments.csv \
  --execution data/private/2027_execution.csv \
  --followup data/private/2027_mature_fruit.csv \
  --out-qc data/private/2027_flower_qc.csv \
  --out-summary data/private/2027_pilot_arm_summary.csv
```

1株につき複数花の情報を保持するが、`pp_genets_with_followup` と `pp_genets_with_seed` を**独立した株レベル**で数える。担当研究者は Intention-to-treat（割付けに基づく追跡）と、per-protocol（割付けどおり実施・袋と花粉品質を確認）の両方の集計を保持する。後者だけで適応度の因果比較を行うと、バッグ破れ・観測脱落に依存した選択バイアスが出る。

- `completed_itt_flowers`：割付けられ、成熟結果が記録された花（処理変更なども含み得る）。
- `completed_pp_flowers`：正しい処理・供与親・バッグ品質・花粉品質が確認され、成熟結果が記録された花。
- `lost_or_pending` と `missing_execution`：欠測・操作不実施を**結実率の分母に混ぜない**。ただし流失が果実形成に依存する可能性があるので、両処理間でその差も報告する。
- `pp_seed_bearing_flowers`：成熟種子が1以上の花。果実が成熟しても種子が0の花と区別。
- `pp_germination_assessed_flowers`：発芽判定が実際に記入された花。未測定を0件発芽としない。

統計モデルは今後の本試験で株・集団の階層性を含める。**本スクリプトは品質監査と記述的集計のみ**であり、得られたp値や因果効果を主張しない。初年度は管理された実験操作が成立したか、HとF双方に陽性種子があるかを明確にすることが優先。


## ステップ4：花粉管・胚の破壊採取（別花の非縦断コホート）

成熟果実の追跡用に割り付けた花とは**異なる花ID**で、柱頭・花粉管・胚を採取する。[段階別空テンプレート](../data/templates/2027_destructive_assay_header.csv) の `assay_stage` は `stigma` / `pollen_tube` / `embryo`。

| assay_stage | 必須測定値 | 解釈上の注意 |
|---|---|---|
| `stigma` | `pollen_grains_count`、`pollen_grains_germinated` | 花粉付着と発芽。顕微鏡で自家／他家を区別しない限り、親由来は不明 |
| `pollen_tube` | `pollen_grains_count`、`pollen_tubes_at_ovary` | 花粉管が子房へ伸長した状態。胚珠への受精・子の父性は未証明 |
| `embryo` | `ovules_examined`、`fertilized_ovules`、`intact_embryos`、`aborted_embryos` | 組織染色・胚形成の根拠を標本IDで管理し、近交弱勢と後期SIを単独で決めない |

**同じ花を三段階で測ったと解釈しない。** 違う花・違う時刻・時には違う株に由来する割合の単純な積から「受粉→受精→種子の生存確率」を作らない。採取段階ごとに別の母株数、花数、品質不良と欠測数を報告する。

```bash
python3 analysis/audit_2027_destructive_assays.py --self-test

# 現地データが本当に得られてから実行。破壊採取と成熟果実追跡の花IDを照合
python3 analysis/audit_2027_destructive_assays.py \
  --assays data/private/2027_dissection_stages.csv \
  --mature-assignments data/private/2027_assignments.csv \
  --output data/private/2027_stage_assay_summary.csv
```

**保全予算の全コホート照合**：破壊採取で用いる花は元の蕾在庫台帳にも記録済みで、成熟果実用の割付け花とはIDが異なる必要がある。初期設定では**成熟果実用に割り付けた花＋追加の破壊採取花の合計**が、その株の確認済み飛散前蕾の最大50%かつ最大10花を超えるとエラーで停止する。第一波で既に50%を操作対象にした株は同じ台帳から破壊採取を追加できない。別の株/年/事前に承認した保全方針に基づきデザインを再構成すること（**自動で上限を緩和しない**）。

当コードは `genet`・`flower` の重複、花粉数を超える発芽、見かけの花粉粒数を超える子房花粉管、胚珠数を超える受精数、胚数の非整合、欠測を0と記録した誤り、人工授粉の花粉供与株取り違えを検出する。**自家不和合性・胚死亡の原因を自動推定しない**。原画像・顕微鏡記録の照合と統計推定はその後の作業。

## 全体の判定ゲート
- **Gate 0**：対象がオニシバリと独立に確認され、形態的H/Fが区別されたか
- **Gate 1**：予定処理が実施され、供与花粉の品質・バッグ・追跡が妥当か
- **Gate 2**：Fの他殖陽性対照が成功したか
- **Gate 3**：Hに人工他殖由来の成熟種子があったか、Hの人工自殖と袋のみで差があったか
- **Gate 4**：花粉管・胚・成熟種子・発芽、可能なら子の父性を検証し、自家不和合性と近交弱勢を切り分けられたか

現時点では**Gate 0–4に実際のオニシバリ野外データはない**。今回追加したものはGate 1–3の監査プログラムと、未観測を成功・失敗に数えないルールである。
