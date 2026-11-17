# v0.42–v0.43｜2010年日本産Daphne DNAバーコーディング配列の所在・原著照合

Date: 2026-11-17 JST。**歴史的な原配列の探索結果**。新たなオニシバリ核系統樹、交雑検定、繁殖データは生成されていない。

## 2010年に元研究チームが本当に測定したマーカー

Kawahara et al. (2010), 「樹木分類群及び産地識別におけるバーコード情報の可能性－ジンチョウゲ属の実例から－」第121回日本森林学会大会, J-STAGE DOI https://doi.org/10.11519/jfsc.121.0.332.0 （発表B08、2010年3月）。

原著の方法にあるマーカー名を確認：
- **葉緑体の6領域**：rbcL, matK, trnH–psbA, psaB–rpl4, rpl14–rps16, trnL–trnF
- **核リボソーム反復配列1領域**：ITS2

これを「独立した核遺伝子7座位」や「葉緑体の独立な6つの遺伝子樹」として扱うのは誤り。6葉緑体領域は同じオルガネラ系譜で連鎖し、ITS2は多コピー反復領域に属する。独立低コピー核遺伝子座は、この発表のマーカー記述だけでは**0**。公開原文には比較対象の候補マーカー結果が書かれているが、GenBankアクセッション・個体別voucherの確実な対応表は記載されていない。

別途、2011年科研費20380093成果報告書 https://kaken.nii.ac.jp/file/KAKENHI-PROJECT-20380093/20380093seika.pdf のp.2–3には、オニシバリ・ナニワズ等4種の雌性両全異株判定と自家和合性の低い分類があるが、花／株ごとの人工授粉・成熟子の表は載らず、研究成果の雑誌記事は査読なしの研究紹介・一般論2件、学会発表7件。2010年の分類学会「日本産ジンチョウゲ属の分子系統とフェノロジー進化」もその学会発表のひとつ。**公開書誌で見つかった発表・マーカー名を実際のFASTAアーカイブが得られたことと混同しない**。

## 2010年の元配列をNCBIから回収できるか：直接の索引検索

NCBI ESearch+ESummary (nuccore), 年次2008–2013（チーム検索は2008–2015も）。元著者名やDaphne分類名を厳密に含む索引検索と、同時期のマーカー別広域検索を別々に実行：

| 検索対象 | NCBI索引ヒット数 | 厳密な結論の上限 |
|---|---:|---|
| Daphne×Kawahara著者×2008–2013 | **0** | この年・著者名の索引に該当なし |
| Daphne×(Kawahara/Yoshimura/Yamashita)著者×2008–2015 | **0** | チーム名で登録した配列を確認できず |
| P/jezoensis/kamtschatica名×2008–2013 | **0** | この学名範囲の索引ヒットはないが、旧学名・種名混同・登録年外の可能性 |
| Daphne×葉緑体スペーサー名×2008–2013 | **2** | 種属とマーカー名の汎用一致のみ |
| Daphne×ITS2×2008–2013 | **60** | 実際の発表研究のアクセッションではなく、属・マーカーの一般的な検索結果 |

5つの検索は同義でも独立な研究でもなく、索引のクエリ・登録日・分類名に依存する。この「0」は**元データが存在しない**と示すものではない。掲載年と登録年のズレ、登録記録に著者名がないケース、他の研究機関での提出、ITS2の古い分類名、研究室保管の配列は残る。

## v0.43: 14個の候補GenBankレコードを実際の記録本文で確認

広域検索で得た標本の最初の計14件をNCBI EFetch GenBankフラットファイルで読み、記載された原著REFERENCE著者・タイトル・分類名・収集者を照合した。

- 葉緑体候補 **GU736360 (D. tangutica)** と **GU736359 (D. mezereum)** は、最初の引用文献が **Zhang, Volis & Sun** の *Stellera chamaejasme* チベット高原系統研究に対応し、Kawahara・Yamashitaらの2010年日本産Daphne研究ではない。
- ITS2候補 **GQ167547 (D. rodriguezii)**、**GQ167544 (D. oleoides)**、および GQ167536、GQ167535、GQ167534、GQ167533、GQ167532、GQ167525、GQ167524、GQ167523、GQ167522、GQ167521（D. laureola）は、最初の文献が **Alonso & Herrera**, *Back-and-forth hermaphroditism: phylogenetic context of reproductive system evolution in subdioecious Daphne laureola*。やはり日本の2010年Daphneバーコード発表の元配列ではない。

検査した14件のうち、**元プロジェクトの研究者名が最初の登録文献の著者情報に一致したものは0件**。残る検索に載った未検査のITS2記録や登録年をずらしたレコードが、元データと無関係だと断言するものではない。

この監査は**「関係ない欧州Daphne/中国Stellera系統の配列を日本産オニシバリの未発表原配列として誤使用する」ことを防ぐため**のもので、系統機構を新しく解明したものではない。

## 2025年の核・葉緑体不一致との接続

Han et al. (2025) DOI https://doi.org/10.11110/kjpt.2025.55.3.157 の P=日本産 D. pseudomezereum、C=韓国産旧var. koreana、K=D. kamtschatica では、葉緑体CKとnrITS PCの一貫したマーカー系譜の違いが再現済み。しかし既存のC葉緑体 ON244034 (2023 Gangneung) とC ITS PX056848 (2025 Pyeongchang) は**別個体**であり、低コピー核遺伝子座の親子対応はまだない。2010年の発表マーカーが将来回収されても、核ITS2+連鎖葉緑体だけでILSと古い浸透性交雑を分離するのは不可能。

次の最優先の回収候補：
1. Kawaharaら2010年の実験元の**ITS2、葉緑体6領域のFASTA/クロマトグラム**と、voucher、集団、個体、元の種同定対応。
2. 2008–2010年科研費の**オニシバリ雄的花・雌花の自家/他家手動授粉対照別の花数・成熟果実・種子・発芽・花粉管データ**（元報告には数表がない）。
3. Hanらの2025年研究とYooら2023年研究の**P/C/Kから同じvoucherの未選別全DNA読み取り**、および不足する独立低コピー核遺伝子座。既に依頼文の草案は作成したが**未送信**。

## 自動再現

- scripts/daphne_2010_barcode_genbank_recovery_v42.py
- scripts/daphne_2010_genbank_provenance_v43.py
- Source v42 CI SUCCESS https://github.com/zuizui0223/onishibari/actions/runs/44011060026
- Source v43 CI SUCCESS https://github.com/zuizui0223/onishibari/actions/runs/44011233441
- JSON artifacts: analysis_results_daphne_2010_barcode_archive_v42/daphne_2010_barcode_genbank_recovery_index_v42.json; analysis_results_daphne_2010_barcode_archive_provenance_v43/historical_daphne_2010_barcode_record_provenance_v43.json

検索索引の「0」と同一論文の複数マーカーを独立研究の実験数へ追加しない。元配列の新取得・分子系統の新推定・オニシバリの自家不和合性実験は未実施。Draft PR#2のみ、main不変更。
