# HSat2-selective exact-match candidates

広いhapを対象とするASOの標的候補として、STEP 2の既存カウントを再評価しました。最終ASO配列・化学修飾・KD効率を決定したものではありません。STEP 1/2の元ファイルは変更していません。

## 選抜

出発集合は既存のpooled_and_broad90_E_gt_10候補9,740配列です。573非参照hapで次の条件を評価しました。CHM13は合算から除外しています。

1. 合算HSat2 / (HSat2 + HSat3) count >= 0.999。
2. 合算log2((HSat2 density + 1e-9)/(HSat3 density + 1e-9)) > 10。
3. 合算log2((HSat2 density + 1e-9)/(background density + 1e-9)) > 10。
4. HSat2に存在するhap数 >= 516。

以上で2,484配列。全2,484配列が、516 hap以上で個別にもHSat2 count >=10かつHSat2対背景E >10を満たしました。

さらに、同じ516 hap以上でHSat2対背景の条件、HSat2対HSat3 E >10、HSat2 count fraction >=0.999をすべて同時に満たすものは1,797配列でした。こちらを厳しい候補集合として保存しました。各hapのHSat3有効開始位置数が0で密度を計算できない場合は、個別family条件を合格とみなしません。

背景は既存HSat2/3 BEDの外側です。densityは領域ごとの有効16-mer開始位置数で正規化し、補正値はper-millionに換算する前の密度に加えます。今回のHSat2対背景Eは、従来のHSat2/3 union対背景Eとは別に計算しています。

## 出力

- `HSat2_selective_pooled.tsv`: 2,484配列。
- `HSat2_selective_broad90.tsv`: 個別HSat2 enrichmentも516 hap以上で合格する2,484配列。
- `HSat2_selective_strict_broad90.tsv`: 個別family特異性も含めて516 hap以上で合格する1,797配列。
- `summary.json`: 閾値、件数、入力ハッシュ、検証結果。

並び順は合算HSat2対背景Eの降順です。背景の絶対カウント・HSat2カウント分位点・GC%も併記しています。配列類似性による重複整理やASO化学修飾に応じたスコアは適用していません。近接・重複する16-merを含み、上位配列同士が独立した標的とは限りません。

配列は5′→3′のcanonical DNA表記と、そのreverse complementです。両鎖のDNAカウントを合算しているため、この表だけからRNAを標的にするASOの向きは確定できません。転写鎖・発現、1–2 mismatch、他の転写産物への一致を確認して最終候補を選びます。HSat3カウント0は今回の注釈領域における完全一致0であり、交差作用がないことの証明ではありません。

再実行:

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/select_hsat2_candidates.py
```

参考: HSATIIを標的にしたLNAでRNA量の増加が観察された報告もあります。DNA配列の濃縮とRNAのKD効果を分けて評価する必要があります。Porter et al., JCI (2022), https://www.jci.org/articles/view/155931 (Figure 7 / Results)。この研究は今回抽出した配列の効果を検証したものではありません。
