# STEP 1 enrichment threshold sweep

既存の正確な16-merカウントを再利用し、E > 2, 4, 6, 8, 10, 12, 14, 16, 18, 20を比較しました。CHM13を除く573 hapが対象です。元のSTEP 1の設定・出力は保存しています。

E = log2((target_density + 1e-9) / (background_density + 1e-9))。targetはHSat2/3のunion、背景は対象T2T染色体上の既存HSat2/3 BED外側です。densityの分母は有効16-mer開始位置数です。全CenSatの除外を保証する背景ではありません。

- pooled: 合算target count >=100、合算E > 閾値。
- hap_only: 少なくとも1 hapでtarget count >=10かつE > 閾値、pooled条件は不成立。
- union: pooledまたはhap条件を満たす候補。
- broad90: 516/573 hap以上で個別にtarget count >=10かつE > 閾値。単に配列が存在するprevalenceとは異なります。対象染色体の利用可能性の影響を受けます。

## 結果と使い方

E > 10のpooled候補は34,332配列（E > 2の1,095,175配列から96.9%減）。そのうちbroad90を満たす9,740配列を、広いhapを対象にする後続解析の出発点として保存しました。E > 12ではpooled 12,111配列、pooledかつbroad90 3,101配列となります。最適なASO閾値はこの解析だけでは決定できず、family/染色体特異性、ミスマッチ、KD評価が必要です。

Eは相対密度比です。閾値を上げても背景の絶対カウントが減るとは限りません。E > 10のpooled候補の背景カウント中央値は492、背景ゼロは2,931配列です。またepsilonにより、背景ゼロでもE > 10にはtarget density > 1.023e-6が必要です。このため高閾値ではtarget abundanceの強い選抜も起こります。

## ファイル

- threshold_summary.tsv: 全閾値・経路の候補数、背景ゼロ、prevalence、カウント中央値。
- threshold_sweep.png: 比較図。
- candidates_E_gt_8.tsv / candidates_E_gt_10.tsv / candidates_E_gt_12.tsv: union候補。pooled_pass/hap_pass列で判別できます。
- candidates_E_gt_10_pooled.tsv: 全hap合算でE > 10の34,332配列。
- candidates_E_gt_10_pooled_broad90.tsv: 上記のうち516 hap以上で個別にも条件を満たす9,740配列。
- local_threshold_counts.npz: 各閾値で個別条件を満たすhap数。codesはSTEP 1のpreliminary candidatesと同順。
- validation.json: 閾値増加に対するpooled/union候補数と個別合格hap数の単調性検証。

再実行（プロジェクト直下）:

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/sweep_step1.py
Rscript scripts/plot_step1_sweep.R
```

元のcount下限を維持するため、STEP 1 preliminary candidates全体を評価すれば今回の閾値比較は網羅できます。count下限を緩和する比較は含みません。
