# Satellite ASO maker

既存の VallePrep/Altemose HSat2・HSat3 注釈と T2T 通過染色体から、HSat2/3 領域に濃縮する canonical 16-mer を見つけ、hap・family・染色体別に数え、配列パターン群と 0/1/2 mismatch によって ASO 標的候補を優先順位付けした解析です。最後に、評価した **72 配列**について CHM13 上の全 exact ヒット座標、CenSat annotation の内外、hap 群別・染色体別の分布を可視化しました。この README は同じ入力と実行環境で解析を再現するための入口です。実装の細部は [STEP 0](step0/README.md)、[STEP 1](step1/README.md)、[閾値スイープ](step1/threshold_sweep/README.md)、[STEP 2](step2/README.md)、[STEP 3](step3/README.md)、[STEP 4](step4/README.md)、[72 配列の位置解析](step4/all72_landscape/README.md)を参照してください。

**保存済み結果と新規ワークフローの区別:** このページの 159,988 配列・72 配列などの数値と各 STEP の手順は、旧設定（pooled count ≥100 または 1 hap count ≥10、後段の threshold sweep、複数染色体の共通 core を合算した主クラスタリング）による保存済み解析を記録したものです。現在の [Snakemake ワークフロー](workflow/README.md)は **pooled count ≥10 または 1 hap count ≥500、`E > 10`、sweep なし、染色体ごとに独立した STEP 3 クラスタリング**をデフォルトとします。新条件の件数は再計算するまで不明です。旧数値を再計算する場合は、変更前の GitHub commit `11fe4a5` のコードを使用してください。

## 解析の全体像

| 段階 | 処理 | 今回の結果・次段階への入力 |
| --- | --- | --- |
| STEP 0 | 入力の検証、HSat2/3 と背景の BED、有効 16-mer 開始位置数 | 574 sequence set（573 hap + CHM13）、9 set を除外 |
| STEP 1 | Jellyfish で exact canonical 16-mer を発見、HSat23 vs 背景の濃縮度 | 予備候補 1,195,753、元の `E >= 2` 採用 1,131,847 |
| 閾値スイープ | 同じ count で `E > 2, 4, ..., 20` を比較 | 後続に `E > 10` の pooled または hap 条件を満たす **159,988** 配列を採用 |
| STEP 2 | 159,988 配列の family × chromosome × hap exact count | Matrix A/B/C、family・染色体特異性、hap prevalence |
| STEP 3 | 共通領域の regional pattern クラスタリング、代表 24 配列の 0–2 mismatch | HSat2 は 3 群、HSat3 は 7 群 |
| STEP 4 | 31,335 exact 候補を統合、72 代表配列の 0–2 mismatch と順位付け | category robustness: A=4、B=3、C=63、D=2 |
| 位置解析 | 72 配列の CHM13 exact ヒットと CenSat 内外、hap 群・染色体別図 | CHM13 に配列別 559,297 ヒット、広い CenSat 外 55 ヒット（21 配列） |

候補とスコアは**ゲノム DNA の配列一致・注釈領域内の分布**を評価したものです。RNA 発現、転写鎖、ASO の化学修飾・送達・KD 効率は測定していません。

## 作業ディレクトリ、入力、環境

以下のコマンドはこの `sat-aso-maker` ディレクトリをカレントディレクトリとして実行します。スクリプトは `scripts/` の親を project root、入力元をその隣の `VallePrep_v0.0.0/` として解決します。採用・除外判定と実際に使用した入力の**絶対パス**は [step0/manifest.tsv](step0/manifest.tsv)に記録しています。別の計算機へ移す場合は同じ配置を再現するか、スクリプト内のパス設定を合わせてください。

`../VallePrep_v0.0.0/results/data/{sample}/` から hap 名を `{hap}` として以下を使用します。

| ファイル | 役割 |
| --- | --- |
| `centro/{hap}.t2t.pass.bed` | 解析対象の T2T 通過染色体 |
| `filfa/{hap}.t2t.filtered.fasta.gz` | 座標・向きの基準となる filtered/oriented FASTA |
| `filfa/{hap}.t2t.filtered.fasta.gz.fai` | FASTA index |
| `filfa/{hap}.t2t.filtered.fasta.gz.gzi` | 圧縮 FASTA のランダムアクセス用 index（元データにある場合） |
| `hsat23/{hap}.HSat2and3_Regions.lifted.bed` | 既存 Altemose HSat2/3 領域。先頭 3 列を座標として使用 |
| `hsat23/{hap}.HSat2and3_Regions.HSat2and3.log` | 注釈処理の完了確認 |

CHM13 の位置解析には追加で `../VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed`（**広い CenSat 内外**の判定）と `../VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed`（alphaSat/HSat2/HSat3/ct の細かい注釈）を使用します。初期調査は [valleprep-hsat23-inventory/](valleprep-hsat23-inventory/) にありますが、後続解析の採用集合は STEP 0 の manifest が正です。

今回の実行環境は Python 3.13（NumPy 等を含む `/home/senescence/miniconda3/envs/kmer_validation_env/bin/python`）、同環境の `samtools`、ローカル Jellyfish 2.3.1（`.local/bin/jellyfish` と exact prefix filter 版 `.local/bin/jellyfish-prefix`）、`g++`、`Rscript` です。スクリプトにはこれらのパスを直接参照するものがあります。C++ scanner のソースと実行ファイルは `scripts/` に保存済みです。必要なら次でビルドできます。STEP 0 はゲノムを複製せず元 FASTA/index への symlink を作るため、**元ファイルが同じ場所で読めること**が再現条件です。

```bash
g++ -O3 -std=c++17 scripts/step0_scan.cpp -lz -o scripts/step0_scan
g++ -O3 -std=c++17 scripts/step1_pack.cpp -lz -o scripts/step1_pack
g++ -O3 -std=c++17 scripts/step2_count.cpp -o scripts/step2_count
g++ -O3 -std=c++17 scripts/step3_count.cpp -o scripts/step3_count
```

Jellyfish の prefix filter 版は [patch](scripts/jellyfish-prefix12.patch)を適用したローカルビルドです。元の Jellyfish と exact count が一致することを検証しています。ビルド済みバイナリを別環境へ移せない場合は同じ Jellyfish 2.3.1 と patch から再ビルドし、[STEP 1 の検証](step1/README.md)を実行してください。

既存結果を閲覧するだけなら再計算は不要です。各段階の `config.json`、`validation_summary.json`、`COMPLETE.json`（該当する段階）に設定と検証結果を保存しています。checkpoint は固定入力を前提に再利用されます。**入力・スクリプト・閾値を変更して解析し直す場合は、既存 checkpoint と混ぜず新しい出力先を用意**してください。574 set の FASTA 走査、Jellyfish 集計、mismatch 再走査は重い処理です。

## STEP 0: 領域の準備と有効窓数

`step0_preprocess.py` は T2T pass BED を起点に FASTA/index、Altemose BED、完了 log を探し、座標の整合性を検査します。完了していない注釈や不正な座標を持つ hap は 0 件として扱わず**set 全体を除外**します。今回 583 入力中 574 が採用、9 が除外されました。理由は [invalid_annotations.tsv](step0/invalid_annotations.tsv)、[manifest.tsv](step0/manifest.tsv)、[chromosome_status.tsv](step0/chromosome_status.tsv)で追えます。

同一 family 内の重複区間を統合し、family 間の重なりを `ambiguous` に分けます。`HSat23` は HSat2/3 の包含 union、`background` は T2T 通過染色体上の `HSat23` の補集合です。BED は `step0/regions/{hap}/{HSat2,HSat3,HSat23,ambiguous,background}.bed`、参照 FASTA のリンクは `step0/sequences/{hap}/source.fasta.gz` にあります。座標は **0-based・半開区間**で filtered/oriented FASTA が基準です。別 interval は連結せず、非 ACGT も窓を切ります。有効開始位置数は各連続 ACGT run の `max(0, 長さ - 15)` の和です。union には family 境界をまたぐ窓が入り得るので、HSat2 と HSat3 の k-mer count の和と一致するとは限りません。

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step0.py
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step0_preprocess.py --workers 4
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/verify_step0.py
```

主な出力は [region_stats.tsv](step0/region_stats.tsv)（hap × 染色体 × 領域の bp・有効窓数）、[sequence_views.tsv](step0/sequence_views.tsv)、[validation_summary.json](step0/validation_summary.json)です。今回の 8,531 set×chromosome と 42,655 region 行の分割・総塩基数が検証済みです。STEP 0 は k-mer を選びません。`background` は**既存 HSat2/3 注釈の外側**であり、全 CenSat の外側や衛星配列を含まない領域とは同義ではありません。

## STEP 1: canonical 16-mer 発見と濃縮閾値

Jellyfish `-m 16 -C` で HSat23 の全 canonical 16-mer を数え、reverse complement を同一配列へまとめます。非参照 573 hap の pooled HSat23 count ≥100、または少なくとも 1 hap の HSat23 count ≥10 を予備候補にし、HSat2、HSat3、HSat23、background を両方向で再カウントします。background は候補とその reverse complement に限定して検索するため、background 全種の表ではありません。CHM13 は単独で集計しますが、pooled 候補選定・prevalence の分母から除外します。

`target_density = target_count / target_valid_16mer_starts`、`background_density = background_count / background_valid_16mer_starts` として、濃縮度を次で計算します。

```text
E = log2((target_density + 1e-9) / (background_density + 1e-9))
```

元の STEP 1 は pooled `E >= 2` かつ count ≥100、**または**少なくとも 1 非参照 hap の `E >= 2` かつ count ≥10 を採用しました。この探索段階の出力は保持しています。後から同じ正確な count で `E > 2, 4, 6, ..., 20` を比較し、下流用に **`E > 10`** を選びました。下流集合は pooled `E > 10` かつ pooled count ≥100、**または**少なくとも 1 hap の個別 `E > 10` かつその hap の count ≥10 の和集合です。`>` は厳密な不等号で、元の `>= 2` と異なります。`E > 10` は補正密度比が 1,024 倍を超えることを意味しますが、それだけで off-target が少ないと決める条件ではありません。

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step1.py
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step1_pipeline.py --phase discovery --workers 8
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step1_pipeline.py --phase counts --workers 8
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/finish_step1.py
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/sweep_step1.py
Rscript scripts/plot_step1_sweep.R
```

`finish_step1.py` は score、plot、最終検証を実行します。結果は [config.json](step1/config.json)、`step1/candidate_kmers.tsv`、`kmer_catalog.tsv.gz`、`enrichment_summary.tsv.gz`、`counts/{hap}.npz`、[COMPLETE.json](step1/COMPLETE.json)です。後続への**確定入力**は [candidates_E_gt_10.tsv](step1/threshold_sweep/candidates_E_gt_10.tsv)の 159,988 配列です。pooled pass は 34,332 配列で、このうち 516/573 hap 以上が個別条件も通過する broad90 は 9,740 配列です。少数 hap でのみ強く濃縮する配列も和集合に残ります。[threshold_summary.tsv](step1/threshold_sweep/threshold_summary.tsv)と [threshold_sweep.png](step1/threshold_sweep/threshold_sweep.png)で閾値感度を確認できます。

`kmer_id` は `k16_` + canonical 16-mer を 2-bit（A=0,C=1,G=2,T=3）で符号化した 8 桁 hex です。canonical 配列の向きは転写鎖や ASO の発注方向ではありません。Jellyfish count が STEP 0 の有効窓数と合うこと、両方向 count の和が canonical count と合うことを検証します。高速化と counter 精度の詳細は [STEP 1 README](step1/README.md)を参照してください。

## STEP 2: family × chromosome × hap の正確な count

`E > 10` 和集合の 159,988 配列を、574 set の T2T 通過染色体ごとに HSat2、HSat3、HSat23 union、ambiguous BED 内で完全一致カウントします。逆相補鎖を合算した canonical count と方向別 count を保持します。αSat 専用の count matrix はありません。

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step2.py
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step2_pipeline.py --workers 6
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/finish_step2.py
```

`finish_step2.py` は matrix、指標、PNG、最終検証を作ります。

| 出力 | 単位・用途 |
| --- | --- |
| [candidates.tsv](step2/candidates.tsv)、[haplotypes.tsv](step2/haplotypes.tsv) | 候補行順、hap 列順、参照/非参照の対応 |
| `matrices/matrix_A_counts.npy` | 159,988 × 574、HSat23 union の hap 合計 |
| `matrices/matrix_B_counts.npy` | 159,988 × 1,148、hap × HSat2/HSat3 合計 |
| `counts/{hap}.npz` と `matrices/matrix_C_columns.tsv` | 染色体 × family × 方向の Matrix C。単一 TSV には展開していない |
| `matrices/availability_and_denominators.npz` | 染色体の可用性と有効 16-mer 開始位置数 |
| `candidate_metrics.tsv.gz`、`chromosome_specificity.tsv.gz` | family specificity、hap prevalence、染色体 concentration / recurrence |
| `plots/*.png`、[COMPLETE.json](step2/COMPLETE.json) | 概観図、最終検証 |

各 count の density per million は `count / valid_starts × 1e6` です。染色体が T2T を通らず**未測定**なら availability=0、測定した領域で k-mer が見つからない場合のみ 0 count です。測定済みでも family BED がなく分母 0 の場合の density は NA です。Matrix A/B は測定可能な染色体上の合計なので、可用性が異なる hap 間で raw count だけを比較しないでください。STEP 2 の chromosome specificity は注釈された HSat2/3 **内部**の分布で、ゲノム全域の結合特異性ではありません。全 hap の染色体別 count 和は STEP 1 の両方向 count と照合済みです。

## STEP 3: regional pattern 群と 24 代表配列の mismatch

family 別に、非参照 hap の 85% 以上で当該 family の有効窓を持つ染色体を共通領域とします。HSat2 は chr7+chr10 が揃う 489 hap、HSat3 は chr5+chr7+chr10+chr20 が揃う 442 hap を主解析に使用しました。共通領域の count を有効窓数で正規化して `log1p` 変換し、変動の大きい候補から相関 `|r| >= 0.95` の冗長性を落とした代表特徴量を作ります。R の Ward.D2 階層クラスタリングで k=2–8、silhouette と特徴量再抽出 20 回の安定性を比較しました。採用条件は平均 silhouette ≥0.25、最小群 ≥10、ARI 中央値 ≥0.75、最小対応 Jaccard 中央値 ≥0.75 です。結果は HSat2 3 群、HSat3 7 群です。

この群は**共通染色体での配列パターン**です。全利用可能染色体で再解析すると群分けは大きく変わるため、全ゲノム・祖先集団の固定的な hap group と解釈しません。共通領域が欠ける hap は `NA`（未割当）です。CHM13 は学習から除外して別に投影しています。

全候補について pan、group、chromosome の exact 条件を比較し、冗長性を下げた **24 代表配列**だけを全 574 set で再走査して Hamming 距離 0、ちょうど 1、ちょうど 2 のヒット数を別々に記録します。1 配列の半径 2 近傍は最大 1,129 種で、reverse complement を含む重複開始位置は query ごとに一度だけ数えます。indel は含みません。pan の厳密判定は target vs background `E > 10`、target vs other family `E > 10`、target family 割合 ≥99.9%、target count ≥10 を pooled と 516/573 hap 以上で要求します。group / chromosome 条件の詳細は [STEP 3 README](step3/README.md)を参照してください。

```bash
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step3_prepare.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step3.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/run_step3.py
```

`run_step3.py` が clustering、exact 候補分類、近傍列挙、pilot、本計算、集計、最終検証を順に実行します。出力の起点は [hap_groups.tsv](step3/hap_groups.tsv)、`clustering/`、`candidates/`、`mismatch/pooled_mismatch_scores.tsv`、`mismatch/candidate_robustness.tsv`、[COMPLETE.json](step3/COMPLETE.json)です。24 配列以外の 0–2 mismatch 性質は、この段階では未評価です。

## STEP 4: 用途別 ASO 標的候補ランキング

STEP 1–3 から重複を除いた 31,335 exact 候補を統合し、`pan`、`regional_group`、`chromosome` の用途別に順位付けします。exact composite score は target abundance 20%、family specificity 20%、背景対比 E 20%、hap coverage 15%、背景の絶対数 10%、GC・homopolymer・base entropy・自己相補性 proxy 15% の family 内 percentile 合成です。物理化学的な ASO 効力予測ではありません。多様性を考慮した shortlist から STEP 3 の 24 配列を再利用し、新規 48 配列を加えた **72 配列**について全 574 set で 0/1/2 mismatch を評価します。

```bash
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_rank.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_neighbors.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_mismatch.py --workers 8
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_finalize.py
```

最初に見るファイルは、31,335 件の `ranked/all_exact_ranked.tsv.gz`、72 件の [ranked/mismatch_evaluated_ranked.tsv](step4/ranked/mismatch_evaluated_ranked.tsv)、用途別 `shortlists/{HSat2,HSat3}_{pan,regional_group,chromosome}_candidates.tsv`、半径別 `mismatch/pooled_mismatch_scores.tsv`、[COMPLETE.json](step4/COMPLETE.json)です。`category_robustness_tier` は用途別 screen を距離 ≤2 まで通過する A、≤1 までの B、exact のみの C、exact 用途別 screen も通らない D です。`strict_pan_robustness_tier` は pan のより厳密な別判定なので混同しないでください。`canonical_target_5to3` はゲノム標的の表記です。最終表は両転写方向に対応する ASO 候補を併記しており、発注方向には転写鎖の確認が必要です。

## CHM13 上の位置と CenSat 内外: 72 配列の図

まず上位 pan-HSat3 の 3 配列で位置図と CenSat annotation を分離した図を作り、[top 3 の方法と結果](step4/top_kmer_landscape/README.md)に保存しました。その後、STEP 4 の 72 配列全件へ拡張しました。72 配列を見るときは [図の索引](step4/all72_landscape/PLOT_INDEX.md)、[CenSat 外ヒット一覧図](step4/all72_landscape/plots/all72_censat_outside_overview.png)、[72 配列の数値要約](step4/all72_landscape/chm13_candidate_censat_summary.tsv)から始められます。

CHM13 の 24 filtered chromosome を直接走査し、canonical 配列と reverse complement の**完全一致**開始位置を記録します。広い CenSat BED に 16-mer 区間が完全に含まれるかを判定し、細かい family BED の重なりには `HSat3 → HSat2 → alphaSat → ct` の優先順で排他的なラベルを付けます。広い CenSat 内外判定と細かい family ラベルは独立です。個別 `plots/{kmer_id}_landscape.png` は、上段が全染色体の CenSat 境界・family annotation・100 kb bin の exact ヒット、中段が target family BED 内の hap 群別 count、下段が群 × 染色体の中央値です。`plots/{kmer_id}_dominant_zoom.png` は最多ヒット染色体の拡大図です。両図とも 72 枚あります。CHM13 に exact ヒットのない 2 配列も図に明記されています。

```bash
g++ -O3 -std=c++17 scripts/scan72_chm13.cpp -lz -o scripts/scan72_chm13
scripts/scan72_chm13 step4/mismatch/queries.tsv step0/sequences/chm13v2.0/source.fasta.gz step4/all72_landscape/raw/chm13_all72_exact_hits.tsv.gz
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/all72_landscape.py
Rscript scripts/plot_all72_landscape.R
Rscript scripts/plot_all72_zoom.R
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/finalize_all72_landscape.py
```

座標の全行は `step4/all72_landscape/chm13_all72_exact_hits_with_censat.tsv.gz`（`start_0based`/`end_0based` は 0-based 半開）、染色体別の内外数は [chm13_censat_by_chromosome.tsv](step4/all72_landscape/chm13_censat_by_chromosome.tsv)、hap × 染色体の数は `hap_chromosome_exact_counts.tsv.gz`、群別要約は [group_total_summary.tsv](step4/all72_landscape/group_total_summary.tsv)と [group_chromosome_summary.tsv](step4/all72_landscape/group_chromosome_summary.tsv)にあります。CHM13 exact ヒットは配列別に合計 559,297 行、うち広い CenSat 外は 55 行（21 配列）でした。**55 は配列別ヒット行の合計で、55 個の異なる座標や hap 全体の CenSat 外数ではありません。** hap 側の `target_family_exact_hits` は対象 HSat2/3 BED 内で完全一致した数で、hap 全ゲノムのヒット数ではありません。CHM13 座標図の内外は 0 mismatch のみで、1/2 mismatch の位置別 off-target はこの図からは分かりません。

## 再現性の確認と解釈上の注意

Snakemake で新条件の独立した解析を実行する場合は [ワークフローの設定・実行方法](workflow/README.md)を参照してください。[設定例](workflow/config.example.yaml)をローカルの `workflow/config.yaml` にコピーし、run ID、入力元、並列数、`E >` 選択閾値とカウント下限を調整できます。初めに `--dry-run` で依存関係を確認し、実行結果は `runs/<run_id>/` に保存します。新ワークフローでは sweep 工程はなく、STEP 2 は STEP 1 の `candidate_kmers.tsv` を直接受け取ります。

各段階の `COMPLETE.json` は当時の全処理完了・照合結果です。STEP 0 は領域分割と窓数、STEP 1 は Jellyfish count と STEP 0 の窓数・向き、STEP 2 は全候補の染色体和と STEP 1、STEP 3/4 は距離 0 の全 574 set count と STEP 1、および mismatch 近傍の網羅性を検証しました。72 配列の [COMPLETE.json](step4/all72_landscape/COMPLETE.json)は CHM13 の fine HSat2/3 位置 count と STEP 4 の距離 0 行列、先行 top 3 の座標、72 枚ずつの PNG、55 件の CenSat 外 count を監査します。ファイル内容やコードを変えた場合、古い `COMPLETE.json` が新しい実行の保証にはなりません。

本解析の `background` は Altemose の HSat2/3 annotation の補集合で、dna-brnn などの事前選択領域外も含みます。CenSat 全域を正確に引いた領域ではありません。CHM13 の「CenSat 外」は別の広い CenSat BED で判定しています。HSat2/3 注釈が不完全なら family 特異性・背景濃縮も変わり得ます。欠測染色体と 0 count、count と density、HSat23 union と二つの family の和、exact と mismatch、CHM13 全ゲノム座標と hap の family 内 count を区別してください。αSat は STEP 0–4 の独立した候補 family ではなく、CHM13 位置図の fine annotation にだけ現れます。
