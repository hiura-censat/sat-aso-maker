# STEP 2: family × chromosome × haplotype counts

> このディレクトリの159,988候補と以下の行列サイズは旧設定で得た保存済み結果です。現在の Snakemake ワークフローは STEP 1 の `candidate_kmers.tsv` を直接入力し、pooled ≥10 または 1 hap ≥500、各経路 `E > 10` の新条件で候補数を決めます。新しい行列は `runs/<run_id>/step2/` に作成され、行数は再計算まで不明です。

## 対象と入力

STEP 1で全573 hap合算、または少なくとも1つのhapでE > 10を満たした159,988個のcanonical 16-merを集計します。合算で合格した34,332配列、さらに516/573 hap以上で個別にも合格した9,740配列を`candidates.tsv`のフラグで識別できます。CHM13を含む574セットを数え、pooled指標からCHM13を除外します。

入力FASTAは`haplotypes.tsv`に記録した既存ファイルです。BEDは`../step0/regions/{hap}/{HSat2,HSat3,HSat23,ambiguous}.bed`、候補の出所は`../step1/threshold_sweep/candidates_E_gt_10.tsv`です。実体を複製せず、samtools faidxで必要領域のみ読みます。使用する座標はSTEP 0で検証済みのfiltered/oriented FASTA上の座標です。αSatの注釈は今回の対象にありません。

## カウントの定義

16塩基の完全一致を、各BEDレコード内で1塩基ずつずらして数えます。小文字ACGTを受け入れ、非ACGTやレコード境界でリセットします。reverse complementはcanonical配列に合算し、両方向のカウントも保存します。回文配列はforward側で一度だけ数えます。方向はfiltered FASTAに対する方向であり、転写方向やASO鎖を意味しません。

STEP 0の同family内の重複統合とfamily間の曖昧区間分離を引き継ぎます。今回のambiguous領域は0 bpです。HSat23はinclusive unionで、HSat2とHSat3の境界をまたぐ16-merも含むため、HSat2 + HSat3と常に同じにはなりません。その差は`hap_qc.tsv`のjunction countに記録します。Matrix Aの定義を維持するためunionも独立にカウントします。

カウンタは`step2_count.cpp`の正確なrolling 16-mer実装です。カウントはuint64。ランダム配列・小文字・N・レコード境界・回文・両鎖・複数染色体・高カウントのfixtureをPython総当たりと比較します。さらに全hapで染色体別カウントを合算し、STEP 1のJellyfish結果と候補・family・方向ごとに完全一致することを検証します。各染色体の全有効開始位置数はSTEP 0と照合します。

## 行列と並び順

すべての候補行は`candidates.tsv`の`row_index`順です。kmer_idはSTEP 1と共通です。

| 出力 | 形状・列 | 内容 |
| --- | --- | --- |
| `matrices/matrix_A_counts.npy` | 159988 × 574 | HSat23 unionのhap別総カウント |
| `matrices/matrix_A_density_per_M.npy` | 同上 | 上記 / hapのunion有効開始位置数 × 1e6 |
| `matrices/matrix_B_counts.npy` | 159988 × 1148 | hap × familyの総カウント |
| `matrices/matrix_B_density_per_M.npy` | 同上 | 上記 / hap・familyの有効開始位置数 × 1e6 |
| `counts/{hap}.npz` | 159988 × 観測染色体数 × 4 × 2 | Matrix Cのhap別シャード、family別に加えてunionとambiguous、方向別カウントを保持 |

Matrix A列は`haplotypes.tsv`、Matrix B列は`matrices/matrix_B_columns.tsv`に対応します。Matrix Cの論理列は`matrices/matrix_C_columns.tsv`にfamily × chromosome × hapを記録し、シャード内位置に対応付けています。巨大な単一TSVに展開する代わりに圧縮シャードを使用します。

NPZシャードのキー:

- `counts`: uint64、[candidate, chromosome, region, orientation]。
- `chromosomes`: このhapでT2T通過した染色体の順序。
- `regions`: HSat2, HSat3, HSat23, ambiguous。
- `valid_starts`: [chromosome, region]の有効16-mer開始位置数。
- `bp`: [chromosome, region]の領域長。

Matrix Cのcanonical countは最後の方向軸の和、density per millionは `counts.sum(axis=3) / valid_starts[None,:,:] * 1e6` で得られます。分母0の密度はNaNにします。両方向と分母を保持し、密度の重複保存を避けています。

### 0と欠測

対象外染色体はシャードに存在せず、Matrix C列の`available=0`となります。解析した染色体で候補が存在しない場合は0です。染色体を解析したが当該familyの注釈がない場合はcount=0、valid_starts=0、density=NAです。

`matrices/availability_and_denominators.npz`は[hap, chromosome]のavailabilityと[hap, chromosome, region]の分母を保持します。A/Bは**利用可能な染色体上の合計**であり、未観測染色体を0とみなした全ゲノム合計ではありません。hap間のcoverage差はクラスタリング時に考慮する必要があります。

## 指標

`candidate_metrics.tsv.gz`には以下を保存します（573非参照hap対象）。

- HSat2/HSat3/unionの合算カウント。
- HSat2 count fraction = HSat2 / (HSat2 + HSat3)。分母0はNaN。
- 各familyの合算density、およびlog2((HSat2 density + 1e-9)/(HSat3 density + 1e-9))。この補正値はper-millionにする前の密度に適用します。
- family別・unionの出現hap数。unionについてcount ≥10、≥100のhap数も出力。
- hap別union countとdensityの10/50/90パーセンタイル、densityのCV。
- 両familyに有効開始位置があり当該候補が少なくとも一方に存在するhapをinformativeとして、どちらのfamily densityが高いかを集計。密度が同値の場合はいずれにも加えません。

`chromosome_specificity.tsv.gz`はHSat2/HSat3/unionごとに:

- 合算カウントが最大の染色体と、そのカウント割合。
- 合算count / 各染色体の合算有効開始位置数が最大の染色体と、全染色体の密度の和に対する割合。
- 染色体別カウント分布のShannon entropy（bits、正規化なし）。
- 合算で最多の染色体が、個別hapでも一意に最多となる割合。

最後のrecurrenceの分母は、候補が当該family/unionに存在し、合算最多染色体に当該領域の有効開始位置があり、かつ少なくとも2染色体に有効開始位置があるhapです。合算最多染色体に同率首位があればrecurrenceはNAにします。個別hap内の同率首位は支持数に加えません。最大染色体名の同率首位は`*_tied`列で区別し、表示名は染色体順で最初のものです。カウント皆無はNAです。

`matrices/pooled_chromosome_counts.npz`には、非参照hapを合算した染色体別count/分母/densityを保存します。可用性が異なるhapを含むので、count割合・density割合・recurrence・informative hap数を併せて解釈してください。90%等の集計は記述的指標で、family/染色体特異的ASOの確定ラベルではありません。

染色体特異性は既存HSat2/3注釈の内部での分布です。注釈外にある同一配列や他familyへの結合まで含めた全ゲノムでの染色体特異性は、この指標だけでは判定できません。

## 図

- `plots/01_family_and_haplotype_distribution.png`: family density、familyカウント割合、hap出現率、染色体への集中度。
- `plots/02_chromosome_density_heatmap.png`: 合算E > 10の候補から各family密度上位15件を選び、unionの染色体密度を行の最大値で規格化。例示用で、クラスタリングはしていません。
- `plots/03_specificity_and_recurrence.png`: hap出現率と染色体集中度、最多染色体のhap間再現率。

ヒートマップの灰色は当該染色体の合算有効開始位置数が0で密度を計算できない場合です。白（密度0）と区別します。

図の散布図はseed=42で最大40,000候補を抽出します。集計表は全候補を使用します。

## 再実行・参照

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step2.py
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step2_pipeline.py --workers 6
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/finish_step2.py
```

`finish_step2.py`は全hapカウント終了を待ち、A/B行列・指標・PNG・最終検証を順に実行します。`COMPLETE.json`は最終検証PASS後に書き出します。

任意の候補1配列について、可読TSVを出力できます（候補に含まれない配列はエラーにし、0回とは扱いません）。

```bash
/home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/export_step2.py --kmer AAAAAAAACGGAATTA --hap HG00096_hap1 > /tmp/step2_example.tsv
```

`--hap`を省略すると全574セットを出力します。背景領域のカウントはSTEP 1の結果を継承し、STEP 2では背景の染色体別再カウントを行っていません。1–2 mismatch、hapクラスタリング、KD効果の評価は後続STEPです。
