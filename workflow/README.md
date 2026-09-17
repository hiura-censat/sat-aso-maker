# Snakemake workflow

このワークフローは既存の STEP 0–4 と top 3 / 72 配列 CHM13 位置解析を、依存関係を持つ 19 ジョブとして実行します。閾値 sweep ジョブはありません。結果は `runs/<run_id>/` に作成します。ルートの既存 `step0/`〜`step4/` は変更しません。実行時には現在の `scripts/` を run directory にコピーし、その SHA256 と設定ファイルの SHA256 を `.bootstrap.json` に記録します。同じ `run_id` で設定やコードが変わった場合、混在を防ぐため停止します。新しい `run_id` を使ってください。

## 設定

[config.example.yaml](config.example.yaml) を `workflow/config.yaml` にコピーし、自分の環境に合わせて編集します。`workflow/config.yaml` はローカル専用で Git 管理しません。以下の項目は計算に反映されます。

実行前に Python 環境へ Snakemake、NumPy、PyYAML、samtools を用意し、`g++`、zlib 開発ライブラリ、R と `cluster` パッケージを使えるようにしてください。Jellyfish 2.3.1 の実行ファイルを `.local/bin/jellyfish` に配置し、`scripts/jellyfish-prefix12.patch` を適用した版を `.local/bin/jellyfish-prefix` に配置します。C++ scanner は最初の bootstrap ジョブがソースからビルドします。元 FASTA/BED は `source_root` に配置し、GitHub には含めません。

| 項目 | 効果 |
| --- | --- |
| `run_id` | 独立した出力先 `runs/<run_id>/` |
| `source_root` | VallePrep データディレクトリ。project root からの相対パスまたは絶対パス |
| `python`, `samtools` | 使用する実行ファイルの絶対パス |
| `threads.*` | STEP 0/1/2 と STEP 3/4 mismatch の同時 worker 数 |
| `selection.selected_E` | STEP 1 本選抜の HSat23/background の strict `E > selected_E` |
| `selection.pooled_min_count` | pooled HSat23 count 下限（予備候補・本選抜） |
| `selection.hap_min_count` | 少なくとも1 hapの HSat23 count 下限（予備候補・本選抜） |
| `selection.broad_hap_fraction` | broad 候補フラグの個別合格 hap 割合 |

新しい初期値は `E > 10`、pooled count ≥10、少なくとも1 hapの count ≥500、broad 90% です。STEP 1 の予備候補は **pooled count ≥10 OR 1 hap count ≥500**、本選抜は **(pooled count ≥10 AND pooled E >10) OR (1 hap count ≥500 ANDそのhapの E >10)** です。E は `log2((HSat23 density + 1e-9)/(background density + 1e-9))` です。CHM13 は pooled/個別 hap 判定に含めません。STEP 2 は `step1/candidate_kmers.tsv` を直接読みます。

同じ573 hapを合算するため、予備候補の count 条件では「1 hap ≥500」は必ず「pooled ≥10」に含まれます。一方、本選抜では pooled E と個別 hap E が異なるため、二つの経路を残します。新設定の候補数は再計算まで未確定です。ルートの `step0/`〜`step4/` にある数値は旧設定の結果です。

STEP 3 の主解析は **HSat2/HSat3 × 染色体ごとに独立したクラスタリング**です。各染色体で注釈領域の有効窓がある非参照 hap だけを用い、測定可能 hap が100未満ならその組み合わせをスキップします。特徴量候補も染色体別 pooled count から最大5,000配列を選びます（高カウント側の半数と残りのカウント分布から均等抽出した半数）。その後の変動特徴の選択、正規化、PCA、階層クラスタリング、群判定も染色体ごとに行います。`step3/hap_chromosome_groups.tsv` に全 hap × family × 染色体の測定状態、群、未割当理由を記録します。同じ `G1` でも染色体が異なれば別群です。群特異的候補と mismatch 判定は、同一染色体で測定可能な hap の群内外だけを比較します。pan と染色体特異性の全染色体評価は維持します。これは HOROSCOPE の**染色体別解析単位**の採用であり、同論文の61-merや類似度式の再実装ではありません。

STEP 4 の top 3 図は新しい ranking から pan-HSat3 配列を3件選びます。mismatch評価は最大72配列で、候補数や多様性条件により72未満の場合も、その実数で後続処理と検証を行います。群別集計は染色体内の群だけを表示します。旧解析で選んだ固定の3配列を強制しません。

`invariants.k=16` と `invariants.max_mismatches=2` は現在の C++ rolling counter、2-bit ID、0/1/2 mismatch 集計に埋め込まれた方法上の固定値です。別の値は起動時にエラーにします。STEP 3 の group 選択基準、STEP 4 の ASO score の重みは現状の科学的手法として固定です。`broad_hap_fraction` はSTEP 2だけでなくSTEP 3/4のpan screenにも反映され、必要hap数は `ceil(non-reference hap数 × broad_hap_fraction)` で計算します。候補選択を極端に厳しくすると下流のpan/top 3候補が不足し、明示的なエラーで停止する可能性があります。

STEP 3のcluster作成とgroup marker探索は、各family × chromosomeで選んだ最大5,000特徴を対象とします。これは計算量を制限する探索的解析であり、全STEP 2候補をcluster確定後に再スコアする独立検証ではありません。clusterとmarkerは同じhap cohortから推定されるため、群特異性は独立検証済み性能ではなく探索的な効果量として扱います。

## 実行と確認

プロジェクト直下で、設定した Python 環境の Snakemake を使って実行します。`snakemake` が PATH にある例です。PATH にない場合はその実行ファイルの絶対パスに置き換えます。

```bash
cp -n workflow/config.example.yaml workflow/config.yaml
snakemake --snakefile Snakefile --cores 8 --dry-run
snakemake --snakefile Snakefile --cores 8
```

`--dry-run --summary` で不足・更新予定の出力を一覧できます。最終成果物は `runs/<run_id>/step4/all72_landscape/COMPLETE.json` と同ディレクトリの `PLOT_INDEX.md` です。途中までの実行例は `.../snakemake --snakefile Snakefile --cores 8 runs/<run_id>/step2/COMPLETE.json` ではなく、ステージの marker `runs/<run_id>/.workflow/step2_finish.done.json` を target にしてください。各 marker はその段階のスクリプトが正常終了し、主要出力がある場合だけ書かれます。最終段階は top 3 座標の独立照合にも依存します。

同じ `run_id`・設定・コード・入力を維持した再実行では、完了済みルールを再実行しません。処理の途中で失敗した場合は原因を直して同じコマンドを再実行できます。ただし元 VallePrep FASTA/BED の内容を変更した場合は、新しい run ID で作り直してください。元データは symlink 経由で参照しており、ワークフローは元データを複製しません。

## 段階と検証

`step0` は `verify_step0.py`、`step1_finish` は `verify_step1.py`、`step2_finish` は `verify_step2.py`、`step3_run` と `step4_finish` は各 finalize、`landscape_finish` は全配列図と CHM13 位置の finalize を実行します。スクリプト内の既存 checkpoint に加え、STEP 1/2/4 の hap 別 shard の存在を stage driver で確認します。重要な入出力の内容照合は各 STEP の既存検証スクリプトが行います。

このワークフローは**段階単位**です。STEP 1/2 と mismatch スクリプトが内部で hap を並列処理するため、Snakemake 側では一段階を一ジョブとして扱い、`threads` を内部 worker 数に渡します。hap 単位ジョブへ分割するには既存スクリプトに hap 指定 CLI を追加する必要があります。
