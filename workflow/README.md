# Snakemake workflow

このワークフローは既存の STEP 0–4 と top 3 / 72 配列 CHM13 位置解析を、依存関係を持つ 20 ジョブとして実行します。結果は `runs/<run_id>/` に作成します。ルートの既存 `step0/`〜`step4/` は変更しません。実行時には現在の `scripts/` を run directory にコピーし、その SHA256 と設定ファイルの SHA256 を `.bootstrap.json` に記録します。同じ `run_id` で設定やコードが変わった場合、混在を防ぐため停止します。新しい `run_id` を使ってください。

## 設定

[config.example.yaml](config.example.yaml) を `workflow/config.yaml` にコピーし、自分の環境に合わせて編集します。`workflow/config.yaml` はローカル専用で Git 管理しません。以下の項目は計算に反映されます。

| 項目 | 効果 |
| --- | --- |
| `run_id` | 独立した出力先 `runs/<run_id>/` |
| `source_root` | VallePrep データディレクトリ。project root からの相対パスまたは絶対パス |
| `python`, `samtools` | 使用する実行ファイルの絶対パス |
| `threads.*` | STEP 0/1/2 と STEP 3/4 mismatch の同時 worker 数 |
| `selection.selected_E` | STEP 2 に渡す HSat23/background の strict `E > selected_E` |
| `selection.sweep_E` | 比較・出力する strict E 閾値の昇順リスト。`selected_E` を含める |
| `selection.pooled_min_count` | 閾値スイープの pooled HSat23 count 下限 |
| `selection.hap_min_count` | 少なくとも 1 hap の HSat23 count 下限 |
| `selection.broad_hap_fraction` | broad 候補フラグの個別合格 hap 割合 |

初期値は元の解析に対応する `E > 10`、pooled count ≥100、hap count ≥10、broad 90% です。STEP 1 の全 k-mer 発見と元の `E >= 2` 採用表は保存し、**STEP 2 以降に渡す候補集合**を上記 selection で変更します。pooled/hap count 下限を 100/10 より小さくできないのは、STEP 1 の preliminary universe がこの下限で作られているためです。下限を小さくするには発見段階から作り直す実装変更が必要です。

`invariants.k=16` と `invariants.max_mismatches=2` は現在の C++ rolling counter、2-bit ID、0/1/2 mismatch 集計に埋め込まれた方法上の固定値です。別の値は起動時にエラーにします。STEP 3 の group 選択基準、STEP 4 の ASO screen と score の重みは現状の科学的手法として固定です。上記 `broad_hap_fraction` はスイープ表と STEP 2 のフラグに作用し、STEP 3/4 の pan screen が要求する 90% を変更するものではありません。候補選択を極端に厳しくすると下流の pan/top 3 候補が不足し、明示的なエラーで停止する可能性があります。

## 実行と確認

プロジェクト直下で、設定した Python 環境の Snakemake を使って実行します。`snakemake` が PATH にある例です。PATH にない場合はその実行ファイルの絶対パスに置き換えます。

```bash
cp workflow/config.example.yaml workflow/config.yaml
snakemake --snakefile Snakefile --cores 8 --dry-run
snakemake --snakefile Snakefile --cores 8
```

`--dry-run --summary` で不足・更新予定の出力を一覧できます。最終成果物は `runs/<run_id>/step4/all72_landscape/COMPLETE.json` と同ディレクトリの `PLOT_INDEX.md` です。途中までの実行例は `.../snakemake --snakefile Snakefile --cores 8 runs/<run_id>/step2/COMPLETE.json` ではなく、ステージの marker `runs/<run_id>/.workflow/step2_finish.done.json` を target にしてください。各 marker はその段階のスクリプトが正常終了し、主要出力がある場合だけ書かれます。最終段階は top 3 座標の独立照合にも依存します。

同じ `run_id`・設定・コード・入力を維持した再実行では、完了済みルールを再実行しません。処理の途中で失敗した場合は原因を直して同じコマンドを再実行できます。ただし元 VallePrep FASTA/BED の内容を変更した場合は、新しい run ID で作り直してください。元データは symlink 経由で参照しており、ワークフローは元データを複製しません。

## 段階と検証

`step0` は `verify_step0.py`、`step1_finish` は `verify_step1.py`、`step2_finish` は `verify_step2.py`、`step3_run` と `step4_finish` は各 finalize、`landscape_finish` は全配列図と CHM13 位置の finalize を実行します。スクリプト内の既存 checkpoint に加え、STEP 1/2/4 の hap 別 shard の存在を stage driver で確認します。重要な入出力の内容照合は各 STEP の既存検証スクリプトが行います。

このワークフローは**段階単位**です。STEP 1/2 と mismatch スクリプトが内部で hap を並列処理するため、Snakemake 側では一段階を一ジョブとして扱い、`threads` を内部 worker 数に渡します。hap 単位ジョブへ分割するには既存スクリプトに hap 指定 CLI を追加する必要があります。
