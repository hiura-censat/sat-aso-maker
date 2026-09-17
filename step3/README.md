# STEP 3: regional haplotype patterns and mismatch evaluation

STEP 2の159,988候補を出発点に、HSat2・HSat3別のクラスタリング、候補のグループ比較、24代表配列の0–2 mismatch評価を行います。元のSTEP 1/2は変更していません。`COMPLETE.json`は全工程と最終検証が完了した場合だけ作成します。

## 結果の解釈

主解析のグループは**共通の注釈領域における配列パターン群**です。全ゲノムを通じた固定的なhapグループとは解釈しません。

- HSat2: chr7 + chr10が測定可能な489非参照hap、3群（158/42/289 hap）。主にchr10の分割と一致します（ARI約0.888）。
- HSat3: chr5 + chr7 + chr10 + chr20が測定可能な442非参照hap、7群（59/54/70/39/47/47/126 hap）。
- 全利用可能領域での分割とのARIはHSat2約0.0014、HSat3約0.0027で、共通した分割は確認できません。
- PC1と利用可能染色体数との相関は両familyとも小さい一方、共通領域の長さとの相関はHSat2で約0.685、HSat3で約0.327です。領域の長さ・染色体構成と配列分布の関係を含むパターンです。

共通領域が欠けるhapは主解析に割り当てずNAにします。CHM13は学習・合算から除外し、PCAへの投影と最も近い群を別に保存します。人口集団・祖先集団等のラベルは推定していません。

## 入力と前処理

- `../step2/candidates.tsv`, `candidate_codes.u32`: 159,988 canonical 16-mer。
- `../step2/counts/{hap}.npz`: family × chromosome × orientationの正確なカウント。
- `../step2/matrices/matrix_B_{counts,density_per_M}.npy`: 全利用可能領域の比較用。
- `../step2/matrices/availability_and_denominators.npz`: T2T可用性と領域分母。
- ミスマッチ再カウント: STEP 0の既存BEDとmanifestに記載したFASTA。FASTAのコピーは作りません。

573非参照hapの85%以上で当該familyの有効16-mer開始位置が存在する染色体を共通領域に採用します。そのすべてが揃うhapだけで主解析を行います。測定可能な染色体で候補がない場合の0と、染色体欠測を区別します。

主解析は共通染色体のcount合計 / 有効開始位置数合計 × 1e6を`log1p`変換します。染色体を同じ重みで平均した密度ではありません。各解析でcount >=10が5%以上のhapにあり、変換後分散 >1e-4の特徴量を対象に、分散上位3,000を選びます。分散順に、Pearson相関の絶対値 >=0.95の特徴量を同じ代表にまとめ、最大500代表を採用します。主解析ではHSat2 285代表、HSat3 131代表で上位3,000特徴をカバーしました。

クラスタリングでは特徴量を平均中心化相当で扱い、分散を1にそろえるスケーリングは行いません。図のヒートマップだけはz-score表示です。全159,988配列のカウントは保持し、候補評価には全配列を使用します。

## PCA・クラスタリング・安定性

PCAはNumPyのSVD、階層クラスタリングはRの`hclust(dist(X), method='ward.D2')`を使用します。PCAは可視化用で、クラスタリングの距離は代表特徴量すべてのEuclidean距離です。

k=2–8について平均silhouette、最小群サイズを計算します。主解析では代表特徴量の80%を非復元抽出する20回の感度解析（seed=42）を行い、元の分割とのARIと、各群の最良対応Jaccardの最小値を計算します。

以下をすべて満たす分割の中からsilhouette最大（同値なら小さいk）を採用します。

- 平均silhouette >=0.25。
- 最小群サイズ >=10。
- 特徴量再抽出のARI中央値 >=0.75。
- 最も再現しにくい群の対応Jaccard中央値 >=0.75。

これらは今回の探索的な採用基準で、普遍的な最適値ではありません。合格する分割がなければ`unresolved`とします。20回の感度解析は**特徴量選択に対する再現性**であり、独立個体を再サンプリングした信頼区間や独立検証ではありません。

HSat2主解析: silhouette約0.596、ARI中央値1.0。HSat3主解析: silhouette約0.299、ARI中央値約0.970。安定性が高くても、他染色体を含めた同じ分割を保証しません。

補助解析は全利用可能領域（特徴量再抽出5回）と、100 hap以上で測定可能な染色体ごとの解析です。染色体別解析では全領域データで選んだfamily別分散上位5,000候補の中から特徴量を選びます。したがって、少数染色体に固有の低頻度特徴を網羅する解析ではありません。十分な変動特徴がない場合は分類不能とし、変動が生物学的に存在しないとは結論しません。

## 完全一致候補の分類

`candidates/`に全候補評価の結果を保存します。

- `HSat2_pan_exact.tsv`, `HSat3_pan_exact.tsv`: 合算family count割合 >=99.9%、他family対比E >10、背景対比E >10、全573 hapの90%以上でfamily count >=10。さらに各主解析群の90%以上で**共通領域内**のcount >=10を要求した保守的な集合です（HSat2 2,470、HSat3 5,254配列）。
- `*_G*_markers.tsv`: 共通領域で群内count >=10の割合 >=80%、群内/群外の密度中央値比 >=4倍、合算family count割合 >=99%、他family対比E >6、背景対比E >10。群外のcount >=10の割合 <=20%なら`group_specific_exact`、それ以外は`group_enriched_exact`と記録します。
- `*_chromosome_exact.tsv`: 合算countの90%以上が1染色体に集中し、informative hap >=100、最多染色体の再現率 >=80%、family count割合 >=99.9%、他family・背景対比ともE >10。

group候補の密度中央値比にはper-million密度に0.001を加えます（元の密度の1e-9に相当）。グループ作成と候補探索に同じデータを使うため、効果量による探索的候補であり、独立検証済みのマーカーではありません。p値・FDRは算出していません。グループ候補は共通領域での特異性で、全染色体にわたるKD特異性とは区別します。

`pan_exact`はカウント出現率による予備選抜であり、前回のHSat2厳密候補1,797配列とは条件が異なります。各hapで背景・他family対比Eも満たすかは、代表配列のミスマッチ評価で再判定します。

## 24配列の0–2 mismatch評価

`candidates/mismatch_queries.tsv`が今回評価する代表配列です。pan・群別・染色体候補を順番に選び、最大24配列に制限します。各選抜区分内では、既選択配列またはそのRCとHamming距離 <=2、または共通12-merを持つ配列を避けます。これは簡単な冗長性抑制で、全候補の完全な類似配列クラスタリングではありません。

24代表は全1797/159988候補の代わりではありません。未評価の候補にミスマッチ結果を外挿しません。以前提示した`CATCAAACGGAATCAA`は今回の対象に含まれます。

各16-merの0・1・2塩基置換をすべて列挙します。1配列あたり最大1 + 16×3 + C(16,2)×9 = 1,129近傍で、24配列全体のcanonical近傍は27,054種類、query-variant関係は27,096件です。全574セットについてHSat2、HSat3、union、背景を再走査します。

- Hamming距離のみ。挿入・欠失は対象外。
- 各ゲノム16-merとquery/RC-queryの距離の最小値を使用。
- 1つのゲノム開始位置を、同じqueryに対して一度だけ数える。
- 距離0、ちょうど1、ちょうど2を別に保存し、<=1、<=2の累積値も評価。
- 領域境界をまたぐ16-merと非ACGTを含む窓は除外。
- ミスマッチ位置はcanonical queryの5′端から1-based。両方向が同距離で位置が異なる場合、位置集計には等分して加算。全代替位置をvariant catalogに記録。

カウンタのcanonical prefix bitmapは近傍集合にない配列を早く除外するための正確なフィルタです。許容近傍を捨てる近似はしません。途中で標準入力の読み込みを高速化しましたが、HG00096_hap1の全variantカウント・queryスペクトル・位置集計が最適化前後で一致することを確認しています（`mismatch/io_validation/validation.json`）。

### ミスマッチ出力

- `mismatch/pooled_mismatch_scores.tsv`: 許容距離0/1/2ごとのfamily・他family・背景カウント、E、family割合、厳密条件を満たすhap数。
- `mismatch/candidate_robustness.tsv`: strict pan条件を各距離で維持するか。合算でfamily割合 >=99.9%、他family/背景E >10、かつ516/573 hap以上で個別にも同条件とtarget count >=10を満たす必要があります。
- `mismatch/group_contrast_by_mismatch.tsv`: 主解析と同じ共通染色体での群内外比較。
- `mismatch/group_contrast_all_available_by_mismatch.tsv`: 利用可能な全染色体上での群内外比較。染色体可用性の差が残ることに注意。
- `mismatch/selected_group_robustness.tsv`: 選択元のspecific/enriched条件とfamily/background条件を各距離で維持するか。
- `mismatch/chromosome_robustness.tsv`: 最多染色体・集中度・hap間再現率が各距離でどう変わるか。
- `mismatch/mismatch_position_counts.tsv`: 非参照hap合算の位置別カウント。
- `mismatch/hap_counts_and_positions.npz`: query × hap × region × exact distanceのcount、位置情報、共通領域カウント、染色体別カウントと可用性・分母。
- `mismatch/counts/{hap}.npz`: query × chromosome × region × exact distanceの各hap結果。
- `mismatch/variants/{hap}.npz`: variant × chromosome × regionのcanonicalカウント。
- `mismatch/variant_catalog.tsv`: variantとquery、最小距離、ミスマッチ位置の対応。

region順はHSat2, HSat3, HSat23, background、距離軸は0,1,2です。NPZの可用性マスクを参照し、未測定染色体を生物学的な0と解釈しないでください。

これらは既存HSat2/3注釈とその外側に対する**DNA配列一致**の評価です。RNAの発現・転写鎖・ASO化学修飾・結合強度・KD効率を表すスコアではありません。αSatは独立familyとして注釈していません。ASOの最終方向と発現転写産物へのoff-target評価にはRNA側の情報が必要です。

## ファイルと検証

- `hap_groups.tsv`: 全574セットの主解析群、欠測、参照状態。
- `clustering/{case}/`: PCA、全kの分割、silhouette/安定性、代表特徴、R階層木、モデル。
- `cluster_concordance.tsv`, `coverage_associations.tsv`: 領域間の一致度とcoverage/領域長の関係。
- `plots/`: coverage、HSat2/HSat3のPCA・階層木・ヒートマップ、領域間一致度、ミスマッチ特異性のPNG計5枚。
- `validation_summary.json`, `COMPLETE.json`: 完了時の検証結果と実行スクリプトハッシュ。

カウンタは、小文字・N・境界・回文・RC・高カウントを含むfixtureでPython総当たりと比較します。全574セットで距離0がSTEP 1の両方向カウントと一致し、有効開始位置数がSTEP 0と一致することを検証します。位置集計の和が距離×一致数に等しいことも確認します。最終検証では近傍集合の網羅性・重複排除、代表特徴の再構築、欠測・CHM13の扱い、12 hapの保存variantからのスペクトル再構築を確認します。

8 hapの集約試験では任意の仮ラベルを使用して出力形状・計算を確認しました。その仮ラベルや試験結果は科学的なグループ解析には使用せず、一時ディレクトリ内だけで実行しています。

再実行（プロジェクト直下）:

```bash
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step3_prepare.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step3.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/run_step3.py
```

入力や手法を変更して再解析する場合は新しい出力先を使うか、対応するcheckpointを明示的に更新してください。

方法の参照: [R hclust / Ward.D2](https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html)、[R cluster silhouette](https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/silhouette.html)。
