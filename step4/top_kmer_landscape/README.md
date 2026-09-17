# Top 3 HSat3 k-mer landscape

STEP 4で`A_le2mm`に分類されたpan-HSat3上位3配列を対象に、CHM13上のexactヒット座標と、573 non-reference hapのHSat3 annotation内exact countを可視化する。

対象配列は以下のcanonical target 5′→3′である。実際のASO発注配列はRNAの転写orientationで決める。

| k-mer ID | canonical target |
|---|---|
| `k16_3e476af8` | `ATTGCACTCGGGTTGA` |
| `k16_a0d0562c` | `GGAATCAACCCGAGTA` |
| `k16_4118ba0e` | `CAACACGAGTGGAATG` |

## 図

`plots/kmer_id.png`はk-merごとに3段で表示する。

1. CHM13の24 filtered chromosome上で、1 Mb bin別exact hit位置。色はSTEP 0のannotation分類。点の大きさはそのbinのヒット数に対応する。
2. STEP 3のHSat3 common-core regional sequence-pattern group（G1–G7）別、hapあたりHSat3 exact count分布。`NA`はcore chromosomeが不完全な131 hapであり、独立した生物学的群ではない。
3. group × chromosomeの、measurable hapあたりHSat3 exact-hit数の中央値。空白はそのgroupで当該染色体に測定できるHSat3 regionがない場合。数字はraw countの中央値で、色は`log10(1 + count)`を使う。densityは表にも保存する。

`plots/chr9_100kb_zoom_all_three.png`はCHM13 chr9の48–77 Mbを100 kb binに拡大し、3配列のexactヒット数とHSat3 annotationを重ねた。

ユーザーの指摘に応じて、CenSat annotationとヒットを分離した図を追加した。

- `plots/kmer_id_censat_overlay.png`: 24染色体ごとに、広いCenSat境界、細かいalphaSat、HSat2、HSat3、exactヒットを独立したtrackで示す。右側に各染色体のCenSat内外ヒット数を直接表示する。
- `plots/chr9_censat_and_hits_separate_tracks.png`: 3配列の100 kbヒット曲線と、広いCenSat境界・細かいHSat3 annotationを別々のtrackで比較する。

## 表

- `chm13_all_exact_hits.tsv.gz`: CHM13の全exactヒット。各行は1座標・1orientation。`start_0based`/`end_0based`はBEDと同じ0-based半開区間。
- `chm13_chromosome_summary.tsv`: CHM13の各染色体・annotation分類別のexactヒット数。
- `chm13_1Mb_bins.tsv`: 上段図に使った1 Mb bin別ヒット数。
- `hap_chromosome_exact_counts.tsv.gz`: 3 k-mer × 574 sequence set × 24 chromosome。HSat3 exact count、HSat3 valid start数、density、染色体のavailabilityを含む。
- `hap_total_exact_counts.tsv`: hapごとのHSat3 exact総数とdensity。
- `group_total_summary.tsv`: group別のhap数、HSat3 countとdensityの平均・中央値・四分位範囲。
- `group_chromosome_summary.tsv`: group × chromosome別のmeasurable hap数、raw count、density、10ヒット以上の保有率。
- `summary.json`: 検証結果。
- `chm13_all_exact_hits_with_censat.tsv.gz`: 各CHM13ヒットに広いCenSat境界内外と細かいCenSat family分類を追加した座標表。
- `chm13_censat_inside_outside_by_chromosome.tsv`: 3配列 × 24染色体のCenSat内外・family別ヒット数。
- `chm13_censat_inside_outside_totals.tsv`: k-mer別の全染色体合計。
- `censat_overlay_summary.json`: BED source SHA256と区間一致・count恒等式の検証。

## 定義と範囲

CHM13位置は、STEP 0のfiltered CHM13 FASTAの24染色体を直接走査し、canonical targetとreverse complementを両方数えた。位置図は全filtered染色体を扱う。hap側の染色体countは、STEP 4のradius 0（exact match）countから取り出した。`HSat3_exact_hits`はHSat3 BED内に完全に収まる16-merだけであり、ゲノム全域countとは異なる。CHM13のHSat3直接走査countとSTEP 4のmatrix countは、3配列 × 24染色体で一致した。

グループ別summaryと図はHSat3 annotation内exact countに基づく。染色体欠損はavailability maskで管理し、生物学的な0 countとして扱わない。群はSTEP 3のregional patternで、whole-genome haplotype groupと解釈しない。

## CenSat内外の直接確認

広いCenSat境界は`VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed`の24染色体区間、細かいfamily annotationは`VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed`を使用した。細かいBEDでは短いalphaSat区間がHSat2などと重なるため、各ヒットの排他的なfine categoryはHSat3 → HSat2 → alphaSat → ctの優先順で完全包含判定した。広いCenSat内外はfine familyとは独立に判定した。細かいBEDのHSat3 247区間はSTEP0のHSat3 BEDと完全一致した。

3配列のCHM13全exactヒット62,603件はすべて広いCenSat境界内、かつ細かいHSat3 annotation内に完全包含された。CHM13でのCenSat外exactヒットは各配列0件である。この0件は**exact 16-merに限る**。1または2 mismatchを許した場合のCenSat外ヒット数や、RNA transcriptome上のoff-targetを0と示すものではない。

## 再実行

```bash
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/top_kmer_landscape.py
Rscript scripts/plot_top_kmer_landscape.R
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/annotate_top_kmer_hits.py
Rscript scripts/plot_top_kmer_censat_overlay.R
```
