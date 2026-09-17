# 72 k-mer: CHM13 CenSat・hap group・chromosome landscape

STEP 4で0/1/2 mismatch評価した全72配列について、CHM13のexactヒット座標、CenSat annotation内外、各hapのtarget family（HSat2またはHSat3）内exact count、regional group別・染色体別のcountを作成した。

## 最初に見るファイル

- `PLOT_INDEX.md`: 72枚の個別PNGにリンクする一覧。CenSat外exact hitの多い順。
- `plots/all72_censat_outside_overview.png`: 72配列のCHM13 CenSat外ヒットを一度に比較。
- `chm13_candidate_censat_summary.tsv`: 配列、用途、robustness tier、CHM13のCenSat内外・family別countを一行ずつ記載。
- `plots/kmer_id_landscape.png`: k-merごとの3段図。
- `plots/kmer_id_dominant_zoom.png`: CHM13で最もヒットが多い染色体について、100 kbヒット曲線とCenSat/family annotationを別trackで拡大。CHM13にexactヒットがない2配列はその旨を表示。

## 個別図の読み方

1. 上段: CHM13の全24 filtered chromosome。広いCenSat境界、細かいalphaSat・HSat2・HSat3区間、exactヒットを別trackとして表示する。ヒットは100 kb binに集約し、黒は広いCenSat内、赤は外。右端に各染色体の内外件数を記す。
2. 中段: 対応するtarget familyについて、STEP 3のregional common-core group別、hapあたりfamily BED内exact countの分布。HSat2はG1–G3、HSat3はG1–G7。`NA`はcore領域が不完全でgroup assignmentがないhapで、独立した生物学的groupではない。
3. 下段: group × chromosome別の、measurable hapあたりtarget-family BED内exact countの中央値。数字はraw count、色は`log10(1 + count)`。空白はそのgroupのその染色体で測定可能なfamily regionがない場合。

## 数値出力

- `chm13_all72_exact_hits_with_censat.tsv.gz`: CHM13の全559,297 exactヒット座標とannotation判定。`start_0based`/`end_0based`はBEDと同じ0-based半開区間。canonicalとreverse complementの両方を走査。
- `chm13_censat_by_chromosome.tsv`: 72 k-mer × 24染色体、広いCenSat内外と細かいfamily別count。
- `chm13_censat_totals.tsv`: k-mer別CHM13総計。
- `chm13_100kb_bins.tsv`: 個別図上段のヒット数。
- `hap_chromosome_exact_counts.tsv.gz`: 72 k-mer × 574 sequence set × 24染色体 = 991,872行。target-family BED内exact count、valid 16-mer start数、density、染色体availability、HSat2/3 unionおよびその外側のexact countを含む。
- `hap_total_exact_counts.tsv.gz`: k-mer × hapごとのtarget-family exact総数。
- `group_total_summary.tsv`: regional group別のhap数、target-family countの平均・中央値・四分位範囲。
- `group_chromosome_summary.tsv`: regional group × 染色体のmeasurable hap数、median count、density、10ヒット以上の割合。
- `COMPLETE.json`: 入力と出力の監査記録。

## CHM13 CenSat外exactヒット

72配列中21配列が少なくとも1件、広いCenSat境界外にexactヒットした。72配列の配列別ヒット件数を合計するとCenSat外は55件。最大はHSat2候補`k16_00c14d34`の16件である。これは**配列別のexact 16-mer出現数**であり、55個の異なるゲノム座標や、573 hapにおけるCenSat外ヒット数を意味しない。CHM13で0ヒットの候補も2配列ある。

細かいannotation BEDでは短いalphaSat区間がHSat2などに重なる場合がある。排他的なfine family分類はHSat3 → HSat2 → alphaSat → ctの優先順で完全包含判定し、広いCenSat内外は独立に判定した。CHM13のfine HSat2・HSat3区間はSTEP 0のBEDとそれぞれ完全一致し、全72配列×24染色体の直接座標countはSTEP 4のradius 0行列と一致した。先に解析した上位3配列の座標も新しい走査と一致した。

## 対象範囲と解釈

CHM13位置表はfiltered CHM13 FASTAの24染色体の**全ゲノムexactヒット**を直接走査する。hap側の`target_family_exact_hits`は既存のHSat2またはHSat3 BED区間に完全に入るヒット数で、hap全ゲノムのexact countではない。染色体欠損やfamily region欠損を生物学的0とみなさず、availabilityとvalid start数を分けて保存した。

広いCenSat境界は`VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed`、細かいfamily annotationは`VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed`を使用した。個別図とCenSat内外countは**0 mismatch（完全一致）だけ**を示す。1/2 mismatchを許したCenSat外位置やRNA transcriptome上のoff-targetはここでは測定していない。グループはSTEP 3のregional patternで、whole-genome haplotype groupと解釈しない。

## 再実行

```bash
g++ -O3 -std=c++17 scripts/scan72_chm13.cpp -lz -o scripts/scan72_chm13
scripts/scan72_chm13 step4/mismatch/queries.tsv step0/sequences/chm13v2.0/source.fasta.gz step4/all72_landscape/raw/chm13_all72_exact_hits.tsv.gz
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/all72_landscape.py
Rscript scripts/plot_all72_landscape.R
Rscript scripts/plot_all72_zoom.R
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/finalize_all72_landscape.py
```
