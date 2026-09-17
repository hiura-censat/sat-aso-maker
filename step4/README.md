# STEP 4: ASO candidate classification and prioritization

## 目的

STEP 1–3の結果を統合し、HSat2/HSat3を標的とする16-merを、次の用途別に優先順位付けする。

- `pan`: 多数のhapで同じfamilyを広く標的化する候補
- `regional_group`: STEP 3の共通core領域で定義したregional sequence-pattern group候補
- `chromosome`: 特定染色体のannotated target領域に集中する候補

## 入力

- `step2/candidate_metrics.tsv.gz`: family count、family specificity、hap prevalence
- `step1/pooled_counts_and_scores.npz`: canonical backgroundを含むpooled exact count
- `step3/candidates/*.tsv`: pan、regional group、chromosome候補
- `step3/mismatch/*`: STEP 3で評価済みの24候補
- `step3/clustering/*`: regional group assignmentとcore chromosome定義
- `step0/manifest.tsv`から参照される574 sequence setと既存HSat2/3 BED由来領域

## 処理

### 1. Exact-match候補の統合

重複をまとめた31,335個の候補をHSat2/HSat3別に評価した。`exact_composite_score`は各family内の経験percentileを次の重みで統合する。

| 成分 | 重み |
|---|---:|
| target abundance | 0.20 |
| HSat2/HSat3 family specificity | 0.20 |
| target family vs background enrichment | 0.20 |
| hap coverage | 0.15 |
| absolute backgroundの少なさ | 0.10 |
| 配列品質proxy | 0.15 |

配列品質proxyはGC%、最大homopolymer、base entropy、自己相補stem proxyだけを使う。Tm、化学修飾、RNA構造、転写量、細胞内取り込み、毒性は含まない。

### 2. 多様な代表候補の選択

用途別に高score候補を選び、同じcategory内でHamming距離2以下、または12-merを共有する配列を除いて冗長性を下げた。STEP 3の24候補を再利用し、新規48候補を加えた計72候補をmismatch評価対象とした。

### 3. 0/1/2 mismatch評価

16-merの両orientationについて、Hamming radius 0、1、2の全置換近傍を重複排除して数えた。indelは扱わない。574 sequence set、eligible T2T chromosome、既存HSat2/3 annotationとそのbackground complementが評価範囲である。

panのstringent判定は、pooledで以下をすべて満たし、かつhap内条件を573 non-reference hap中516以上（90%）で満たす場合とした。

- target vs background `E > 10`
- target vs other family `E > 10`
- target family count fraction `>= 0.999`
- hap内target count `>= 10`

Regional groupはcore領域におけるgroup内presenceとdensity ratio、chromosome候補はannotated target領域内のcount concentrationとhap間recurrenceを評価する。

### 4. robustness tier

- `A_le2mm`: 2 mismatchを含めても用途別screenを通過
- `B_le1mm`: 1 mismatchまで通過
- `C_exact`: exact matchだけ通過
- `D_fails_exact`: 代表候補のexact用途別screenを通過しない

`D`は入力段階のexact候補選定が誤りという意味ではない。STEP 4の目的別screenがより厳しい、または同じ配列が複数categoryに属し別categoryから代表選定された場合を含む。

## 主な出力

- `ranked/all_exact_ranked.tsv.gz`: 31,335 exact候補の統合ranking
- `shortlists/exact_diverse_shortlist.tsv`: mismatch展開前の用途別・非冗長shortlist
- `ranked/mismatch_evaluated_ranked.tsv`: 72候補の最終ranking、robustness tier、ASO orientation欄
- `shortlists/HSat2_pan_candidates.tsv`, `HSat3_pan_candidates.tsv`
- `shortlists/HSat2_regional_group_candidates.tsv`, `HSat3_regional_group_candidates.tsv`
- `shortlists/HSat2_chromosome_candidates.tsv`, `HSat3_chromosome_candidates.tsv`
- `mismatch/pooled_mismatch_scores.tsv`: radius別pooled specificity
- `mismatch/candidate_robustness.tsv`: stringent pan判定
- `mismatch/group_robustness.tsv`: regional group contrast
- `mismatch/chromosome_robustness.tsv`: chromosome concentration
- `mismatch/mismatch_position_counts.tsv`: mismatch位置別count
- `plots/*.png`: ranking、mismatch specificity、tier分布
- `COMPLETE.json`: 完了検証記録

## 配列orientationの読み方

`canonical_target_5to3`はcanonical化したgenomic target配列で、転写鎖を意味しない。実際のASO配列は標的RNAに相補的である必要があるため、最終表には両方の可能性を明記した。

- canonical targetが転写される場合: `aso_5to3_if_canonical_target_is_transcribed`
- reverse-complement targetが転写される場合: `aso_5to3_if_reverse_complement_target_is_transcribed`

転写orientationと発現を実験系で確認してから発注配列を決める。

## 実行方法

```bash
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_rank.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_neighbors.py
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_mismatch.py --workers 8
OPENBLAS_NUM_THREADS=2 /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step4_finalize.py
```

## 解釈上の範囲

このSTEPはgenome配列上の標的候補rankingである。regional groupはSTEP 3のcommon-core regional patternで、確立したwhole-genome haplotype groupではない。RNAへの到達性、発現量、ASO chemistry、RNase H活性、毒性、transcriptome-wide off-targetはwet-lab選定前に別途評価する。
