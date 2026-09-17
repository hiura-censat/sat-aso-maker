# STEP 1: exact canonical 16-mer discovery

> このディレクトリ内の数値と以下の旧手順は旧設定（pooled ≥100 または 1 hap ≥10、初回 E ≥2）の保存済み結果です。現在の Snakemake ワークフローは pooled ≥10 または 1 hap ≥500 の予備候補から、各経路で `E > 10` を直接評価し、sweep を実行しません。新しい結果は別の `runs/<run_id>/` に保存します。詳細は [workflow README](../workflow/README.md)を参照してください。

This run uses Jellyfish 2.3.1 (installed under `.local/`) and the validated STEP 0
inputs. CHM13 is evaluated separately and excluded from pooled candidate
selection. Only the 573 eligible non-reference haplotypes contribute to pooled
counts, pooled denominators and hap prevalence.

## Algorithm

1. Extract each HSat23 union interval as a separate FASTA record using samtools
   faidx and the existing FASTA/FAI/GZI files. BED [start,end) becomes faidx
   start+1 through end. Run Jellyfish `-m 16 -C` without a count cutoff other than
   excluding zero entries. Check the total against STEP 0's valid 16-mer starts.
2. Retain a preliminary canonical candidate if its non-reference pooled target
   count is >=100 OR any non-reference hap target count is >=10. These are
   necessary conditions for the final selection rule; discarded rare target
   k-mers cannot pass that rule. Full target-discovery counts remain available.
3. For every hap/reference, count HSat2, HSat3 and HSat23 without canonical
   collapse, then select the preliminary candidates. Count background with
   Jellyfish `--if`, seeded with each preliminary candidate and its reverse
   complement. A self reverse-complementary sequence is seeded only once.
   No background-only k-mer database is retained. All temporary Jellyfish
   databases are converted to compressed counts and removed after success.
4. Calculate density = count / valid 16-mer starts, and
   E = log2((target_density + 1e-9)/(background_density + 1e-9)).
   Final selection is (pooled E>=2 AND pooled target count>=100) OR
   (at least one non-reference hap with E>=2 AND target count>=10).
   No GC/complexity/prevalence filter is applied in STEP 1.

The union is counted independently: HSat2 + HSat3 counts may not equal HSat23
counts because a 16-mer can cross a family boundary inside the union. N and
other non-ACGT bases break windows. Separate intervals are never concatenated.
Orientation is relative to the filtered, strand-oriented FASTA and the canonical
representative, not the transcriptional direction or an ASO strand decision.

Background means the complement of the existing HSat2/3 annotations on eligible
T2T-passed chromosomes. It can include other satellites and missed annotations;
the original Altemose search used preselected regions. Missing chromosomes are
not treated as zeros. Hap prevalence is conditional on the available STEP 0
chromosomes and is not yet adjusted for chromosome-specific missingness.

## Counter precision and checks

The locally built Jellyfish 2.3.1 produced zero counts with an 8-byte output
counter in the pilot. This run uses its default 4-byte output counters instead.
Every per-hap region has fewer than 2^32 valid windows (checked before counting),
so no possible per-k-mer count can overflow. All conversions, hap aggregates and
pooled sums use uint64. Internal Jellyfish counters use `-c 16`, which is not a
count cap; expansion beyond 65535 is checked by a synthetic high-count test.
Synthetic tests also compare Jellyfish against brute-force counts with N,
lowercase, record boundaries, reverse complements, palindromes and seeded
candidate counting. Each target region's full count sum must match STEP 0.
The independently counted canonical union must agree exactly with the sum of
its two orientation counts. Background counts are candidate-restricted and
therefore do not sum to the full background window denominator.

## Files and formats

- `config.json`: exact thresholds, software and region/orientation definitions.
- `discovery/<hap>.bin.gz`: all canonical HSat23 counts, repeated little-endian
  uint32 encoded k-mer + uint64 count (12 bytes per record).
- `preliminary_candidates.u32`: sorted preliminary canonical codes, little-endian
  uint32. This is the shared row index for every count/enrichment NPZ file.
- `counts/<hap>.npz`: `counts`, shape [preliminary_candidates, 4 regions, 2
  orientations], uint64. Region order: HSat2, HSat3, HSat23, background.
  Orientation 0 is the canonical sequence; 1 is its reverse complement.
  Palindromic k-mers have their count in orientation 0 and zero in orientation 1.
- `counts/<hap>.json`: QC results and source/candidate fingerprints.
- `enrichment_by_hap/<hap>.npz`: target_density, background_density,
  log2_enrichment and the four valid-start denominators, in shared row order.
- `pooled_counts_and_scores.npz`: non-reference aggregate and selection masks.
- `kmer_catalog.tsv.gz`: sequence identities and unfiltered GC/complexity metrics
  for every preliminary candidate. `forward_sequence` is the canonical
  representative; actual orientation-specific occurrences are in the counts.
- `enrichment_summary.tsv.gz`: scores for all preliminary candidates, including
  rejected candidates and separate CHM13 comparison columns.
- `candidate_kmers.tsv`: selected candidates with selection reasons.
- `candidate_codes.u32`: sorted selected codes for downstream exact queries.
- `qc_summary.tsv`: count reconciliation per hap and region.
- `pooled_threshold_sensitivity.tsv`: alternative pooled selection cutoffs,
  restricted to the preliminary universe (not a new unrestricted discovery).
- `logs/`: Jellyfish/samtools diagnostics and `/usr/bin/time -v` resource reports.

A code uses A=0, C=1, G=2, T=3 and 2 bits per base, left to right. The stable ID is
`k16_` followed by the eight-digit hexadecimal code. Reverse complements are
lexicographically canonicalized. Binary arrays do not require Python pickles.

## Reproduction

Run from the project root with the NumPy-enabled Python environment:

    /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/test_step1.py
    /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step1_pipeline.py --phase discovery --workers 8
    /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step1_pipeline.py --phase counts --workers 8
    /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/step1_pipeline.py --phase score

All reference FASTA inputs remain in their original locations. No whole-genome
FASTA copies or genome-wide GenMap indexes are created. Completion is recorded
only after all counting, scoring and final verification succeed.

## Background acceleration

The production background counter is the locally compiled `jellyfish-prefix`
variant. The patch is preserved in `scripts/jellyfish-prefix12.patch`; the stock
binary remains available as `.local/bin/jellyfish`. With
`SAT_ASO_PREFIX_FILTER=1`, for noncanonical 16-mer counting with `--if`, a 2 MB
bitmap rejects 12-base prefixes absent from every supplied query. The original
full-key `--if` lookup still decides every accepted count. The gate cannot reject
a listed query because its own prefix is inserted in the bitmap. Seeds must be
plain FASTA with one uppercase 16-mer per sequence line; other formats fail
explicitly rather than producing partial results. This is an exact prefilter,
not approximate counting. It is disabled for ordinary Jellyfish usage.

Stock and prefix-gated counts were compared on an entire chromosome, and on all
background intervals of HG00096_hap1, including both orientations. Stock pilot
counts for the first three sequence sets are retained. The Bloom-counter variant
was tested separately but was slower and is not used for production.
See prefix_benchmark.log, prefix_validation.json and runtime_implementation.json.

To export a human-readable count table after completion:

    /home/senescence/miniconda3/envs/kmer_validation_env/bin/python scripts/export_step1_counts.py --hap HG00096_hap1 > HG00096_hap1.step1.tsv

An optional `--kmer ACGTACGTACGTACGT` limits the export to one sequence, accepting
either orientation. A sequence absent from the preliminary universe is reported
as uncounted in background, not as an observed zero.

The prefix bitmap may be shared read-only across counting processes using an
adjacent `.prefix12.bitmap` file. The counter reconstructs the expected bitmap
from all query sequences and verifies the shared file byte-for-byte before using
it. A stale or mismatched bitmap causes an explicit error. The first production
jobs used private bitmaps; later jobs use the validated shared representation.
Both implement the same prefix predicate. `prefix_shared_validation.txt` records
the comparison. The latest patch and binary checksum are in the runtime metadata.
