# STEP 0: HSat2/3 input preparation

Inputs are existing VallePrep Altemose annotations, T2T pass BEDs and filtered,
strand-oriented FASTA files. `manifest.tsv` records absolute source paths and
eligibility. Chromosome coordinates refer to the filtered, oriented FASTA, not
to the unoriented original assembly. Annotation lift-over adds the extraction
start offset; extraction was from this same filtered FASTA without strand reversal.
BED coordinates are zero-based, half-open. Source BED9 is preserved verbatim;
its thickStart/thickEnd fields were not lifted by the original pipeline and must
not be used as genome coordinates. Use the first three columns and strand column.

Only completed Altemose runs with valid coordinates are included. Empty or
incomplete runs are not biological zeros. An invalid annotation excludes the
whole hap rather than silently correcting its boundaries. See manifest.tsv and
invalid_annotations.tsv. CHM13, if eligible, is retained under its original ID.

Region BEDs contain merged HSat2 and HSat3 intervals, with cross-family overlaps
removed from those family-specific sets and saved as ambiguous. HSat23 is their
inclusive union. Background is its complement on all T2T-passed chromosomes.
No label for non-CenSat is inferred: background can include other satellite
families and HSat2/3 missed by the existing annotation. The Altemose pipeline
searched dna-brnn-selected, padded regions, not every base of each chromosome.
Consequently the background definition is annotation-based, not proof of absence.
Original per-occurrence orientation remains in source_annotation.bed; normalized
region views use the forward orientation of the filtered FASTA.

Sequence storage uses symlinks to existing FASTA and FAI plus region BEDs, avoiding
roughly 319 GB of redundant genomes. The links depend on the source remaining in
place. To materialize a region view (records are never concatenated):

    python3 scripts/extract_step0_sequences.py --hap HG00096_hap1 --region HSat2 | gzip > step0/HG00096_hap1.HSat2.fa.gz

region_stats.tsv counts bp, ACGT bp, other bp, and valid 16-mer starts separately
for every eligible hap/chromosome/region. Each maximal ACGT run contributes
max(0, length - 15); runs stop at region boundaries. Lowercase ACGT is accepted.
No reverse-complement collapsing or k-mer counting occurs in STEP 0.
The union's valid-window count need not equal the sum of the family counts
because windows crossing a family boundary can occur within the union.

Reproduce: python3 scripts/step0_preprocess.py --workers 4
Per-hap stats under qc/ are resumable checkpoints for this fixed input set.
For changed inputs or code, use a fresh output directory/remove old checkpoints
before rerunning. COMPLETE.json is written only after all sequence scans succeed.

## Completed run

- 574 eligible sequence sets: 573 haplotypes plus the CHM13 reference.
- 8,531 sequence-set/chromosome pairs; 42,655 region-statistics rows.
- Nine excluded inputs: eight incomplete/failed annotations and one invalid BED
  coordinate (`HG01786_hap2`, chr7:59222563-59222541).
- All HSat23/background partitions sum to their chromosome lengths; all family
  bp and ACGT totals agree with the union. Input sizes and modification times
  were unchanged during processing. See validation_summary.json.
- Disk usage is approximately 36 MB because genome sequences are source-linked.
- sequence_views.tsv explicitly maps each region view to its FASTA and BED.
- The preprocessing script uses the compiled scanner if available, otherwise
  the Python implementation. The two scanners were checked against synthetic
  boundary, ambiguous-base, lowercase and randomized-sequence fixtures.
- Validate a completed run with: python3 scripts/verify_step0.py
