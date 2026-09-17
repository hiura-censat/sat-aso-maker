STEP 0 plots (PNG, 180 dpi)

01_step0_overview.png: input selection, target/background length, valid 16-mer starts and non-ACGT fraction. Includes CHM13.
02_hsat23_by_chromosome.png: per-chromosome HSat2/3 length distributions and eligible-haplotype denominators. CHM13 is a separate diamond.
03_haplotype_chromosome_heatmaps.png: eligibility and family lengths, with the same row order and color scale for both families.
heatmap_row_order.tsv: exact row identities from top to bottom.

No clustering or k-mer abundance analysis is performed here. Missing chromosomes are never treated as biological zeros.
Background means outside existing HSat2/3 annotations on eligible chromosomes, not confirmed non-CenSat sequence.
Reproduce from the project directory: Rscript scripts/plot_step0.R
