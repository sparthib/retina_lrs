# processed_data/short_reads

Processed outputs of the short-read RNA-seq pipeline (BioProject PRJNA754196),
and the platform comparison against the long-read bambu quantification.

## Layout

- `01_correlation/` — short-read (salmon) vs long-read (bambu) sample-level
  Spearman correlation over all common protein-coding genes, matrix + heatmap.
- `02_platform_DE/` — platform differential expression at matched developmental
  stages (edgeR QLF): per-stage tables (W_S1/S2/S3) + summary.
- `03_length_gc/` — gene length and GC content of platform-DE up/down genes:
  per-gene annotation, distribution plot, and summary.

## Full count matrices (not version-controlled here, size)

  /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/
    gene_counts_matrix.csv        (62,266 genes x 28 samples)
    transcript_counts_matrix.csv  (251,955 transcripts x 28 samples)
    txi_salmon.rds, txi_salmon_transcript.rds

## Code and methods

- Pipeline: `code/11_short_reads_processing/` (numbered 01-06).
- DE + correlation scripts: `code/05_visualization/` (08b, 08d).
- Methods: `code/11_short_reads_processing/short_read_preprocessing_methods.md`.
