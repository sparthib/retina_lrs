# 01_counts

The full short-read count matrices are not version-controlled here (size).
They live on the data filesystem:

  /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/
    gene_counts_matrix.csv        (62,266 genes x 28 samples)
    transcript_counts_matrix.csv  (251,955 transcripts x 28 samples)
    txi_salmon.rds, txi_salmon_transcript.rds

See code/11_short_reads_processing/ for the pipeline that produces them and
short_read_preprocessing_methods.md for the methods.
