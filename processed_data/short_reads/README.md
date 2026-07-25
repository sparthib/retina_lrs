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
- `04_isoforms/` — isoform (transcript) detection in the short reads, the
  distribution of detected isoforms per gene, the short-read vs long-read
  isoform overlap (with the isoforms found only in long reads), and GO:BP
  enrichment of the long-read-only isoform genes (`go/`).

## Isoform (transcript) matrices

Both platforms produce a transcript-level counts matrix. Neither matrix is
version-controlled here (size); both live on the data filesystem:

  Short read (salmon + tximport, txOut=TRUE):
    /dcs04/.../10_short_reads/salmon/transcript_counts_matrix.csv
      251,955 transcripts x 28 samples  (full GENCODE release-46 annotation)

  Long read (bambu, extended annotation, retinal-organoid samples):
    /dcs04/.../06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS
      277,401 transcripts x 7 samples   (common-isoform-filtered detection matrix)

### How the number of isoforms was identified

The row count of a transcript matrix is the size of the *annotation* it was
quantified against, NOT the number of isoforms actually expressed. "Detected"
means an annotated transcript with expression evidence. Detection rule used
throughout (both platforms): >= 1 estimated count in >= 1 sample. Transcript
IDs are version-stripped (ENST...) before any comparison.

Short read (28 samples), from `transcript_counts_matrix.csv`:
  - in annotation (salmon index)          : 251,955
  - detected (>=1 count in >=1 sample)    : 203,124  (80.6%)
  - >=1 count in >=3 samples (reproducible): 169,671
  - >=5 count in >=3 samples (confident)  : 124,701
  Short reads were quantified against a PLAIN GENCODE release-46 transcriptome
  index, so only known ENST isoforms exist on this side — no novel isoforms.

Long read (7 retinal-organoid samples), traced through the bambu pipeline:
  - raw bambu counts_transcript.txt (annotation) : 278,023  (276,905 ENST + 1,118 novel Bambu)
  - after common-isoform filter (filtered/)      : 277,401  (276,905 ENST +   496 novel)
  - detected (>=1 count in >=1 of 7 samples)     : 160,887  (160,392 ENST +   495 novel)
  - protein-coding + expression filter (DE set)  :  55,708  ( 55,366 ENST +   342 novel)
  Long reads used an EXTENDED annotation = known GENCODE isoforms PLUS novel
  isoforms bambu discovered (Bambu... IDs), which short reads cannot detect by
  construction. ~117k annotated transcripts had zero reads in the 7 RO samples,
  which is why detected (160,887) is far below the annotation (277,401).

Note: short reads "detect" more total isoforms (203k vs 161k) partly because
there are 28 SR samples vs 7 LR samples (more libraries = more chances to clear
>=1 count), and because short reads spread multi-mapping reads across a gene's
annotated isoforms. This is isoform QUANTIFICATION against the annotation, not
isoform DISCOVERY; the long-read data remains the authority for true isoform
structure and for genuinely novel isoforms.

### Short-read vs long-read comparison (known ENST isoforms)

Only known GENCODE (ENST) isoforms are comparable between platforms. LR-only =
detected in long reads but not in the short-read quant.

  all known isoforms      : SR 203,124  LR 160,392  both 141,906  LR-only 18,486  SR-only 61,218
  protein-coding isoforms : SR  79,503  LR  63,974  both  60,310  LR-only  3,664  SR-only 19,193
  plus 495 novel (Bambu) isoforms detected only in long reads (no ENST equivalent).

The LR-only isoforms are annotated with gene_id / gene_name and summarized per
gene (`isoform_LRonly_*_by_gene.csv`): 13,801 genes carry an LR-only known
isoform, 2,776 genes an LR-only protein-coding isoform.

### 04_isoforms/ contents

- `transcript_detection_per_sample.csv` / `_per_day.csv` — detected isoforms per sample / per timepoint (union).
- `isoforms_per_gene_detected.csv` — detected vs annotated isoforms per gene, and detected fraction.
- `isoforms_per_gene_distribution.csv` / `.png` — distribution of detected isoforms per gene.
- `isoform_SR_LR_overlap_summary.csv` — the SR vs LR overlap table above.
- `isoform_LRonly_all_known.csv` / `isoform_LRonly_protein_coding.csv` — LR-only isoforms (transcript_id, biotype, gene_id, gene_name).
- `isoform_LRonly_all_known_by_gene.csv` / `isoform_LRonly_protein_coding_by_gene.csv` — per-gene counts of LR-only isoforms.
- `go/GO_BP_LRonly_{all_known,protein_coding}.{png,csv}` — GO:BP over-representation (clusterProfiler enrichGO) of LR-only isoform genes.

## Full count matrices (not version-controlled here, size)

  /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/
    gene_counts_matrix.csv        (62,266 genes x 28 samples)
    transcript_counts_matrix.csv  (251,955 transcripts x 28 samples)
    txi_salmon.rds, txi_salmon_transcript.rds

## Code and methods

- Pipeline: `code/11_short_reads_processing/` (numbered 01-08:
  download, manifests, integrity, salmon index, salmon quant, tximport,
  isoform detection, SR-LR isoform comparison).
- DE + correlation + GO scripts: `code/05_visualization/` (08b, 08d, 08e).
- Methods: `code/11_short_reads_processing/short_read_preprocessing_methods.md`.
