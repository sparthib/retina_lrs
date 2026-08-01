# processed_data/short_reads

Processed outputs of the short-read RNA-seq pipeline (BioProject PRJNA754196,
Wahlin lab; 28 retinal-organoid samples across D0-D280) and its comparison
against the long-read bambu quantification.

## Layout

- `01_correlation/` — SR vs LR sample-level Spearman correlation, matrix + heatmap.
- `02_platform_DE/` — platform differential expression at matched stages, tables + summary.
- `03_length_gc/` — length/GC of platform-DE genes.
- `04_isoforms/` — isoform detection and SR-vs-LR isoform comparison (see that folder's outputs).
  - `04_isoforms/salmon_comparison/` — quantifier-controlled SR-salmon vs LR-salmon isoforms-per-gene comparison.

## Short-read preprocessing and quantification

1. Download. Paired FASTQs for the 28 samples pulled from ENA (SRR2118 re-deposit
   series, Oligo-dT SMART-Seq v4 libraries). Non-UMI, so no deduplication.
2. Reference. GENCODE release 46 primary-assembly transcriptome — the same
   annotation as the long-read bambu quantification, so Ensembl gene IDs are
   directly comparable across platforms.
3. Index. Plain (decoy-free) salmon index, `salmon index -k 31`. The index holds 251,955 targets (transcript sequences) = every isoform in the release-46 annotation.
4. Quantification. Each sample: `salmon quant -l A --validateMappings --gcBias
   --seqBias` (no dedup). Mean mapping rate 72.3% (range 65-80%).
5. Counts matrix. `tximport` aggregates the 28 quants to a gene matrix
   (62,266 genes x 28) and, with `txOut=TRUE`, a transcript matrix
   (251,955 x 28). Column names are the sample IDs (e.g. D100_5).

Pipeline scripts: `code/11_short_reads_processing/` (01 download, 02 manifests,
03 integrity, 04 salmon index, 05 salmon quant, 06 tximport).
Full methods: `code/11_short_reads_processing/short_read_preprocessing_methods.md`.

## Correlation (`01_correlation/`)

SR (salmon) and LR (bambu) gene expression are correlated at the sample level.
Both matrices are subset to the protein-coding genes found in common between the
two platforms (16,970 genes), CPM-normalized, and a Spearman correlation is
computed between every SR sample and every LR sample. The result is a 28 x 7
correlation matrix (`SRsalmon_vs_LR_spearman_corr_allPC.csv`) and a heatmap
(`SRsalmon_LR_correlation_heatmap.png`). SR samples correlate most strongly with
the developmentally matched LR timepoint (early SR -> LR D45, mid -> LR D100,
late -> LR D200), rho ~0.61-0.85.

## Platform differential expression (`02_platform_DE/`)

Platform DE contrasts SR vs LR at three matched developmental stages, reported
as a technology difference (not biology):

  W_S1  LR D45  vs  SR D25 + D65
  W_S2  LR D100 vs  SR D65 + D100
  W_S3  LR D200 vs  SR D180 + D280

Method: edgeR quasi-likelihood F-test (QLF) with a platform factor. TMM
normalization and the CPM>=1 expression filter are computed on the protein-coding
gene subset within each platform before testing.

Significance: a gene is called differentially expressed at **FDR < 0.05 AND
|log2FC| >= 1**. Genes with log2FC >= +1 are up in long reads, log2FC <= -1 are
up in short reads. `platform_DE_stage_summary.csv` reports both the FDR-only
count (`sig_FDR05`) and the FDR + fold-change count (`sig_FDR05_absLFC1`),
plus the directional split (`up_in_LR_LFC1`, `up_in_SR_LFC1`):

  W_S1  FDR05 6893   FDR05+|LFC|>=1 5160   (up LR 3535 / up SR 1625)
  W_S2  FDR05 9182   FDR05+|LFC|>=1 5878   (up LR 4094 / up SR 1784)
  W_S3  FDR05 9364   FDR05+|LFC|>=1 7229   (up LR 4680 / up SR 2549)

Per-stage full tables: `platform_DE_stage_W_S1/S2/S3.csv`.
Scripts: `code/05_visualization/` (08b DE, 08d length/GC, 08e isoform GO).

## Isoform quantification and detection (`04_isoforms/`)

Isoforms are quantified at the transcript level on both platforms and compared
over the known GENCODE (ENST) isoforms. Detection filter: >= 1 estimated count
in >= 1 sample.

Detected transcripts after the >= 1 count filter:

| | Long read (7 samples) | Short read (28 samples) |
|---|---|---|
| Detected transcripts | 160,887 | 203,124 |
| — known ENST | 160,392 | 203,124 |
| — novel (Bambu) | 495 | 0 |

Overlap of detected known (ENST) isoforms between platforms:

| Set | Count |
|---|---|
| Detected in both | 141,906 |
| LR-only (detected in LR, not SR) | 18,486 |
| SR-only (detected in SR, not LR) | 61,218 |
| Novel (Bambu), LR-only by construction | 495 |

Only known ENST isoforms are comparable: the salmon index holds known GENCODE
transcripts only, so short reads cannot report novel isoforms, whereas the bambu
long-read annotation adds novel isoforms (Bambu... IDs). Outputs (LR-only isoform
lists, per-gene summaries, GO enrichment) are in `04_isoforms/`.

## LR-salmon quantification and the salmon-vs-salmon isoform comparison (`04_isoforms/salmon_comparison/`)

The isoform-per-gene difference between platforms is largely a **quantifier**
effect, not a read-length effect. To isolate it, the 7 retinal-organoid long-read
samples were re-quantified with salmon in **alignment mode** (`salmon quant -a`)
against the *same* release-46 transcriptome as the short reads, using the existing
minimap2 transcriptome BAMs (`05_bams/transcriptome/ver_46`, `map-ont`). Salmon's
EM then apportions long reads across compatible isoforms exactly as it does for
short reads, so the only variable left is platform (SR vs LR), not the algorithm.

Both platforms are then put through the identical isoform pipeline: protein-coding
subset + `filterByExpr(min.count=2)` + TMM (the same `01c` recipe used for gene
matrices). Detection filter throughout: >= 1 count in >= 1 sample.

Provenance from salmon quantification to the plotted per-gene distribution
(`isoforms_per_gene_salmon_vs_bambu_provenance.csv`):

| Stage | SR-salmon (28s) | LR-salmon (7s) | LR-bambu (7s) |
|---|---|---|---|
| 1. Annotation / index target | 251,955 | 252,835 | 278,023 (GTF, +1,118 novel) |
| 2. Quantifier | salmon quasi-map + EM | salmon aln-mode + EM | bambu 1-read to 1-isoform |
| 3. Isoforms detected (>=1 in >=1) | 203,124 | 226,256 | 160,887 |
| 4. Genes with those isoforms | 45,533 | 58,248 | 48,088 |
| 5. Mean isoforms/gene (unfiltered) | 4.46 | 3.88 | 3.35 |
| 6. Isoforms after PC + min.count=2 | 117,762 | 103,949 | 55,708 |
| 7. Genes in filtered matrix | 17,738 | 17,585 | 15,683 |
| 8. Mean isoforms/gene (filtered) | 6.64 | 5.91 | 3.55 |
| 9. Median / q90 / max (filtered) | 5 / 14 / 116 | 4 / 12 / 106 | 3 / 7 / 39 |
| 10. Genes >=10 isoforms (filtered) | 3,908 (22%) | 3,043 (17%) | 788 (5%) |

Held constant (salmon on both), SR and LR give near-identical gene universes
(17,738 vs 17,585 filtered genes) and closely tracking isoform-per-gene
distributions; the bambu-vs-salmon gap (mean 3.55 vs 6.64) closes ~76% when the
long reads are quantified with salmon (to 5.91). The apparent "short reads detect
more isoforms per gene" was mostly the bambu-vs-salmon assignment difference.

Files:
- `isoforms_per_gene_salmon_SR_vs_LR_persample.png` — per-sample dot plot (SR = red,
  LR = blue), isoforms/gene binned 1-10 and 10+, filtered matrices.
- `isoforms_per_gene_salmon_SR_vs_LR_persample_plotdata.csv` — exact plotted data
  (platform, sample, stage, bin, n_genes).
- `isoforms_per_gene_salmon_vs_bambu_provenance.csv` — the funnel table above.

Scripts: `code/11_short_reads_processing/11_salmon_quant_longreads.sh` (LR salmon
alignment-mode quant), `12_LRsalmon_isoform_matrix_and_comparison.R` (filtering +
comparison).

## Full count matrices (not version-controlled, size)

  /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/
    gene_counts_matrix.csv        (62,266 genes x 28 samples)
    transcript_counts_matrix.csv  (251,955 transcripts x 28 samples)
