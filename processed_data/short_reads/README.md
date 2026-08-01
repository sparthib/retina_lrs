# processed_data/short_reads

Processed outputs of the short-read RNA-seq pipeline (BioProject PRJNA754196,
Wahlin lab; 28 retinal-organoid samples across D0-D280) and its quantifier-
controlled comparison against the long reads. The long-read samples are quantified
two ways: with **bambu** (the original full-length assignment) and — to remove the
quantifier confound — re-quantified with **salmon** in alignment mode against the
same release-46 transcriptome. Unless noted, the comparisons below use **salmon on
both platforms** (SR-salmon vs LR-salmon), so the only variable is platform, not
the counting algorithm.

## Layout

- `01_correlation/` — SR-salmon vs LR-salmon sample-level Spearman correlation, matrix + heatmap.
- `02_platform_DE/` — SR-salmon vs LR-salmon differential expression at matched stages, tables + summary + GO (`go/`).
- `03_length_gc/` — length/GC of the platform-DE genes.
- `04_isoforms/` — isoform-level comparison.
  - `04_isoforms/salmon_comparison/` — SR-salmon vs LR-salmon isoforms-per-gene comparison + provenance.
  - `04_isoforms/go/` — GO:BP of genes with isoforms detected only in LR-salmon (filtered).

## Short-read preprocessing and quantification

1. Download. Paired FASTQs for the 28 samples pulled from ENA (SRR2118 re-deposit
   series, Oligo-dT SMART-Seq v4 libraries). Non-UMI, so no deduplication.
2. Reference. GENCODE release 46 primary-assembly transcriptome — the same
   annotation as the long-read quantification, so Ensembl gene IDs are directly
   comparable across platforms.
3. Index. Plain (decoy-free) salmon index, `salmon index -k 31`, 251,955 targets.
4. Quantification. Each sample: `salmon quant -l A --validateMappings --gcBias
   --seqBias` (no dedup). Mean mapping rate 72.3% (range 65-80%).
5. Counts matrix. `tximport` aggregates the 28 quants to a gene matrix
   (62,266 genes x 28) and, with `txOut=TRUE`, a transcript matrix
   (251,955 x 28). Column names are the sample IDs (e.g. D100_5).

The long reads are re-quantified with salmon in alignment mode
(`salmon quant -a`) on the existing minimap2 transcriptome BAMs
(`05_bams/transcriptome/ver_46`, `map-ont`) against the same release-46 FASTA,
then `tximport`ed to gene (62,700 x 7) and transcript (252,835 x 7) matrices.

Pipeline scripts: `code/11_short_reads_processing/` (01 download, 02 manifests,
03 integrity, 04 salmon index, 05 salmon quant, 06 tximport, 11 LR salmon quant,
12 LR-salmon isoform matrix + comparison).

## Correlation (`01_correlation/`)

SR-salmon and LR-salmon gene expression are correlated at the sample level. Both
matrices are CPM-normalized and subset to the protein-coding genes found in common
between the two platforms (19,993 genes); a Spearman correlation is computed
between every SR sample and every LR sample -> a 28 x 7 matrix
(`SRsalmon_vs_LRsalmon_spearman_corr_allPC.csv`) and a heatmap
(`SRsalmon_vs_LRsalmon_correlation_heatmap.png`, white->pink->red, midpoint 0.5).
With the quantifier held constant, cross-platform agreement is high: rho ~0.76-0.90
(vs ~0.61-0.85 for the earlier SR-salmon vs LR-bambu comparison). SR samples
correlate most strongly with the developmentally matched LR stage.

## Platform differential expression (`02_platform_DE/`)

Platform DE contrasts SR-salmon vs LR-salmon at three matched developmental stages,
reported as a technology difference (not biology):

  W_S1  LR D45  vs  SR D25 + D65
  W_S2  LR D100 vs  SR D65 + D100
  W_S3  LR D200 vs  SR D180 + D280

Method: edgeR quasi-likelihood F-test (QLF) with a platform factor. TMM
normalization and the CPM>=1 expression filter are computed on the protein-coding
gene subset within each platform before testing. LR-salmon is subset to protein-
coding genes to match the biotype scope.

Significance: **FDR < 0.05 AND |log2FC| >= 1**. log2FC >= +1 = up in long reads,
<= -1 = up in short reads. `platform_DE_stage_summary.csv` reports the FDR-only
count (`sig_FDR05`), the FDR + fold-change count (`sig_FDR05_absLFC1`), and the
directional split (`up_in_LR_LFC1`, `up_in_SR_LFC1`):

  W_S1  FDR05 6323   FDR05+|LFC|>=1 4524   (up LR 2800 / up SR 1724)
  W_S2  FDR05 8836   FDR05+|LFC|>=1 5063   (up LR 3209 / up SR 1854)
  W_S3  FDR05 8592   FDR05+|LFC|>=1 6318   (up LR 3847 / up SR 2471)

Fewer genes are called than in the SR-salmon vs LR-bambu version (5160/5878/7229),
as expected once both platforms use the same quantifier. Per-stage full tables
(with gene_name + median tx length/GC): `platform_DE_stage_W_S1/S2/S3.csv`.
GO:BP of the up-in-LR / up-in-SR sets per stage: `go/`. up-in-LR is enriched for
cilium/sensory/visual-perception terms (W_S3: visual perception, sensory perception
of light stimulus); up-in-SR for carbohydrate/lipid metabolic terms.
Scripts: `code/05_visualization/08h_platform_DE_LRsalmon.R`, `08f_platform_DE_GO_analysis.R`.

## Length / GC of platform-DE genes (`03_length_gc/`)

Distributions of median gene GC% and median transcript length for genes up- vs
down-in-LR at each stage (same FDR<0.05 & |log2FC|>=1 rule). Genes higher in long
reads are consistently **longer and lower-GC** (median GC ~49%, length ~1346-1529 bp)
than genes higher in short reads (~53%, ~1181-1235 bp). Files: `length_gc_distributions.png`,
`length_gc_distribution_summary.csv`, `gene_length_gc.csv`.
Script: `code/05_visualization/08i_correlation_and_lengthgc_LRsalmon.R` (builds 01 + 03).

## Isoforms (`04_isoforms/`)

Isoforms are compared with salmon on both platforms. The long reads are
re-quantified with salmon (alignment mode) and both platforms go through the same
isoform pipeline: protein-coding subset + `filterByExpr(min.count=2)` + TMM.
Detection filter throughout: >= 1 count in >= 1 sample.

Provenance from salmon quantification to the plotted per-gene distribution
(`salmon_comparison/isoforms_per_gene_salmon_vs_bambu_provenance.csv`):

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
long reads are quantified with salmon (to 5.91). `salmon_comparison/` holds the
per-sample dot plot (SR = red, LR = blue), its exact plot data, and the provenance
table.

LR-salmon-only isoforms: isoforms in the LR-salmon filtered matrix but not the
SR-salmon filtered matrix = 15,316 isoforms across 7,500 genes (854 entirely
LR-specific). GO:BP of those genes (universe = 17,585 LR-salmon filtered genes) is
in `go/` (`GO_BP_LRsalmonOnly.csv/.png`, 99 terms — small-GTPase signaling,
peptidyl-serine phosphorylation, axonogenesis, and olfactory/chemosensory
detection). Gene + isoform lists: `go/LRsalmonOnly_genes.csv`, `LRsalmonOnly_isoforms.csv`.
Script: `code/05_visualization/08g_LRsalmon_only_isoform_GO_analysis.R`.
Per-sample SR transcript detection: `transcript_detection_per_sample.csv`, `_per_day.csv`.

## Full count matrices (not version-controlled, size)

  /dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/
    gene_counts_matrix.csv            (62,266 genes x 28 samples, SR-salmon)
    transcript_counts_matrix.csv      (251,955 transcripts x 28 samples, SR-salmon)
    LRsalmon_gene_counts_matrix.csv   (62,700 genes x 7 samples, LR-salmon)
    LRsalmon_transcript_counts_matrix.csv (252,835 transcripts x 7 samples, LR-salmon)
