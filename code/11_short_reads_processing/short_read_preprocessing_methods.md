# Short-read RNA-seq preprocessing: from raw FASTQ to counts matrix

**Data acquisition.** Short-read RNA-seq data for the retinal organoid
developmental time course were obtained from NCBI BioProject PRJNA754196
(Wahlin lab) — 28 biological samples across seven timepoints (D0, D10, D25,
D65, D100, D180, D280).

**Reference.** GENCODE release 46 primary assembly transcriptome — matching the
long-read bambu quantification so Ensembl gene IDs are directly comparable
across platforms.

**salmon index.** Decoy-free transcriptome index with salmon v1.10.1
(`salmon index -k 31`), 251,955 targets; tx2gene map from the pipe-delimited
FASTA headers.

**salmon quantification.** Each of the 28 samples quantified with
`salmon quant -l A --validateMappings --gcBias --seqBias`.

**Counts matrix (tximport).** Transcript estimates summarized to gene level
with tximport 1.37.2 (R 4.5.x) using the release-46 tx2gene map; version
suffixes stripped, columns renamed to D&lt;day&gt;_&lt;rep&gt; and
timepoint-ordered.

## Platform differential expression (08b)

**Gene universe.** The short-read gene-count matrix was intersected with the
long-read bambu gene-count matrix (`filtered_gene_counts.RDS`), which had already
been restricted to protein-coding genes (biomaRt `gene_biotype == "protein_coding"`)
and expression-filtered. Ensembl version suffixes were stripped on both sides so
the intersection is over shared ENSG identifiers; all normalization therefore
operated only on the shared protein-coding set.

**Matched developmental windows.** Long-read and short-read samples were compared
at matched developmental stages: Stage 1 = long-read D45 vs short-read D25+D65;
Stage 2 = long-read D100 vs short-read D65+D100; Stage 3 = long-read D200 vs
short-read D180+D280.

**Normalization and filtering.** Within each window, TMM normalization factors
were computed **separately** for the long-read and short-read samples
(edgeR `calcNormFactors` on each platform's sub-matrix) and combined into a single
`DGEList`. A depth-independent expression filter was then applied in **both**
platforms — a gene was retained only if CPM ≥ 1 in at least *n* samples of each
platform (*n* = the number of long-read samples in that window).

**Testing.** Differential expression on the platform factor
(`factor(platform, levels = c("SR","LR"))`) was tested with the edgeR
quasi-likelihood F-test: `estimateDisp` → `glmQLFit(~platform)` →
`glmQLFTest(coef = "platformLR")`. Positive log2 fold change denotes higher
expression in long reads. Significant platform-DE genes were called at
FDR < 0.05 and |log2FC| ≥ 1 (up in long read: log2FC ≥ +1; up in short read:
log2FC ≤ −1).

**Interpretation.** Platform is fully confounded with sequencing technology, so
this comparison measures technical concordance/discordance between the two
platforms, not biology, and cannot be batch-corrected (removing the platform
effect removes the comparison).

## Short-read vs long-read correlation matrix

Sample-level concordance was quantified over all protein-coding genes common to
both platforms (16,970 genes = the long-read filtered matrix intersected with the
short-read matrix, excluding the mitochondrial rRNA gene ENSG00000210082).
Short-read expression was TMM-normalized CPM computed on the protein-coding
subset (edgeR `calcNormFactors` → `cpm`); long-read expression was the bambu CPM
matrix (`filtered_gene_cpm.RDS`). For every short-read × long-read sample pair,
the Spearman rank correlation across the 16,970 shared genes was computed,
yielding a 28 × 7 correlation matrix (short-read samples × long-read samples).
