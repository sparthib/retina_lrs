library(edgeR)
library(sessioninfo)
# Set directories
isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE/"
dte_output_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DTE/"

#### DGE Analysis ####
# Load data
counts <- read.table(file.path(isoquant_dir, "OUT.gene_grouped_counts.tsv"), header = TRUE)
tpm <- read.table(file.path(isoquant_dir, "OUT.gene_grouped_tpm.tsv"), header = TRUE)

# Define sample names and subset data for RO samples
sample_names <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                  "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
                  "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")
colnames(counts) <- c("gene_id", sample_names)
colnames(tpm) <- c("gene_id", sample_names)

# Select RO samples and define group factors
RO_counts <- counts[,1:8]
RO_tpm <- tpm[,1:8]
group <- factor(c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"))

# Create DGEList object and filter
y <- DGEList(counts=RO_counts, group=group, genes=RO_counts$gene_id, samples=sample_names[1:7])
keep <- filterByExpr(y, min.count = 3)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalize and estimate dispersion
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- estimateDisp(y, design, robust=TRUE)

# Fit the model and create contrasts
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(D200_vs_D100 = RO_D200 - RO_D100, 
                       D200_vs_D45 = RO_D200 - RO_D45, 
                       D100_vs_D45 = RO_D100 - RO_D45, 
                       levels=design)

# Create output directory for DGE results
dir.create(dge_output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through contrasts and save DGE results
for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$gene_id <- gsub("\\..*", "", tt$gene_id)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
                       filter="ensembl_gene_id", 
                       values=tt$gene_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype")
  
  tt <- merge(tt, annotLookup, by="gene_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == 1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == -1]
  
  file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  write_tsv(tt, file)
}

# DTE Analysis (similar to DGE)
dir.create(dte_output_dir, showWarnings = FALSE, recursive = TRUE)

counts <- read.table(file.path(isoquant_dir, "OUT.transcript_grouped_counts.tsv"), header = TRUE)
tpm <- read.table(file.path(isoquant_dir, "OUT.transcript_grouped_tpm.tsv"), header = TRUE)

colnames(counts) <- c("isoform_id", sample_names)
colnames(tpm) <- c("isoform_id", sample_names)

RO_counts <- counts[,1:8]
RO_tpm <- tpm[,1:8]

y <- DGEList(counts=RO_counts, group=group, genes=RO_counts$isoform_id, 
             samples=sample_names[1:7])
keep <- filterByExpr(y, min.count = 1)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$isoform_id <- gsub("\\..*", "", tt$isoform_id)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "transcript_biotype"), 
                       filter="ensembl_transcript_id", 
                       values=tt$isoform_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "isoform_id", "gene_name", "transcript_biotype")
  
  tt <- merge(tt, annotLookup, by="isoform_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == 1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == -1]
  
  file <- paste0(dte_output_dir, colnames(contr)[i], "_DTEs.tsv")
  write_tsv(tt, file)
}

sessionInfo()




