library(edgeR)
library(sessioninfo)
library(biomaRt)
library(readr)

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)


method <- "bambu"
comparison <- "ROs"

# Set directories
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
dge_output_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE/")
dte_output_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE/")

#### DTE Analysis ####
# Load data
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison)
counts <- file.path(matrix_dir, "isoform_counts.RDS") 
counts <- readRDS(counts)


isoformFeatures <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                      method, comparison, "isoformFeatures.tsv"))

counts <- counts[rownames(counts) %in% isoformFeatures$isoform_id,]
nrow(counts)

groups  = c("C_RO_D45", "C_RO_D45", "B_RO_D100","B_RO_D100","B_RO_D100", "A_RO_D200","A_RO_D200" )
y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

y <- normLibSizes(y)

# Normalize and estimate dispersion
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- estimateDisp(y, design, robust=TRUE)

# Fit the model and create contrasts
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(D200_vs_D100 = A_RO_D200 - B_RO_D100, 
                       D200_vs_D45 = A_RO_D200 - C_RO_D45, 
                       D100_vs_D45 = B_RO_D100 - C_RO_D45, 
                       levels=design)

# Create output directory for DGE results
dir.create(dte_output_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$isoform_id <- tt$genes
  
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

# DGE Analysis
dir.create(dge_output_dir, showWarnings = FALSE, recursive = TRUE)

counts <- file.path(matrix_dir, "gene_counts.RDS") 
counts <- readRDS(counts)

y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

keep <- filterByExpr(y, )
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)

design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)


# Loop through contrasts and save DGE results
for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$gene_id <- gsub("\\..*", "", tt$genes)
  
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


sessionInfo()
