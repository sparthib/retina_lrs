library(edgeR)
library(sessioninfo)
library(biomaRt)
library(readr)

options(timeout = 999999)
us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)


method <- "bambu"
comparison <- "RO_vs_RGC"

# Set directories
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
dge_output_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "protein_coding", "DGE/")
dte_output_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison,"protein_coding",  "DTE/")
dir.create(dge_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dte_output_dir, showWarnings = FALSE, recursive = TRUE)
#### DTE Analysis ####
# Load data
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered_by_counts_and_biotype")
counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)


isoformFeatures <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                      method, comparison,"protein_coding",  "isoformFeatures.tsv"))

counts <- counts[rownames(counts) %in% isoformFeatures$isoform_id,]
nrow(counts)

groups  = c("Stage_1", "Stage_1", "Stage_2","Stage_2","Stage_2", 
            "Stage_3","Stage_3", "RGC", "RGC" )
y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

y <- normLibSizes(y)
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "protein_coding", "plots/")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file = file.path(plots_dir, "isoform_MDS.pdf"))
plotMDS(y, labels = y$samples$group, )
dev.off()
# Normalize and estimate dispersion
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))


y <- estimateDisp(y, design, robust=TRUE)
y <- estimateTagwiseDisp(y)

pdf(file = file.path(plots_dir, "isoform_BCV.pdf"))
plotBCV(y, col.common="red", col.trend="blue", col.tagwise="black")
dev.off()

# Fit the model and create contrasts
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(D45_vs_RGC = Stage_1 - RGC, 
                       D100_vs_RGC = Stage_2 - RGC, 
                       D200_vs_RGC = Stage_3 - RGC, 
                       levels=design)

# Create output directory for DGE results
dir.create(dte_output_dir, 
           showWarnings = FALSE, recursive = TRUE)


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
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == -1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == 1]
  
  file <- paste0(dte_output_dir, colnames(contr)[i], "_DTEs.tsv")
  write_tsv(tt, file)
}

# DGE Analysis

counts <- file.path(matrix_dir, "filtered_gene_counts.RDS") 
counts <- readRDS(counts)

y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))
y <- normLibSizes(y)

design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
y <- estimateDisp(y, design, robust=TRUE)
y <- estimateTagwiseDisp(y)

pdf(file = file.path(plots_dir, "gene_BCV.pdf"))
plotBCV(y, col.common="red", col.trend="blue", col.tagwise="black")
dev.off()

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
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == -1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == 1]
  
  file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  write_tsv(tt, file)
}

sessionInfo()