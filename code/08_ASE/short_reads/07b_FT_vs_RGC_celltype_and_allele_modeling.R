library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

samples <- c("H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2") 
gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gene_counts_all_samples"

# load all counts matrix in the directory
files <- list.files(gene_counts_dir,
                    pattern = "_counts.txt", 
                    full.names = TRUE)

# create a list of data frames

# Load all the count data into a list of data frames
counts_list <- lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1)  # Adjust if needed (e.g., use sep = "," for CSV)
})

# Combine all data frames into a single matrix by column binding them
counts_matrix <- do.call(cbind, counts_list)

## remove rows with 0 values

counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]

colnames(counts_matrix) <- gsub(".bam", "", colnames(counts_matrix))


#remove version number from row names
rownames(counts_matrix) <- gsub("\\..*", "", rownames(counts_matrix))

# only keep PTC genes

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( 
    "external_gene_name",
    "gene_biotype", "ensembl_gene_id"),
  filter="ensembl_gene_id",
  values=rownames(counts_matrix),
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id")
annotLookup <- annotLookup |> dplyr::distinct()
annotLookup |> nrow()
annotLookup <- annotLookup |> 
  dplyr::filter(gene_biotype == "protein_coding") 
annotLookup |> nrow()

# filter for protein coding genes
counts_matrix <- counts_matrix[rownames(counts_matrix) 
                               %in% annotLookup$gene_id, ]

#FT vs RGC samples 
counts_matrix <- counts_matrix[, 9:16]

groups <- c("H1_FT", "H2_FT", "H1_FT", "H2_FT",
            "H1_RGC", "H2_RGC", "H1_RGC", "H2_RGC")

library(edgeR)

y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             group = groups,
             genes = rownames(counts_matrix))

y <- normLibSizes(y)

design <- model.matrix(~ 0 + group,data = y$samples)

colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

# Fit the model and create contrasts
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(H1_FT_vs_H2_FT = H2_FT - H1_FT, 
                       H1_RGC_vs_H2_RGC = H2_RGC - H1_RGC, 
                       H1_FT_vs_H1_RGC = H1_RGC - H1_FT,
                       H2_FT_vs_H2_RGC = H2_RGC - H2_FT,
                       levels=design)

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/"


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
  tt$neg_log10_FDR <- -log10(tt$FDR)
  tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
  tt <- tt |>
    arrange(desc(abs(logFC)), FDR)
  tt$label <- NA
  tt$label[1:20] <- tt$gene_name[1:20]
  file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  write_tsv(tt, file)
}

i = 4
for (i in seq_len(ncol(contr))) {
  read_file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  tt <- read_tsv(read_file)
  pdf(file = paste0(dge_output_dir, colnames(contr)[i], "_volcano_plot.pdf"),
      width = 8, height = 6)
  ggplot(tt, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = colnames(contr)[i],
         x = "log2 Fold Change",
         y = "-log10 FDR") +
    theme_minimal()
  dev.off()
}
