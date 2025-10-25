library(dplyr)
library(ggrepel)
library(ggplot2)
library(edgeR)
library(grid)
library(patchwork)
library(pheatmap)


source("/users/sparthib/retina_lrs/code/05_visualization/helper.R")

gene_counts <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/H9_vs_EP1/filter_counts/ASE_H1_vs_H2_PTC_gene_counts.txt"
gene_counts <- read.table(gene_counts, header = TRUE, row.names = 1, sep = "\t")
colnames(gene_counts)

alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")

groups <- c("Stage_3", "Stage_3", 
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "Stage_2", "Stage_2",
            "Stage_3", "Stage_3",
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "FT", "FT",  
            "FT", "FT",
            "RGC", "RGC",
            "RGC", "RGC")

groups <- paste0(groups,"_", alleles)

stage_group <- gsub("_(H1|H2)$", "", groups)  # Remove _H1 or _H2 suffix
haplotype <- ifelse(grepl("H1$", groups), "H1", 
                    ifelse(grepl("H2$", groups), "H2", "Unknown"))


# > groups
# [1] "Stage3 H1" "Stage3 H2" "Stage2 H1" "Stage2 H2" "Stage1 H1" "Stage1 H2"
# [7] "Stage2 H1" "Stage2 H2" "Stage3 H1" "Stage3 H2" "Stage2 H1" "Stage2 H2"
# [13] "Stage1 H1" "Stage1 H2" "FT H1"     "FT H2"     "FT H1"     "FT H2"    
# [19] "RGC H1"    "RGC H2"    "RGC H1"    "RGC H2"   

samples <- c("Stage3 H1 1", "Stage3 H2 1", 
             "Stage2 H1 1", "Stage2 H2 1", 
             "Stage1 H1 1", "Stage1 H2 1", 
             "Stage2 H1 2", "Stage2 H2 2", 
             "Stage3 H1 2", "Stage3 H2 2", 
             "Stage2 H1 3", "Stage2 H2 3", 
             "Stage1 H1 2", "Stage1 H2 2", 
             "FT H1 1",     "FT H2 1",     
             "FT H1 2",     "FT H2 2",    
             "RGC H1 1",    "RGC H2 1",    
             "RGC H1 2",    "RGC H2 2")

# convert counts to cpm
dge <- DGEList(counts = gene_counts, group = groups, samples = samples)
dge <- calcNormFactors(dge, method = "TMM")
cpm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)


filter_gene_counts <- function(counts_matrix, group ){ 
  min_counts <- 10
  dge <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(dge, group = group, min.count = min_counts)
  dge <- dge[keep, ]
  dge <- calcNormFactors(dge)
  # Convert to CPM (Counts Per Million)
  cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
  return(list(filtered_counts = counts_matrix[keep, ], cpm = cpm_matrix))
  
}

filtered_results <- filter_gene_counts(gene_counts, groups)

gene_cpm <- filtered_results$cpm


RO_gene_counts <- gene_counts[,1:14]
RO_gene_cpm <- gene_cpm[,1:14]
RO_groups <- groups[1:14]
RO_samples <- samples[1:14]

pca_plots_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/pca_plots/"
dir.create(pca_plots_dir, showWarnings = FALSE)
plot_pca(gene_cpm, samples, groups, "gene", pca_plots_dir)
plot_pca(RO_gene_cpm, RO_samples, RO_groups, "RO_gene", pca_plots_dir)

# Plot correlation heatmap
colnames(gene_cpm) <- samples
cor_mat <- cor(gene_cpm, method = "spearman")

colSums(gene_cpm)

pdf(paste0(pca_plots_dir, "sample_correlation_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  cor_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  display_numbers = TRUE,
  color = colorRampPalette(c("white", "pink", "red"))(100),
  main = "Sample Correlation Heatmap"
)
dev.off()







