library(dplyr)
library(ggrepel)
library(ggplot2)
library(grid)
library(patchwork)
library(pheatmap)
library(readr)

source("/users/sparthib/retina_lrs/code/05_visualization/helper.R")

gene_cpm <- read_tsv("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_cpm.tsv")

samples <- colnames(gene_cpm)
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

# samples <- c("Stage3 H1 1", "Stage3 H2 1", 
#              "Stage2 H1 1", "Stage2 H2 1", 
#              "Stage1 H1 1", "Stage1 H2 1", 
#              "Stage2 H1 2", "Stage2 H2 2", 
#              "Stage3 H1 2", "Stage3 H2 2", 
#              "Stage2 H1 3", "Stage2 H2 3", 
#              "Stage1 H1 2", "Stage1 H2 2", 
#              "FT H1 1",     "FT H2 1",     
#              "FT H1 2",     "FT H2 2",    
#              "RGC H1 1",    "RGC H2 1",    
#              "RGC H1 2",    "RGC H2 2")

RO_gene_cpm <- gene_cpm[,1:14]
RO_groups <- groups[1:14]
RO_samples <- samples[1:14]

pca_plots_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/pca_plots/bambu"
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
