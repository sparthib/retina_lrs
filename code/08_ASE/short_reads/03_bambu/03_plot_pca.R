library(dplyr)
library(ggrepel)
library(ggplot2)
library(grid)
library(patchwork)
library(pheatmap)
library(readr)

source("/users/sparthib/retina_lrs/code/05_visualization/helper.R")

gene_cpm <- read_tsv("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_cpm.tsv")
rownames(gene_cpm) <- gene_cpm$gene_id
gene_cpm <- gene_cpm %>% select(-gene_id)

samples <- c("EP1_Stage3_rep1_H1", "EP1_Stage3_rep1_H2", 
             "EP1_Stage2_rep1_H1", "EP1_Stage2_rep1_H2", 
             "EP1_Stage1_rep1_H1", "EP1_Stage1_rep1_H2", 
             "H9_Stage2_rep2_H1", "H9_Stage2_rep2_H2", 
             "H9_Stage3_rep2_H1", "H9_Stage3_rep2_H2", 
             "H9_Stage2_rep3_H1", "H9_Stage2_rep3_H2", 
             "H9_Stage1_rep2_H1", "H9_Stage1_rep2_H2", 
             "H9_FT_rep1_H1",     "H9_FT_rep1_H2",     
             "H9_FT_rep2_H1",     "H9_FT_rep2_H2",    
             "H9_RGC_rep1_H1",    "H9_RGC_rep1_H2",    
             "H9_RGC_rep2_H1",    "H9_RGC_rep2_H2")

colnames(gene_cpm) <- samples


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
groups <- paste0(groups, "_", alleles)

RO_gene_cpm <- gene_cpm[,1:14]
RO_groups <- groups[1:14]
RO_samples <- samples[1:14]

pca_plots_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/pca_plots/bambu"
dir.create(pca_plots_dir, showWarnings = FALSE)
plot_pca(gene_cpm, samples, groups, "gene", pca_plots_dir)
plot_pca(RO_gene_cpm, RO_samples, RO_groups, "RO_gene", pca_plots_dir)
