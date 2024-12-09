library(dplyr)
library(tibble)
library(readr)
library(grid)

counts_matrix_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/"

source("/users/sparthib/retina_lrs/code/05_visualization/helper.R")

load_and_plot_data <- function(method,compare, counts_matrix_dir) {
  isoform_tpm <- readRDS(file.path(counts_matrix_dir,
                                   method,
                                   compare,
                                   "isoform_cpm.RDS"))
  gene_tpm <- readRDS(file.path(counts_matrix_dir,
                                method,
                                compare,
                                "gene_cpm.RDS"))
  samples <- colnames(isoform_tpm)
  if (compare == "FT_vs_RGC"){
    groups <- c("FT", "FT", "RGC", "RGC")
  }else if (compare == "ROs"){
    groups <- c("RO_D45", "RO_D45", "RO_D100","RO_D100", 
                "RO_D100", "RO_D200", "RO_D200")
  }
  
  pca_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                 method, compare, "plots", "PCA")
  if (!dir.exists(pca_plots_dir)) {
    dir.create(pca_plots_dir, recursive = TRUE)
  }
  
  plot_pca(gene_tpm, samples, groups, "gene", pca_plots_dir)
  plot_pca(isoform_tpm, samples, groups, "isoform", pca_plots_dir)
  
}

load_and_plot_data("bambu", "ROs", counts_matrix_dir)
load_and_plot_data("bambu", "FT_vs_RGC", counts_matrix_dir)
load_and_plot_data("Isoquant", "ROs", counts_matrix_dir)
load_and_plot_data("Isoquant", "FT_vs_RGC", counts_matrix_dir)