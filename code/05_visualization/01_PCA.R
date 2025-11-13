library(dplyr)
library(tibble)
library(readr)
library(grid)
library(patchwork)
library(ggrepel)
library(ggplot2)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")


counts_matrix_dir <- file.path(data_dir, "06_quantification/counts_matrices")

source(file.path(code_dir,"code/05_visualization/helper.R"))

load_and_plot_data <- function(method,compare, counts_matrix_dir) {
  isoform_tpm <- readRDS(file.path(counts_matrix_dir,
                                   method,
                                   compare, "filtered_by_counts_and_biotype",
                                   "filtered_isoform_cpm.RDS"))
  gene_tpm <- readRDS(file.path(counts_matrix_dir,
                                method,
                                compare, "filtered_by_counts_and_biotype",
                                "filtered_gene_cpm.RDS"))
  
  isoform_tpm <- remove_zero_var_rows(isoform_tpm)
  gene_tpm <- remove_zero_var_rows(gene_tpm)
  
  if (compare == "FT_vs_RGC"){
    samples <- c("FT_1", "FT_2", "RGC_1", "RGC_2")
    groups <- c("FT", "FT", "RGC", "RGC")
  }else if (compare == "ROs"){
    samples <- c("Stage_1_1", "Stage_1_2", "Stage_2_1", 
                 "Stage_2_2", "Stage_2_3", "Stage_3_1", "Stage_3_2")
    groups <- c("Stage_1", "Stage_1", "Stage_2","Stage_2", 
                "Stage_2", "Stage_3", "Stage_3")  }
    else if (compare == "RO_vs_RGC"){ 
      samples <- c("Stage_1_1", "Stage_1_2", "Stage_2_1", 
                   "Stage_2_2", "Stage_2_3", "Stage_3_1", "Stage_3_2",
                   "RGC_1", "RGC_2")
      groups <- c("Stage_1", "Stage_1", "Stage_2","Stage_2", 
                  "Stage_2", "Stage_3", "Stage_3", "RGC", "RGC")
      }

  
  pca_plots_dir <- file.path(code_dir,"processed_data/dtu/",
                                 method, compare, "protein_coding",  "plots", "PCA")
  if (!dir.exists(pca_plots_dir)) {
    dir.create(pca_plots_dir, recursive = TRUE)
  }
  
  plot_pca(gene_tpm, samples, groups, "gene", pca_plots_dir)
  plot_pca(isoform_tpm, samples, groups, "isoform", pca_plots_dir)
  
}

load_and_plot_data("bambu", "ROs", counts_matrix_dir)
load_and_plot_data("bambu", "FT_vs_RGC", counts_matrix_dir)
load_and_plot_data("bambu", "RO_vs_RGC", counts_matrix_dir)
