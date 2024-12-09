library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

### import functions from helper.R
source("/users/sparthib/retina_lrs/code/04_dtu_dge_dte/helper.R")

comparison <- "ROs"
method <- "bambu"
# Define constants
base_dir <- "/users/sparthib/retina_lrs/processed_data"
output_dir <- file.path(base_dir, "dtu", method, comparison, "plots", "go")
input_file <- file.path(base_dir, "dtu", method, comparison, "DGE_DTU_DTE.tsv")


# Load data for different comparisons
dge_data <- list(
  D100_vs_D45 = load_dge_data(input_file, "B_RO_D100", "C_RO_D45"),
  D200_vs_D45 = load_dge_data(input_file, "A_RO_D200", "C_RO_D45"),
  D200_vs_D100 = load_dge_data(input_file, "A_RO_D200", "B_RO_D100")
)

# Generate gene lists
gene_lists <- lapply(dge_data, get_gene_list)

# Perform analyses and save results/plots
onts <- c("MF", "BP", "CC")
for (comparison in names(gene_lists)) {
  genelist <- gene_lists[[comparison]]
  for (ont in onts) {
    run_ora(genelist, ont, output_dir)
    run_ora_plots(genelist, ont, output_dir)
  }
}


