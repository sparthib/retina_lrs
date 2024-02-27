# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)

# tidyverse-friendly packages
library(plotly)
library(ggrepel)
library(GGally)
library(tidybulk)
library(tidyHeatmap)
library(tidySummarizedExperiment)
colData <- data.frame(group = c("RGC", "RGC", "RGC","RGC"),
                      sample_name=c("DG-WT-hRGC", "hRGC","YZ-15T_hRGC","YZ-3T_hRGC"))

geneCountMatrix <- read.csv("/users/sparthib/retina_lrs/processed_data/de/gene_counts_matrix_Jan13.csv")
rownames(geneCountMatrix) <- geneCountMatrix$gene_name


geneCountMatrix <- geneCountMatrix[,c(4,9, 10,11)]

geneCountMatrix <- as.matrix(geneCountMatrix)

se <- SummarizedExperiment(assays=list(counts=geneCountMatrix),
                           colData=colData)

se |> group_by(`.sample`) |> summarise(n=n())
# # A tibble: 4 × 2
# .sample         n
# <chr>       <int>
#   1 DG.WT.hRGC  55910
# 2 YZ.15T_hRGC 55910
# 3 YZ.3T_hRGC  55910
# 4 hRGC        55910



# We filter out lowly expressed genes using tidybulk `keep_abundant`
# or `identify_abundant`. 
# These functions can use the edgeR filterByExpr function described 
# in (Law et al. 2016) to automatically identify the genes with adequate abundance 
# for differential expression testing.

# Pre-processing


# Filtering counts
counts_scaled <- se |> tidybulk::identify_abundant(factor_of_interest = group) |>  tidybulk::scale_abundance()

# take a look
counts_scaled


######### Clustering ###########
counts_scaled_top_500 <- counts_scaled |> # extract 500 most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) |> as_tibble()

counts_scaled_top_500 |> heatmap( .column = `.sample`,
                                  .row = `.feature`,
                                  .value = counts_scaled,
                                  transform = log1p,
                                  row_names_gp = grid::gpar(fontsize = 0.8),
                                  column_names_gp = grid::gpar(fontsize = 5),
                                  column_title_gp = grid::gpar(fontsize = 7),
                                  row_title_gp = grid::gpar(fontsize = 7)) |> 
  add_tile(group,
           palette = c("#CC6677","#332288")) |>
  save_pdf("/users/sparthib/retina_lrs/plots/eda/tidy_bulk/annotated_cluster_RGC_Jan30.pdf")


