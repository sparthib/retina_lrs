library(dplyr)
library(ggplot2)
library(here)
library(httr)
library(readr)
library(jsonlite)
library(stringr)
library(reshape2)

url <- "https://raw.githubusercontent.com/WahlinLab/Organoid_RNAseq_SciData22/refs/heads/main/Normalized_counts_D0-D280.tabular"
short_read_data <- read.table(url, header = TRUE, sep = "\t", row.names = 1)
# Extract the day numbers from the column names
day_order <- c(0, 10, 25, 65, 100, 180, 280)
day_extracted <- as.numeric(stringr::str_extract(colnames(short_read_data), "\\d+"))
column_order <- order(factor(day_extracted, levels = day_order))
short_read_data <- short_read_data[, column_order]
saveRDS(short_read_data,
        file = here("raw_data/agarwal_short_norm_counts.RDS"))

short_read_data <- readRDS(here("raw_data/agarwal_short_norm_counts.RDS"))

#### DGE Comparisons ####
# 
# # Define the base URL for raw files on GitHub
# raw_base_url <- "https://raw.githubusercontent.com/WahlinLab/Organoid_RNAseq_SciData22/main/DESeq2_TimePoint_Comparisons/"
# 
# # List of tabular files (manually listed as per your input)
# tabular_files <- c(
#   "TimePoint_D0_vs_D10.tabular",
#   "TimePoint_D0_vs_D100.tabular",
#   "TimePoint_D0_vs_D180.tabular",
#   "TimePoint_D0_vs_D25.tabular",
#   "TimePoint_D0_vs_D280.tabular",
#   "TimePoint_D0_vs_D65.tabular",
#   "TimePoint_D100_vs_D180.tabular",
#   "TimePoint_D10_vs_D100.tabular",
#   "TimePoint_D10_vs_D180.tabular",
#   "TimePoint_D10_vs_D25.tabular",
#   "TimePoint_D10_vs_D280.tabular",
#   "TimePoint_D10_vs_D65.tabular",
#   "TimePoint_D25_vs_D100.tabular",
#   "TimePoint_D25_vs_D180.tabular",
#   "TimePoint_D25_vs_D280.tabular",
#   "TimePoint_D25_vs_D65.tabular",
#   "TimePoint_D65_vs_D100.tabular",
#   "TimePoint_D65_vs_D180.tabular",
#   "TimePoint_D65_vs_D280.tabular"
# )
# 
# # Initialize an empty list to store data frames
# dataframes_list <- list()
# 
# # Loop through each file, read its content, and add the 'comparison' column
# for (file_name in tabular_files) {
#   file_url <- paste0(raw_base_url, file_name)
#   
#   # Try reading the file
#   tryCatch({
#     # Read the tabular file
#     df <- read.table(file_url, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#     
#     # Add a 'comparison' column with the file name (excluding the extension)
#     colnames(df) <- c("gene_id", "mean_normalized_counts", "log2_fold_change", 
#                       "log2_fold_change_se", "wald_statistic", "p_value", "adj_p_value")
# 
#     df$comparison <- sub("\\.tabular$", "", file_name)
#     
#     # Store the data frame in the list
#     dataframes_list[[file_name]] <- df
#   }, error = function(e) {
#     message(paste("Error loading file:", file_name))
#   })
# }
# 
# DGE_data <- do.call(rbind, dataframes_list)
# DGE_data$DGE = ifelse(DGE_data$adj_p_value < 0.05 & abs(DGE_data$log2_fold_change) >= 1, 
#                       "TRUE", "FALSE")
# saveRDS(DGE_data, file = here("raw_data/agarwal_DGE_data.RDS"))


##### load DGE_data ####

DGE_data <- readRDS(here("raw_data/agarwal_DGE_data.RDS"))

DGE_genes <- DGE_data |> filter(DGE == TRUE) |> select(gene_id) |> unique()

short_read_data_DGE_counts <- short_read_data[DGE_genes$gene_id, ]

## get rownames 
# remove_zero_var_rows <- function(mat) {
#   mat[apply(mat, 1, function(x) min(x) != max(x)), ]
# }
# get_rownames <- function(method){
#     long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
#                                   method, "ROs", "gene_cpm.RDS")
#     
#     long_read_counts <- readRDS(long_read_counts)
#     
#     rownames(long_read_counts) <-
#       gsub("\\.\\d+$", "", rownames(long_read_counts))
#     
#     long_read_counts <- as.data.frame(long_read_counts)
#     
#     long_read_counts <- long_read_counts |>
#       mutate(rowname = rownames(long_read_counts)) |>
#       filter(str_detect(rowname, "^ENSG"))
#     
#     rownames(long_read_counts) <- long_read_counts$rowname
#     
#     DTU_DTE_DGE <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
#                              method, "ROs", "DGE_DTE_DTU.tsv"))
#     DGE_genes <- DTU_DTE_DGE |> filter(DGE == TRUE) |> select(gene_id) |> unique()
#     DGE_genes <- gsub("\\.\\d+$", "", DGE_genes$gene_id)
#     
#     intersect(rownames(long_read_counts), DGE_genes)
#     
# }
# 
# 
# bambu_genes <- get_rownames("bambu")
# Isoquant_genes <- get_rownames("Isoquant")
# 
# sh_bambu <- intersect(bambu_genes, rownames(short_read_data_DGE_counts))
# sh_Isoquant <- intersect(Isoquant_genes, rownames(short_read_data_DGE_counts))


method <- "bambu"
comparison <- "ROs"
load_data <- function(method) {
  long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                                method, comparison, "filtered_by_counts_and_biotype", "genes_cpm.RDS")
  
  long_read_counts <- readRDS(long_read_counts)
  
  # long_read_counts <- remove_zero_var_rows(long_read_counts)

  rownames(long_read_counts) <- gsub("\\.\\d+$", "", rownames(long_read_counts))
  
  long_read_counts <- as.data.frame(long_read_counts)
  
  long_read_counts <- long_read_counts |>
    mutate(rowname = rownames(long_read_counts)) |>
    filter(str_detect(rowname, "^ENSG"))
  
  rownames(long_read_counts) <- long_read_counts$rowname
  long_read_counts <- long_read_counts |>
    select(-rowname)
  #common rows with short reads matrix 
  
  DTU_DTE_DGE <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                    method, "ROs", "protein_coding" , "DGE_DTE_DTU.tsv"))
  DGE_genes <- DTU_DTE_DGE |> filter(DGE == TRUE) |> select(gene_id) |> unique()
  DGE_genes <- gsub("\\.\\d+$", "", DGE_genes$gene_id)
  
  long_read_DGE_genes <- intersect(rownames(long_read_counts), DGE_genes)
  
  common_genes <- intersect(long_read_DGE_genes, rownames(short_read_data_DGE_counts))

  long_read_counts <- long_read_counts[common_genes, ]
  short_read_counts <- short_read_data[common_genes, ]
  return(list(long_read_counts = long_read_counts, short_read_counts = short_read_counts))
}

bambu_res <- load_data("bambu")




method <- "bambu"

plot_scatter_plot <- function( method) { 

  output_plot_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, "/ROs/protein_coding/plots/short_read_correlation/")
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  long_read_counts <- load_data(method)$long_read_counts
  short_read_counts <-  load_data(method)$short_read_counts
  
  excluded_gene <- "ENSG00000210082"
  long_read_counts <- long_read_counts[!rownames(long_read_counts) %in% excluded_gene, ]
  short_read_counts <- short_read_counts[!rownames(short_read_counts) %in% excluded_gene, ]
  
  
  
  output_pdf <- file.path(output_plot_dir, "short_read_correlation_scatter_plots.pdf")
  pdf(output_pdf, width = 8, height = 6) 
  for (long_sample in colnames(long_read_counts)) {
    for (short_sample in colnames(short_read_counts)) {
      # Extract CPM values for the selected samples
      long_cpm <- long_read_counts[[long_sample]]
      short_cpm <- short_read_counts[[short_sample]]
      
      # Identify points with CPM > 1000 on both ends
      data_points <- data.frame(Long_Read = long_cpm, Short_Read = short_cpm)
     
      
      pearson_corr <- round(cor(long_cpm, short_cpm, method = "pearson"), 2)
      
      # Create scatter plot
      # Create scatter plot
      scatter_plot <- ggplot(data = data_points, 
                             aes(x = Long_Read, y = Short_Read)) +
        geom_point(size = 1, color = "blue") +
        geom_smooth(method = "lm", color = "red", linetype = "dashed") +
        # geom_text(aes(label = annotate), na.rm = TRUE, hjust = -0.2, vjust = -0.2, color = "darkgreen") +
        labs(
          title = paste("Scatter Plot of", long_sample, "vs", short_sample),
          x = paste("CPM - Long Read:", long_sample),
          y = paste("CPM - Short Read:", short_sample)
        ) +
        annotate("text", x = max(long_cpm), y = max(short_cpm), 
                 label = paste("Pearson r =", pearson_corr), 
                 hjust = 1.1, vjust = 1.1, color = "darkred", size = 4) +
        theme_minimal()
      
      # Print the plot to the PDF
      print(scatter_plot)
    }
  }
  
  dev.off()
  
}
  

method <- "bambu"
corr <- "spearman"
plot_correlation_heatmap <- function(method = "bambu", corr = "spearman") {
  
  long_read_counts <- load_data(method)$long_read_counts
  short_read_counts <-  load_data(method)$short_read_counts
  
  
  excluded_gene <- "ENSG00000210082"
  long_read_counts <- long_read_counts[!rownames(long_read_counts) %in% excluded_gene, ]
  short_read_counts <- short_read_counts[!rownames(short_read_counts) %in% excluded_gene, ]
  
  cor_matrix <- cor(short_read_counts, long_read_counts, method = corr)
  # Convert the correlation matrix to a long format for ggplot2
  cor_data <- melt(cor_matrix)
  cor_data$Var1 <- gsub(".tabular.txt", "", cor_data$Var1)
  ordered_levels <- c("D00.1", "D00.2", "D00.3",
                      "D10.3", "D10.4", "D10.5", "D10.6",
                      "D25.1", "D25.2", "D25.3", "D25.4",
                      "D65.2", "D65.3", "D65.5", "D65.6",
                      "D100.3", "D100.4", "D100.5", "D100.6",
                      "D180.1", "D180.2", "D180.3", "D180.4", "D180.5", "D180.6",
                      "D280.A2", "D280.B1", "D280.C1")
  
  # Convert cor_data$Var1 into a factor with the specified levels
  cor_data$Var1 <- factor(cor_data$Var1, levels = ordered_levels)
  
  # Verify the factor levels
  levels(cor_data$Var1)
  
  

  output_plot_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, "/ROs/protein_coding/plots/short_read_correlation/")
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  pdf(paste0(output_plot_dir, corr, "_short_read_correlation_heatmap.pdf"), width = 10, height = 6)
  p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3, angle = 45) + # Add annotations
    
    scale_fill_gradient2(low = "white", mid = "pink", high = "red", midpoint = 0.5, 
                         limits = c(0, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Short Read Samples", y = "Long Read Samples", 
         title = paste0(method, " ", corr," correlation between short read and long read counts"))
  print(p)
  dev.off()
}

plot_correlation_heatmap(method = "bambu", corr ="spearman")


