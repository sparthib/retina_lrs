library(dplyr)
library(ggplot2)
library(here)
library(stringr)

url <- "https://raw.githubusercontent.com/WahlinLab/Organoid_RNAseq_SciData22/refs/heads/main/Normalized_counts_D0-D280.tabular"
short_read_data <- read.table(url, header = TRUE, sep = "\t", row.names = 1)

# Extract the day numbers from the column names
day_order <- c(0, 10, 25, 65, 100, 180, 280)
day_extracted <- as.numeric(stringr::str_extract(colnames(short_read_data), "\\d+"))
column_order <- order(factor(day_extracted, levels = day_order))
short_read_data <- short_read_data[, column_order]
saveRDS(short_read_data,
        file = here("raw_data/agarwal_short_norm_counts.RDS"))

## get rownames 
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}
get_rownames <- function(method){
    long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                                  method, "ROs", "gene_cpm.RDS")
    
    long_read_counts <- readRDS(long_read_counts)
    
    long_read_counts <- remove_zero_var_rows(long_read_counts)
    
    rownames(long_read_counts) <-
      gsub("\\.\\d+$", "", rownames(long_read_counts))
    
    long_read_counts <- as.data.frame(long_read_counts)
    
    long_read_counts <- long_read_counts |>
      mutate(rowname = rownames(long_read_counts)) |>
      filter(str_detect(rowname, "^ENSG"))
    
    rownames(long_read_counts) <- long_read_counts$rowname
    rownames(long_read_counts)
}

bambu_genes <- get_rownames("bambu")
Isoquant_genes <- get_rownames("Isoquant")
length(Isoquant_genes)
length(bambu_genes)

# > length(Isoquant_genes)
# [1] 39564
# > length(bambu_genes)
# [1] 59218
common_genes <- intersect(bambu_genes, Isoquant_genes)
length(common_genes)

sh_bambu <- intersect(bambu_genes, rownames(short_read_data))
sh_Isoquant <- intersect(Isoquant_genes, rownames(short_read_data))
length(sh_bambu)
length(sh_Isoquant)
# > length(sh_bambu)
# [1] 53616
# > length(sh_Isoquant)
# [1] 35568



load_data <- function(method) {
  long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                                method, "ROs", "gene_cpm.RDS")
  
  long_read_counts <- readRDS(long_read_counts)
  
  # long_read_counts <- remove_zero_var_rows(long_read_counts)

  rownames(long_read_counts) <-
    gsub("\\.\\d+$", "", rownames(long_read_counts))
  
  long_read_counts <- as.data.frame(long_read_counts)
  
  long_read_counts <- long_read_counts |>
    mutate(rowname = rownames(long_read_counts)) |>
    filter(str_detect(rowname, "^ENSG"))
  
  rownames(long_read_counts) <- long_read_counts$rowname
  long_read_counts <- long_read_counts |>
    select(-rowname)
  #common rows with short reads matrix 
  common_genes <- intersect(rownames(long_read_counts), rownames(short_read_data))
  # length(common_genes) 56612 Isoquant
  # length(common_genes)  56614 Bambu
  
  #filter out the common genes
  long_read_counts <- long_read_counts[common_genes, ]
  short_read_counts <- short_read_data[common_genes, ]
  return(list(long_read_counts = long_read_counts, short_read_counts = short_read_counts))
}

######
# long_read_counts <- load_data("bambu")$long_read_counts
# short_read_counts <-  load_data("bambu")$short_read_counts
# nrow(long_read_counts)
# nrow(short_read_counts)
# bambu_pearson_corr <- cor(short_read_counts, long_read_counts, method = "pearson") |> melt()
# bambu_spearman_corr <- cor(short_read_counts, long_read_counts, method = "spearman") |> melt()
# 
# long_read_counts <- load_data("Isoquant")$long_read_counts
# short_read_counts <-  load_data("Isoquant")$short_read_counts
# nrow(long_read_counts)
# nrow(short_read_counts)
# Isoquant_pearson_corr <- cor(short_read_counts, long_read_counts, method = "pearson") |> melt()
# Isoquant_spearman_corr <- cor(short_read_counts, long_read_counts, method = "spearman") |> melt()
######


###### get_correlation_data ######
# get_correlation_data <- function(quant, corr) {
#   long_read_counts <- load_data(quant)$long_read_counts
#   short_read_counts <-  load_data(quant)$short_read_counts
# 
#   cor_matrix <- cor(short_read_counts, long_read_counts, method = corr)
#   # Convert the correlation matrix to a long format for ggplot2
#   cor_data <- melt(cor_matrix)
# 
#   return(cor_data)
# }
# 
# bambu_pearson <- get_correlation_data("bambu", "pearson")
# bambu_spearman <- get_correlation_data("bambu", "spearman")
# 
# isoquant_pearson <- get_correlation_data("Isoquant", "pearson")
# isoquant_spearman <- get_correlation_data("Isoquant", "spearman")
# 
# range(bambu_pearson$value)
# range(bambu_spearman$value)
# range(isoquant_pearson$value)
# range(isoquant_spearman$value)


#####

##### make correlation plots ##### 

library(ggplot2)
library(reshape2)
plot_correlation_heatmap <- function(method = "bambu", corr = "pearson") {
  
  long_read_counts <- load_data(method)$long_read_counts
  short_read_counts <-  load_data(method)$short_read_counts
  
  cor_matrix <- cor(short_read_counts, long_read_counts, method = corr)
  # Convert the correlation matrix to a long format for ggplot2
  cor_data <- melt(cor_matrix)

  output_plot_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, "/ROs/plots/short_read_correlation/")
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  pdf(paste0(output_plot_dir, corr, "_short_read_correlation_heatmap.pdf"), width = 10, height = 10)
  p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 2) + # Add annotations
    
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, 
                         limits = c(0, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Short Read Samples", y = "Long Read Samples", 
         title = paste0(method, " ", corr," correlation between short read and long read counts"))
  print(p)
  dev.off()
}

plot_correlation_heatmap(method = "bambu", corr = "pearson")
plot_correlation_heatmap(method = "bambu", corr ="spearman")

plot_correlation_heatmap(method = "Isoquant", corr ="pearson")
plot_correlation_heatmap(method = "Isoquant", corr ="spearman")

