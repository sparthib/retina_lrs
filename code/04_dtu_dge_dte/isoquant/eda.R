library("ggplot2")
library("DGEobj.utils")
library("patchwork")
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library("ggVennDiagram")
library(readxl)
library(stringr)


# Define directories and common variables
plots_dir <- "/users/sparthib/retina_lrs/plots/de/isoquant"

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
             "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
             "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")

# Function to load data
load_expression_data <- function(data_type) {
  counts <- read.table(file.path(isoquant_dir, paste0("OUT.", data_type, "_grouped_counts.tsv")), header = TRUE)
  tpm <- read.table(file.path(isoquant_dir, paste0("OUT.", data_type, "_grouped_tpm.tsv")), header = TRUE)
  colnames(counts) <- c("id", samples)
  colnames(tpm) <- c("id", samples)
  list(counts = counts, tpm = tpm)
}

# Function to filter rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}

# Function to prepare TPM data for PCA
prepare_tpm <- function(tpm, sample_indices) {
  selected_tpm <- tpm[, sample_indices]
  rownames(selected_tpm) <- tpm[, 1]
  remove_zero_var_rows(selected_tpm)
}

# PCA and plotting function
pca_plots_dir <- file.path(plots_dir, "pca")
plot_pca <- function(tpm, samples, groups, output_name, output_dir) {
  pc <- prcomp(t(log2(tpm + 1)), scale = TRUE)
  pcr <- data.frame(pc$x[, 1:2], row.names = samples)  # PC scores
  
  # PCA plot
  p <- ggplot(pcr, aes(PC1, PC2, color = groups)) +  
    geom_point(size = 2) +
    theme_bw() +
    ggtitle(paste("PCA on", output_name, "Expression")) +
    xlab(paste("PC1 (", round(pc$sdev[1]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +
    ylab(paste("PC2 (", round(pc$sdev[2]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +
    geom_label(aes(label = rownames(pcr)), size = 2, fill = "white", alpha = 0.7)
  
  # Scree plot
  var_explained <- pc$sdev[1:10]^2 / sum(pc$sdev^2)
  var_explained_df <- data.frame(
    Principal_Component = factor(1:length(var_explained)),
    Variance_Explained = var_explained
  )
  
  q <- ggplot(var_explained_df, aes(x = Principal_Component, y = Variance_Explained)) + 
    geom_bar(stat = "identity", fill = "skyblue") +
    xlab("Principal Component") +  
    ylab("Variance Explained") +  
    ggtitle("Scree Plot") +
    ylim(0, 1) +
    theme_minimal()
  
  # Save plots
  pdf(file.path(output_dir, paste0(output_name, "_level_pca.pdf")), width = 12, height = 6)
  print(p + q)
  dev.off()
}

# Prepare data and plot for isoforms
isoform_data <- load_expression_data("transcript")
RO_isoform_tpm <- prepare_tpm(isoform_data$tpm, 2:8)
#remove isoform id version 
RO_isoform_tpm <- RO_isoform_tpm |> rownames_to_column(var = "isoform_id") |>
  mutate(isoform_id = gsub("\\..*", "", isoform_id)) |> column_to_rownames(var = "isoform_id")

RO_isoquant_groups <- factor(c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"))
pca_plots_dir <- file.path(plots_dir, "pca")
plot_pca(RO_isoform_tpm, samples[2:8], RO_isoquant_groups, "isoform", pca_plots_dir)

# Prepare data and plot for genes
gene_data <- load_expression_data("gene")
RO_gene_tpm <- prepare_tpm(gene_data$tpm, 2:8)
plot_pca(RO_gene_tpm, samples[2:8], RO_isoquant_groups, "gene", pca_plots_dir)


######## HEATMAPS ########

heatmap_plots_dir <- file.path(plots_dir, "heatmaps")
if (!dir.exists(heatmap_plots_dir)) {
  dir.create(heatmap_plots_dir, recursive = TRUE)
}

input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs"



plot_heatmap <- function(input_data_dir, quant_name, compare, tmm, groups, output_plots_dir) {
  # Read the DEXSeqSwitchList file
    switch_file <- read_tsv(file.path(input_data_dir, "DGE_DTU_DTE.tsv"))
    significant_DTUs <- switch_file |> dplyr::group_by(isoform_id ) |>
      filter(
        DTU_qval == min(DTU_qval) & 
          abs(DTU_dIF) == max(abs(DTU_dIF))
      ) |> filter(DTU_qval < 0.05 & abs(DTU_dIF) > 0.1) |>
      arrange(DTU_qval) |> 
      select(isoform_id, gene_name) |>
      distinct() |> 
      head(n = 200)
  
  # Keep only tmm rows that are in DTU_isoforms based on its rownames
  tpm <- RO_isoform_tpm
  tpm <- as.data.frame(tpm)
  tpm$isoform_id <- rownames(tpm)
  TPM_significant_isoforms <- tpm |>
    inner_join(significant_DTUs, by = "isoform_id")
  
  # Create a new column with gene name and isoform name
  TPM_significant_isoforms <- TPM_significant_isoforms |> 
    mutate(gene_isoform = paste(gene_name, isoform_id, sep = "_")) |>
    column_to_rownames(var = "gene_isoform") |>
    select(-isoform_id, -gene_name)
  
  
  # Convert to matrix and scale
  TPM_matrix <- as.matrix( TPM_significant_isoforms )
  TPM_matrix <- t(scale(t(TPM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Heatmap annotation
  compare = "ROs"
  groups = RO_isoquant_groups
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen")
                            ))
  }

  
  # Generate and draw heatmap
  pdf(file.path(heatmap_plots_dir, paste0(compare, "_DTU_heatmap.pdf")))
  ht_list <- Heatmap(TPM_matrix , name = paste0("scaled TPM expression of top", nrow(TPM_matrix), " DTU isoforms"),  row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
  
}

plot_heatmap(input_data_dir, "isoquant", "ROs", ROs_isoquant_tpm, ROs_isoquant_groups, heatmap_plots_dir)


####### Upset Plot ########

# save venn diagram as pdf 
DGE_DTU_DTE <- read_tsv(file.path(input_data_dir, "DGE_DTU_DTE.tsv"))

gene_overlaps = DGE_DTU_DTE |> dplyr::select( gene_id, condition_1, condition_2,DGE, DTU )

gene_overlaps <- gene_overlaps |>
  mutate(condition = case_when(
    condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45" ~ "RO_D100_vs_RO_D45",
    condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45" ~ "RO_D200_vs_RO_D45",
    condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100" ~ "RO_D100_vs_RO_D200",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))



gene_overlaps <- gene_overlaps |> mutate(DTU_RO_D100_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D200_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D100_vs_RO_D200 = ifelse(DTU == TRUE & condition == "RO_D100_vs_RO_D200", TRUE, FALSE),
                                         DGE_RO_D100_vs_RO_D200 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D200", TRUE, FALSE),
                                         DGE_RO_D200_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DGE_RO_D100_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE))

gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_RO_D100_vs_RO_D45 = any(DTU_RO_D100_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D45 = any(DTU_RO_D200_vs_RO_D45), 
            DTU_RO_D100_vs_RO_D200 = any(DTU_RO_D100_vs_RO_D200), 
            DGE_RO_D100_vs_RO_D200 = any(DGE_RO_D100_vs_RO_D200), 
            DGE_RO_D200_vs_RO_D45 = any(DGE_RO_D200_vs_RO_D45), 
            DGE_RO_D100_vs_RO_D45 = any(DGE_RO_D100_vs_RO_D45)) 

rownames(gene_overlaps) <- gene_overlaps$gene_id




# Prepare data for VennDiagram
venn_data <- list(
  DTU_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D45],
  DTU_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D45],
  DTU_RO_D100_vs_RO_D200 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D200],
  DGE_RO_D100_vs_RO_D200 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D200],
  DGE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D45],
  DGE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D45]
)


library("UpSetR")
upset_path <- file.path(output_plots_dir, "DGE_DTU_ROs_upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()



plot_barplot <- function(df, condition){ 
  
  dge_overlaps = df |> group_by(gene_id) |> dplyr::select( -DTE) |> 
    summarise(DGE=any(DGE), DTU=any(DTU))  
  
  dte_overlaps = df |> group_by(gene_id) |> dplyr::select(-DGE) |> 
    summarise(DTE = any(DTE), DTU=any(DTU)) 
  
  dge_contingency_table <- as_tibble(xtabs(~ DGE + DTU, data = dge_overlaps))
  dge_contingency_table <- plyr::ddply(dge_contingency_table, ~DGE, transform, Prop = n / sum(n))
  
  
  dge_bar <- ggplot(dge_contingency_table, aes(x = DGE, y = n, fill = DTU)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "DGE and DTU genes", x = "DGE", y = "Count") + 
    scale_y_continuous(labels = scales::comma) +  # Display y-axis as percentage
    theme_minimal() +
    geom_text(aes(label = n), 
              position = position_stack(vjust = 0.5), 
              size = 2)
  
  dte_contingency_table <- xtabs(~ DTE + DTU, data = dte_overlaps)
  dte_contingency_table <- plyr::ddply(as_tibble(dte_contingency_table), ~DTE, transform, Prop = n / sum(n))
  
  dte_bar <- ggplot(dte_contingency_table, aes(x = DTE, y = n, fill = DTU)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "DTE and DTU genes", x = "DTE", y = "Count") + 
    scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
    theme_minimal() +
    geom_text(aes(label = n), 
              position = position_stack(vjust = 0.5), 
              size = 2)
  
  library(patchwork)
  p <- dge_bar + dte_bar + plot_annotation(title = condition)
  ggsave(path = output_plots_dir, 
         device = "pdf", 
         plot = p, 
         filename = paste0(condition, "df_barplot.pdf"))
  
}

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "B_RO_D100" &
                                     condition_2 == "C_RO_D45"), "RO_D100_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
                                     condition_2 == "C_RO_D45"), "RO_D200_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
                                     condition_2 == "B_RO_D100"), "RO_D100_vs_RO_D200")



######## RetNet Plots #######

RetNet_gene_list <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")

genes_and_diseases <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")


intersect(significant_DTUs$gene_name, RetNet_gene_list$Symbol)

# Initialize an empty list to store results
results_list <- list()
for (gene in intersect(significant_DTUs$gene_name, RetNet_gene_list$Symbol)) {
  # Filter the data
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  
  # Append results to the list
  results_list[[gene]] <- result$disease_category
}

# Convert the list to a dataframe
results_df <- data.frame(
  Gene = names(results_list),
  Disease_Category = I(results_list),
  stringsAsFactors = FALSE
)

# If you want to view the dataframe, you can use:



######## GO Analysis #########















