library(readr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
# install.packages("DGEobj.utils")
library("DGEobj.utils")
library(grid)


###### ADD ENSEMBL ID TO THE DATA AND SAVE ######
# head(splicing_factors)
# 
# #convert gene symbol to ENSEMBL ID
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ids <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                filters = "hgnc_symbol", values = splicing_factors$`Gene Symbol`,
#                mart = ensembl)
# 
# 
# # Merge the data
# splicing_factors <- merge(splicing_factors, gene_ids, by.x = "Gene Symbol", by.y = "hgnc_symbol")
# 
# # Save the data
# write_csv(splicing_factors, here("raw_data", "GeneCards-Pathway-Splicing.csv"))


# Load required helper functions
source("/users/sparthib/retina_lrs/code/05_visualization//helper.R")

# Handle rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, var) != 0, , drop = FALSE]
}

# Load gene counts matrix function
load_gene_counts_matrix <- function(analysis_type, quant_method, counts_matrix_dir, splicing_factors_path,
                                    table_type = "DTE") {
  # Validate inputs
  if (!analysis_type %in% c("FT_vs_RGC", "ROs")) {
    stop("Invalid analysis_type. Choose 'FT_vs_RGC' or 'ROs'.")
  }
  if (!quant_method %in% c("bambu", "Isoquant")) {
    stop("Invalid quant_method. Choose 'bambu' or 'isoquant'.")
  }
  
  # Load splicing factors
  splicing_factors <- read_csv(splicing_factors_path) |>
    dplyr::select(ensembl_gene_id, `Gene Symbol`) |>
    rename(gene_id = ensembl_gene_id, gene_name = `Gene Symbol`) |>
    distinct(gene_id, .keep_all = TRUE)
  
  gene_file <- file.path(counts_matrix_dir, quant_method, analysis_type, "filtered", "gene_cpm.RDS")
  
  # Load and process gene TPM
  gene_tpm <- readRDS(gene_file)
  rownames(gene_tpm) <- gsub("\\..*", "", rownames(gene_tpm))
  gene_tpm <- gene_tpm[rownames(gene_tpm) %in% splicing_factors$gene_id, ]
  
  input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", quant_method, analysis_type)
  DGE_DTE_DTU <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))
  
  # Filter significant genes
  if (table_type == "DTE") {
    significant_genes <- DGE_DTE_DTU$gene_id[DGE_DTE_DTU$DTE]
    significant_genes <- gsub("\\..*", "", significant_genes)
    significant_genes <- unique(significant_genes)
  } else if (table_type == "DTU") {
    significant_genes <- DGE_DTE_DTU$gene_id[DGE_DTE_DTU$DTU]
    significant_genes <- gsub("\\..*", "", significant_genes)
    significant_genes <- unique(significant_genes)
  } else if (table_type == "DGE") {
    significant_genes <- DGE_DTE_DTU$gene_id[DGE_DTE_DTU$DGE]
    significant_genes <- gsub("\\..*", "", significant_genes)
    significant_genes <- unique(significant_genes)
  }
  
  gene_tpm <- gene_tpm[rownames(gene_tpm) %in% significant_genes, ]
  gene_tpm <- remove_zero_var_rows(gene_tpm)
  
  # Define groups and output directory
  groups <- if (analysis_type == "ROs") {
    c("RO_D45", "RO_D45", "RO_D100", "RO_D100", "RO_D100", "RO_D200", "RO_D200")
  } else {
    c("FT", "FT", "RGC", "RGC")
  }
  
  output_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu", 
                                quant_method, analysis_type, "plots", "splicing_factor_analysis")
  dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(list( gene_tpm = gene_tpm, 
              splicing_factors = splicing_factors, groups = groups, 
              output_plots_dir = output_plots_dir))
}

# Function to plot heatmap
plot_heatmap <- function(tpm, groups, compare, output_plots_dir, splicing_factors_df, table_type) {
 
  if(compare == "FT_vs_RGC"){ 
    colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") }
  else if(compare == "ROs"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") }
  
   tpm <- tpm |>
    as.data.frame() |>
    rownames_to_column(var = "gene_id") |>
    inner_join(splicing_factors_df, by = "gene_id") |>
    mutate(gene = paste(gene_name, gene_id, sep = "_")) |>
    column_to_rownames(var = "gene") |>
    dplyr::select(-gene_id, -gene_name)
  
  tpm_matrix <- tpm |> as.matrix() |> t() |> scale(center = TRUE, scale = TRUE) |> t()
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  ha <- switch(
    compare,
    "FT_vs_RGC" = HeatmapAnnotation(type = groups, annotation_name_side = "left",
                                    col = list(type = c("FT" = "lightgreen", "RGC" = "brown")),
                                    annotation_name_gp = gpar(fontsize = 2)),
    "ROs" = HeatmapAnnotation(type = groups, annotation_name_side = "left",
                              col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen")),
                              annotation_name_gp = gpar(fontsize = 2))
  )
  
  pdf(file.path(output_plots_dir, paste0(compare, "_", table_type, "_splicing_factors_heatmap.pdf")))
  ht_list <- Heatmap(
    tpm_matrix, name = "Scaled TPM Expression of Splicing Factors", row_km = 5, col = col_fun,
    top_annotation = ha, show_row_names = TRUE, show_column_names = TRUE, 
    row_title = "Isoforms", row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 5), show_row_dend = TRUE, show_column_dend = TRUE
  )
  draw(ht_list)
  dev.off()
  
  write.table(rownames(tpm), file = file.path(output_plots_dir, paste0(compare, "_", table_type, "_splicing_factors_heatmap.tsv")), sep = "\t")
}

counts_matrix_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/"
splicing_factors_path <- "/users/sparthib/retina_lrs/raw_data/GeneCards-Pathway-Splicing.csv"

# Wrapper function to plot all heatmaps
plot_all_heatmaps <- function() {
  methods <- c("bambu", "Isoquant")
  comparisons <- c("ROs", "FT_vs_RGC")
  method <- "Isoquant"
  comparison <- "ROs"
  
  for (method in methods) {
    for (comparison in comparisons) {
      DTE_data <- load_gene_counts_matrix(comparison, method, counts_matrix_dir, splicing_factors_path, "DTE")
      DTU_data <- load_gene_counts_matrix(comparison, method, counts_matrix_dir, splicing_factors_path, "DTU")
      DGE_data <- load_gene_counts_matrix(comparison, method, counts_matrix_dir, splicing_factors_path, "DGE")
      
      plot_heatmap(DTE_data$gene_tpm, DTE_data$groups, comparison, DTE_data$output_plots_dir, 
                   DTE_data$splicing_factors, "DTE")
      plot_heatmap(DTU_data$gene_tpm, DTU_data$groups, comparison, DTU_data$output_plots_dir, 
                   DTU_data$splicing_factors, "DTU")
      plot_heatmap(DGE_data$gene_tpm, DGE_data$groups, comparison, DGE_data$output_plots_dir, 
                   DGE_data$splicing_factors, "DGE")
    }
  }
}

plot_all_heatmaps()

