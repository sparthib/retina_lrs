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

######.   ##########

# Handle rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, var) != 0, , drop = FALSE]
}

analysis_type <- "RO_vs_RGC"
quant_method <- "bambu"

counts_matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                               quant_method, analysis_type, "filtered_by_counts_and_biotype")
splicing_factors_path <- "/users/sparthib/retina_lrs/raw_data/GeneCards-Pathway-Splicing.csv"

read_csv(splicing_factors_path) |> nrow()
# Load gene counts matrix function
load_gene_counts_matrix <- function(analysis_type, quant_method, splicing_factors_path,
                                    table_type = "DTE") {
  # Validate inputs
  if (!analysis_type %in% c("FT_vs_RGC", "ROs", "RO_vs_RGC")) {
    stop("Invalid analysis_type. Choose 'FT_vs_RGC', 'ROs' or 'RO_vs_RGC'.")
  }
  if (!quant_method %in% c("bambu", "Isoquant")) {
    stop("Invalid quant_method. Choose 'bambu' or 'isoquant'.")
  }
  
  counts_matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                                 quant_method, analysis_type, "filtered_by_counts_and_biotype")
  
  # Load splicing factors
  splicing_factors <- read_csv(splicing_factors_path) |>
    dplyr::select(ensembl_gene_id, `Gene Symbol`) |>
    rename(gene_id = ensembl_gene_id, gene_name = `Gene Symbol`) |>
    distinct(gene_id, .keep_all = TRUE)
  
  isoform_file <- file.path(counts_matrix_dir, "filtered_isoform_cpm.RDS")
  isoform_tpm <- readRDS(isoform_file)
  rownames(isoform_tpm) <- gsub("\\..*", "", rownames(isoform_tpm))
  
  input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", quant_method, analysis_type, "protein_coding")
  DGE_DTE_DTU <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))
  
  genes_and_isoforms <- DGE_DTE_DTU |> dplyr::select(c(gene_id, isoform_id)) |> distinct()
  genes_and_isoforms$gene_id <- gsub("\\..*", "", genes_and_isoforms$gene_id)
  genes_and_isoforms$isoform_id <- gsub("\\..*", "", genes_and_isoforms$isoform_id)
  nrow(genes_and_isoforms)
  
  genes_and_isoforms <- genes_and_isoforms |> filter(gene_id %in% splicing_factors$gene_id)
  nrow(genes_and_isoforms)
  
  isoform_tpm <- isoform_tpm[rownames(isoform_tpm) %in% genes_and_isoforms$isoform_id, ]
  

  # Filter significant genes
  if (table_type == "DTE") {
    significant_genes <- DGE_DTE_DTU |> filter(DTE == TRUE) |> dplyr::select(gene_id, isoform_id) |> distinct()
    significant_genes$gene_id <- gsub("\\..*", "", significant_genes$gene_id)
    significant_genes$isoform_id <- gsub("\\..*", "", significant_genes$isoform_id)
  } else if (table_type == "DTU") {
    significant_genes <- DGE_DTE_DTU |> filter(DTU == TRUE) |> dplyr::select(gene_id, isoform_id) |> distinct()
    significant_genes$gene_id <- gsub("\\..*", "", significant_genes$gene_id)
    significant_genes$isoform_id <- gsub("\\..*", "", significant_genes$isoform_id) 
  }
  
  isoform_tpm <- isoform_tpm[rownames(isoform_tpm) %in% significant_genes$isoform_id, ]
  isoform_tpm <- remove_zero_var_rows(isoform_tpm)
  
  # Define groups and output directory
  groups <- if (analysis_type == "ROs") {
    c("Stage_1", "Stage_1", "Stage_2", "Stage_2", "Stage_2", "Stage_3", "Stage_3")
  } else if (analysis_type == "FT_vs_RGC") {
    c("FT", "FT", "RGC", "RGC")
  } else if (analysis_type == "RO_vs_RGC") {
    c("Stage_1", "Stage_1", "Stage_2", "Stage_2", "Stage_2", "Stage_3", "Stage_3", "RGC", "RGC")
  }
  
  
  output_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu", 
                                quant_method, analysis_type, "protein_coding", "plots", "splicing_factor_analysis")
  dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(list(isoform_tpm = isoform_tpm,
              genes_and_isoforms = genes_and_isoforms,
              splicing_factors = splicing_factors, 
              groups = groups, 
              output_plots_dir = output_plots_dir))
}


# Function to plot heatmap
plot_heatmap <- function(tpm, genes_and_isoforms, groups, compare, output_plots_dir, splicing_factors_df, table_type) {
 
  if(compare == "FT_vs_RGC"){ 
    colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") }
  else if(compare == "ROs"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") }
  else if(compare == "RO_vs_RGC"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2" )
  } 
   # tpm <- tpm |>
   #  as.data.frame() |>
   #  rownames_to_column(var = "gene_id") |>
   #  inner_join(splicing_factors_df, by = "gene_id") |>
   #  mutate(gene = paste(gene_name, gene_id, sep = "_")) |>
   #  column_to_rownames(var = "gene") |>
   #  dplyr::select(-gene_id, -gene_name)
  
  tpm <- tpm |> 
    as.data.frame() |>
    rownames_to_column(var = "isoform_id") |>
    inner_join(genes_and_isoforms, by = "isoform_id") |> 
    inner_join(splicing_factors_df, by = "gene_id") |> 
    mutate(isoform = paste(gene_name, isoform_id, sep = "_")) |>
    column_to_rownames(var = "isoform") |>
    dplyr::select(-isoform_id, -gene_id, -gene_name)
  
  tpm_matrix <- tpm |> as.matrix() |> t() |> scale(center = TRUE, scale = TRUE) |> t()
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  ha <- switch(
    compare,
    "FT_vs_RGC" = HeatmapAnnotation(type = groups, annotation_name_side = "left",
                                    col = list(type = c("FT" = "skyblue", "RGC" = "brown")),
                                    annotation_name_gp = gpar(fontsize = 2)),
    "ROs" = HeatmapAnnotation(type = groups, annotation_name_side = "left",
                              col = list(type = c("Stage_3" = "purple", "Stage_2" = "orange", "Stage_1" = "seagreen")),
                              annotation_name_gp = gpar(fontsize = 2)),
    "RO_vs_RGC" = HeatmapAnnotation(type = groups, annotation_name_side = "left",
                                    col = list(type = c("Stage_3" = "purple", "Stage_2" = "orange", "Stage_1" = "seagreen", "RGC" = "brown")),
                                    annotation_name_gp = gpar(fontsize = 2))
  )
  
  pdf(file.path(output_plots_dir, paste0(compare, "_", table_type, "_splicing_factors_heatmap.pdf")))
  ht_list <- Heatmap(
    tpm_matrix, name = "Scaled TPM Expression of Splicing Factors", row_km = 5, col = col_fun,
    top_annotation = ha, show_row_names = TRUE, show_column_names = TRUE, 
    row_title = "Isoforms", row_names_gp = gpar(fontsize = 2),
    column_names_gp = gpar(fontsize = 5), show_row_dend = TRUE, show_column_dend = TRUE
  )
  draw(ht_list)
  dev.off()
  
  write.table(rownames(tpm), file = file.path(output_plots_dir, 
                                              paste0(compare, "_", table_type, "_splicing_factors_heatmap.tsv")), sep = "\t")
}



# Wrapper function to plot all heatmaps
plot_all_heatmaps <- function() {
  methods <- c("bambu")
  comparisons <- c("ROs", "FT_vs_RGC", "RO_vs_RGC")
  
  for (method in methods) {
    for (comparison in comparisons) {
      DTE_data <- load_gene_counts_matrix(comparison, method,  splicing_factors_path, "DTE")
      DTU_data <- load_gene_counts_matrix(comparison, method, splicing_factors_path, "DTU")
      # DGE_data <- load_gene_counts_matrix(comparison, method,splicing_factors_path, "DGE")
      
      plot_heatmap(DTE_data$isoform_tpm, DTE_data$genes_and_isoforms, DTE_data$groups, comparison, DTE_data$output_plots_dir, 
                   DTE_data$splicing_factors, "DTE")
      plot_heatmap(DTU_data$isoform_tpm, DTU_data$genes_and_isoforms, DTU_data$groups, comparison, DTU_data$output_plots_dir, 
                   DTU_data$splicing_factors, "DTU")
      # plot_heatmap(DGE_data$isoform_tpm, DGE_data$groups, comparison, DGE_data$output_plots_dir, 
      #              DGE_data$splicing_factors, "DGE")
    }
  }
}

plot_all_heatmaps()

