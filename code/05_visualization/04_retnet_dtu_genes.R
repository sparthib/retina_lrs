library(readr)
library(tidyr)
library(dplyr)
library(readxl)
library(biomaRt)
library(stringr)
library(tibble)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)


data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")

remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}


raw_data_dir <- file.path(code_dir, "raw_data")
RetNet_gene_list <- read_excel(file.path(raw_data_dir, "RetNet.xlsx"),
                               sheet = "genes_and_locations")

genes_and_diseases <- read_excel(file.path(raw_data_dir, "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")

#convert gene name to gene ID and add to the data-frame using bio-maRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "hgnc_symbol",
               values = RetNet_gene_list$Symbol,
               mart = mart)

#from genes and diseases
gene_vector <- readRDS(file = file.path(code_dir, "processed_data/dtu/retnet_gene_vector.RDS"))
disease_genes <- tibble(gene = gene_vector) |> distinct()
# disease_gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                  filters = "hgnc_symbol",
#                  values = disease_genes$gene,
#                  mart = mart)
# # nrow(disease_genes)
# #get rows with duplicate gene ids in gene_id
# duplicated_gene_id <- disease_gene_id[duplicated(disease_gene_id$hgnc_symbol),]

results_list <- list()
for (gene in disease_genes$gene) {
  # Filter the data
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  
  # Append results to the list 
  results_list[[gene]] <- result$disease_category
}

length(results_list)
#290

#convert results_list to a dataframe
results_df <- data.frame(
  gene_name = names(results_list),
  disease_category = as.character(results_list),
  stringsAsFactors = FALSE
)

write.table(results_df, file = file.path(code_dir, "processed_data/dtu/retnet_disease_genes.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

analysis_type <- "FT_vs_RGC"
quant_method <- "bambu"
table_type <- "DTE"

load_gene_counts_matrix <- function(analysis_type, quant_method,
                                    table_type = "DTE") {
  
  # Validate inputs
  if (!analysis_type %in% c("FT_vs_RGC", "ROs", "RO_vs_RGC")) {
    stop("Invalid analysis_type. Choose 'FT_vs_RGC', 'ROs' or 'RO_vs_RGC'.")
  }
  if (!quant_method %in% c("bambu", "Isoquant")) {
    stop("Invalid quant_method. Choose 'bambu' or 'isoquant'.")
  }
  
  counts_matrix_dir <- file.path(data_dir, "06_quantification/counts_matrices/",
                                 quant_method, analysis_type, "filtered_by_counts_and_biotype")
  
  isoform_file <- file.path(counts_matrix_dir, "filtered_isoform_cpm.RDS")
  isoform_tpm <- readRDS(isoform_file)
  rownames(isoform_tpm) <- gsub("\\..*", "", rownames(isoform_tpm))
  
  input_data_dir <- file.path(code_dir, "processed_data/dtu/", quant_method, analysis_type, "protein_coding")
  DGE_DTE_DTU <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))
  
  genes_and_isoforms <- DGE_DTE_DTU |> dplyr::select(c(gene_id, isoform_id, gene_name)) |> distinct()
  genes_and_isoforms$gene_id <- gsub("\\..*", "", genes_and_isoforms$gene_id)
  genes_and_isoforms$isoform_id <- gsub("\\..*", "", genes_and_isoforms$isoform_id)
  nrow(genes_and_isoforms)
  
  genes_and_isoforms <- genes_and_isoforms |> filter(gene_name %in% results_df$gene_name)
  nrow(genes_and_isoforms)
  
  isoform_tpm <- isoform_tpm[rownames(isoform_tpm) %in% genes_and_isoforms$isoform_id, ]
  
  
  # Filter significant genes
  if (table_type == "DTE") {
    
    significant_genes <- DGE_DTE_DTU |> 
      dplyr::filter(DTE == TRUE) |>
      distinct()
    
    # Keep only genes with multiple isoforms
    significant_genes <- significant_genes |>
      group_by(gene_id) |>
      filter(n() > 1) |>
      ungroup() |>
      distinct()
    
    # only keep genes in genes_and_isoforms
    significant_genes <- significant_genes |>
      filter(gene_id %in% genes_and_isoforms$gene_id) |>
      distinct()
    
    # Keep only isoforms from those top 30 genes
    top_isoforms <- significant_genes |>
      arrange(DTE_qval, desc(abs(DTE_log2FC))) |>
      slice_head(n = 30) |>
      pull(isoform_id)
    
    print("length of top DTE isoforms")
    print(length(top_isoforms))
    
    significant_genes <- significant_genes |>
      filter(isoform_id %in% top_isoforms) |>
      dplyr::select(gene_id, isoform_id)
    

  } else if (table_type == "DTU") {
    
    significant_genes <- DGE_DTE_DTU |> 
      dplyr::filter(DTU == TRUE) |>
      distinct()
    
    significant_genes <- significant_genes |>
      group_by(gene_id) |>
      filter(n() > 1) |>
      ungroup() |>
      distinct()
    
    # only keep genes in genes_and_isoforms
    significant_genes <- significant_genes |>
      filter(gene_id %in% genes_and_isoforms$gene_id) |>
      distinct()
    
    top_isoforms <- significant_genes |>
      arrange(DTU_qval, desc(abs(dIF))) |>
      slice_head(n = 30) |>
      pull(isoform_id)
    
    print("length of top DTU isoforms")
    print(length(top_isoforms))
    
    significant_genes <- significant_genes |>
      filter(isoform_id %in% top_isoforms) |>
      dplyr::select(gene_id, isoform_id)
    
    print

  }
  
  isoform_tpm <- isoform_tpm[rownames(isoform_tpm) %in% significant_genes$isoform_id, ]

  print("length of isoform_tpm")
  print(length(rownames(isoform_tpm)))
  # Define groups and output directory
  groups <- if (analysis_type == "ROs") {
    c("Stage_1", "Stage_1", "Stage_2", "Stage_2", "Stage_2", "Stage_3", "Stage_3")
  } else if (analysis_type == "FT_vs_RGC") {
    c("FT", "FT", "RGC", "RGC")
  } else if (analysis_type == "RO_vs_RGC") {
    c("Stage_1", "Stage_1", "Stage_2", "Stage_2", "Stage_2", "Stage_3", "Stage_3", "RGC", "RGC")
  }
  
  output_plots_dir <- file.path(code_dir, "processed_data/dtu", 
                                quant_method, analysis_type, "protein_coding", "plots", "retnet")
  dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(list(isoform_tpm = isoform_tpm,
              genes_and_isoforms = genes_and_isoforms,
              groups = groups, 
              output_plots_dir = output_plots_dir))
}



plot_heatmap <- function(tpm, genes_and_isoforms, groups, compare, output_plots_dir, table_type) {
  
  if(compare == "FT_vs_RGC"){ 
    colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") }
  else if(compare == "ROs"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") }
  else if(compare == "RO_vs_RGC"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2" )
  } 
  
  tpm <- tpm |> 
    as.data.frame() |>
    rownames_to_column(var = "isoform_id") |>
    inner_join(genes_and_isoforms, by = "isoform_id") |> 
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
  
  pdf(file.path(output_plots_dir, paste0(compare, "_", table_type, "_retnet_heatmap_top_30.pdf")))
  ht_list <- Heatmap(
    tpm_matrix, name = paste0("Top 30 Isoform TPM Expression of ", table_type, " RetNet Genes"), row_km = 5, col = col_fun,
    top_annotation = ha, show_row_names = TRUE, show_column_names = TRUE, 
    row_title = "Isoforms", row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 5), show_row_dend = TRUE, show_column_dend = TRUE
  )
  draw(ht_list)
  dev.off()
  
  write.table(rownames(tpm), file = file.path(output_plots_dir, paste0(compare, "_", table_type, "_retnet_heatmap_top_30.tsv")), sep = "\t")
}


method <- "bambu"
comparison <- "FT_vs_RGC"

DTE_data <- load_gene_counts_matrix(comparison, method,  "DTE")
DTU_data <- load_gene_counts_matrix(comparison, method, "DTU")
# DGE_data <- load_gene_counts_matrix(comparison, method,splicing_factors_path, "DGE")

plot_heatmap(DTE_data$isoform_tpm, DTE_data$genes_and_isoforms, DTE_data$groups, 
             comparison, DTE_data$output_plots_dir,  "DTE")
plot_heatmap(DTU_data$isoform_tpm, DTU_data$genes_and_isoforms, DTU_data$groups,
             comparison, DTU_data$output_plots_dir,  "DTU")




