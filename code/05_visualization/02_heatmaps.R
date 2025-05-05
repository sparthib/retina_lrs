library(dplyr)
library(tibble)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

######## HEATMAP FUNCTIONS ########
plot_DTU_heatmap <- function(input_data_dir, quant_name, compare, tpm, groups, output_plots_dir, type) {
    tpm <- as.data.frame(tpm)
    if(compare == "FT_vs_RGC"){ 
      colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") }
    else if(compare == "ROs"){
      colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                         "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") }
    else if(compare == "RO_vs_RGC"){
      colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                         "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2" ) }
    tpm$isoform_id <- rownames(tpm)
    #remove version number 
    tpm$isoform_id <- gsub("\\..*", "", tpm$isoform_id)
    
    # Keep only tmm rows that are in DTU_isoforms based on its rownames
    significant_DTUs$isoform_id <-  ifelse(
      grepl("^ENST", significant_DTUs$isoform_id),  # Check if isoform_id starts with "ENST"
      gsub("\\..*", "", significant_DTUs$isoform_id),  # Remove everything after the first dot
      significant_DTUs$isoform_id  # Keep other isoform_id values unchanged
    )
    TPM_significant_isoforms <- tpm |>
      inner_join(significant_DTUs, by = "isoform_id") 
    
    # Create a new column with gene name and isoform name
    TPM_significant_isoforms <- TPM_significant_isoforms |> 
      mutate(gene_isoform = paste(gene_name, isoform_id, sep = "_")) |>
      tibble::column_to_rownames(var = "gene_isoform") |>
      dplyr::select(-isoform_id, -gene_name)

  # Convert to matrix and scale
  TPM_matrix <- as.matrix( TPM_significant_isoforms )
  TPM_matrix <- t(scale(t(TPM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  

  groups = groups
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen")
                            )) 
    }else if(compare == "RO_vs_RGC"){
      ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                              col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen",
                                                  "RGC" = "brown")
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


plot_DTE_heatmap <- function(input_data_dir, quant_name, compare, tpm, groups, output_plots_dir, type) {
  
  if(compare == "FT_vs_RGC"){ 
    colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") 
    } else if(compare == "ROs"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") 
    }else if(compare == "RO_vs_RGC"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2" ) }
  tpm <- as.data.frame(tpm)
  tpm$isoform_id <- rownames(tpm)
  #remove version number 
  tpm$isoform_id <- gsub("\\..*", "", tpm$isoform_id)
  
  # Keep only tmm rows that are in DTU_isoforms based on its rownames
  significant_DTEs$isoform_id <- ifelse(
    grepl("^ENST", significant_DTEs$isoform_id),  # Check if isoform_id starts with "ENST"
    gsub("\\..*", "", significant_DTEs$isoform_id),  # Remove everything after the first dot
    significant_DTEs$isoform_id  # Keep other isoform_id values unchanged
  )
  TPM_significant_isoforms <- tpm |>
    inner_join(significant_DTEs, by = "isoform_id")
  
  # Create a new column with gene name and isoform name
  TPM_significant_isoforms <- TPM_significant_isoforms |> 
    mutate(gene_isoform = paste(gene_name, isoform_id, sep = "_")) |>
    column_to_rownames(var = "gene_isoform") |>
    dplyr::select(-isoform_id, -gene_name)
  

  # Convert to matrix and scale
  TPM_matrix <- as.matrix( TPM_significant_isoforms )
  TPM_matrix <- t(scale(t(TPM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  
  groups = groups
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen")
                            ))
  }else if(compare == "RO_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen", 
                                        "RGC" = "brown")
                            ))
  }
  
  
  # Generate and draw heatmap
  pdf(file.path(heatmap_plots_dir, paste0(compare, "_DTE_heatmap.pdf")))
  ht_list <- Heatmap(TPM_matrix , name = paste0("scaled TPM expression of top", nrow(TPM_matrix), " DTE isoforms"),  row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
  
}

plot_DGE_heatmap <- function(input_data_dir, quant_name, compare, tpm, groups, output_plots_dir, type) {
  tpm <- as.data.frame(tpm)
  if(compare == "FT_vs_RGC"){ 
    colnames(tpm) <- c("H9_FT_1", "H9_FT_2", "H9_hRGC_1", "H9_hRGC_2") 
    } else if(compare == "ROs"){
    colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                       "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO") 
    } else if(compare == "RO_vs_RGC"){
      colnames(tpm) <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45" ,"EP1_WT_hRO_2" ,  "H9_BRN3B_hRO_2",
                         "H9_CRX_hRO_2" ,  "EP1_BRN3B_RO"  , "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2" ) }
  
  tpm$gene_id <- rownames(tpm)
  #remove version number 
  tpm$gene_id <- gsub("\\..*", "", tpm$gene_id)
  
  # Keep only tmm rows that are in DTU_isoforms based on its rownames
  significant_DGEs$gene_id <- ifelse(
    grepl("^ENST", significant_DGEs$gene_id),  # Check if isoform_id starts with "ENST"
    gsub("\\..*", "",  significant_DGEs$gene_id),  # Remove everything after the first dot
    significant_DGEs$gene_id  # Keep other isoform_id values unchanged
  )
  TPM_significant_genes <- tpm |>
    inner_join(significant_DGEs, by = "gene_id")
  
  # Create a new column with gene name and isoform name
  TPM_significant_genes <- TPM_significant_genes|> 
    mutate(gene_name_id = paste(gene_name, gene_id, sep = "_")) |>
    column_to_rownames(var = "gene_name_id") |>
    dplyr::select(-gene_id, -gene_name)
  
  # Convert to matrix and scale
  TPM_matrix <- as.matrix( TPM_significant_genes)
  TPM_matrix <- t(scale(t(TPM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  
  groups = groups
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen")
                            ))
  } else if(compare == "RO_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("Stage_3" = "purple", "Stage_1" = "orange", "Stage_2" = "seagreen",
                                                "RGC" = "brown"))) } 
  
  
  # Generate and draw heatmap
  pdf(file.path(heatmap_plots_dir, paste0(compare, "_DGE_heatmap.pdf")))
  ht_list <- Heatmap(TPM_matrix , name = paste0("scaled TPM expression of top", nrow(TPM_matrix), " DGE genes"),  row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "genes",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
  
}


counts_matrix_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/"
method <- "bambu"
compare <- "RO_vs_RGC"

load_and_plot_data <- function(method,compare, counts_matrix_dir) {
  isoform_tpm <- readRDS(file.path(counts_matrix_dir,
                                   method,
                                   compare, "filtered_by_counts_and_biotype",
                                   "filtered_isoform_cpm.RDS"))
  gene_tpm <- readRDS(file.path(counts_matrix_dir,
                                method,
                                compare, "filtered_by_counts_and_biotype",
                                "filtered_gene_cpm.RDS"))
  
  samples <- colnames(isoform_tpm)
  if (compare == "FT_vs_RGC"){
    groups <- c("FT", "FT", "RGC", "RGC")
   colnames(gene_tpm) <- samples
  }else if (compare == "ROs"){
    groups <- c("Stage_1", "Stage_1", "Stage_2","Stage_2", 
                "Stage_2", "Stage_3", "Stage_3")
  }else if (compare == "RO_vs_RGC"){
    groups <- c("Stage_1", "Stage_1", "Stage_2","Stage_2", 
                "Stage_2", "Stage_3", "Stage_3", "RGC", "RGC")
  }
    
    heatmap_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                   method, compare, "protein_coding","plots", "heatmaps")
    if (!dir.exists(heatmap_plots_dir)) {
      dir.create(heatmap_plots_dir, recursive = TRUE)
    }
      
    input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                             method, compare,"protein_coding" )
    
    DGE_DTE_DTU <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))
    
    significant_DTUs <- DGE_DTE_DTU |> dplyr::group_by(isoform_id) |>
      filter(DTU == TRUE) |>
      arrange(DTU_qval) |> 
      dplyr::select(isoform_id, gene_name) |>
      distinct() |> 
      head(n = 50)
    
    plot_DTU_heatmap(input_data_dir, method, compare, 
                     isoform_tpm, groups, heatmap_plots_dir)
    
    write_tsv(significant_DTUs, file.path(heatmap_plots_dir, "significant_DTUs.tsv"))
    
    significant_DTEs <- DGE_DTE_DTU |> dplyr::group_by(isoform_id) |>
      filter(DTE == TRUE) |>
      arrange(DTE_qval) |> 
      dplyr::select(isoform_id, gene_name) |>
      distinct() |> 
      head(n = 50)
    plot_DTE_heatmap(input_data_dir, method, compare, 
                     isoform_tpm, groups, heatmap_plots_dir)
    
    write_tsv(significant_DTEs, file.path(heatmap_plots_dir, "significant_DTEs.tsv"))
    
    significant_DGEs <- DGE_DTE_DTU |> dplyr::group_by(gene_id) |>
      filter(DGE == TRUE) |>
      arrange(DGE_qval) |> 
      dplyr::select(gene_name) |>
      distinct() |> 
      head(n = 50)
    
    plot_DGE_heatmap(input_data_dir, method, compare, 
                     gene_tpm, groups, heatmap_plots_dir)
    
    write_tsv(significant_DGEs, file.path(heatmap_plots_dir, "significant_DGEs.tsv"))
    
}

load_and_plot_data("bambu", "ROs", counts_matrix_dir)
load_and_plot_data("bambu", "FT_vs_RGC", counts_matrix_dir)
load_and_plot_data("bambu", "RO_vs_RGC", counts_matrix_dir)


