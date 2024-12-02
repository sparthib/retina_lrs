#####clustering based on only top isoforms that show DTU. 
library(readr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library("DGEobj.utils")
library(grid)

FT_vs_RGC_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/counts_gene.txt",
                                  sep = "\t", row.names = 1, header = TRUE)
ROs_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/counts_gene.txt",
                            sep = "\t", row.names = 1, header = TRUE)
FT_vs_RGC_isoquant_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/OUT.gene_grouped_tpm.tsv",
                                     sep = "\t", row.names = 1)
ROs_isoquant_tpm <-read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/OUT.gene_grouped_tpm.tsv",
                              sep = "\t", row.names = 1)


#remove rows that have the same value across all columns
removeZeroVarRows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}

FT_vs_RGC_bambu_tpm <- convertCounts(as.matrix(FT_vs_RGC_bambu_tpm),
                                     unit = "CPM",normalize = "TMM"
)
ROs_bambu_tpm <- convertCounts(as.matrix(ROs_bambu_tpm),
                               unit = "CPM",normalize = "TMM"
)
FT_vs_RGC_isoquant_tpm <- convertCounts(as.matrix(FT_vs_RGC_isoquant_tpm),
                                        unit = "CPM",normalize = "TMM"
)
ROs_isoquant_tpm <- convertCounts(as.matrix(ROs_isoquant_tpm),
                                  unit = "CPM",normalize = "TMM"
)

FT_vs_RGC_bambu_tpm <- removeZeroVarRows(FT_vs_RGC_bambu_tpm)
ROs_bambu_tpm <- removeZeroVarRows(ROs_bambu_tpm)
FT_vs_RGC_isoquant_tpm <- removeZeroVarRows(FT_vs_RGC_isoquant_tpm)
ROs_isoquant_tpm <- removeZeroVarRows(ROs_isoquant_tpm)


#samples
FT_vs_RGC_bambu_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_bambu_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
                       "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")
FT_vs_RGC_isoquant_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_isoquant_samples <-c("EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
                         "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2")


#reset column names 
colnames(FT_vs_RGC_bambu_tpm) <-  FT_vs_RGC_bambu_samples
colnames(ROs_bambu_tpm) <- ROs_bambu_samples
colnames(FT_vs_RGC_isoquant_tpm) <- FT_vs_RGC_isoquant_samples
colnames(ROs_isoquant_tpm) <- ROs_isoquant_samples

dge_dir <- "/users/sparthib/retina_lrs/processed_data/dge/edgeR"
#bambu ROs
bambu_D100_vs_D45 <- read_tsv(file.path(dge_dir, "bambu", "ROs", "D100_vs_D45_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
bambu_D200_vs_D100 <- read_tsv(file.path(dge_dir, "bambu", "ROs", "D200_vs_D100_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
bambu_D200_vs_D45 <- read_tsv(file.path(dge_dir, "bambu", "ROs", "D200_vs_D45_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
#rbind, sort by FDR and remove duplicate and select top 200
bambu_top_RO_genes <- rbind(bambu_D100_vs_D45, bambu_D200_vs_D100, bambu_D200_vs_D45) |>
  arrange(FDR) |>
  distinct(gene_id, .keep_all = TRUE) |>
  head(n = 200) |> select(gene_id, gene_name)

#isoquant ROs
isoquant_D100_vs_D45 <- read_tsv(file.path(dge_dir, "isoquant", "ROs", "D100_vs_D45_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
isoquant_D200_vs_D100 <- read_tsv(file.path(dge_dir, "isoquant", "ROs", "D200_vs_D100_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
isoquant_D200_vs_D45 <- read_tsv(file.path(dge_dir, "isoquant", "ROs", "D200_vs_D45_DGEs.tsv")) |>
  filter(FDR < 0.05 & abs(logFC) > 0.5) |>
  arrange(FDR) |>
  head(n = 200)
#rbind, sort by FDR and remove duplicate and select top 200
isoquant_top_RO_genes <- rbind(isoquant_D100_vs_D45, isoquant_D200_vs_D100, isoquant_D200_vs_D45) |>
  arrange(FDR) |>
  distinct(gene_id, .keep_all = TRUE) |>
  head(n = 200) |> select(gene_id, gene_name)




output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/heatmaps"
quant_name <-"isoquant"
compare <-  "ROs"
tmm <- ROs_isoquant_tpm
groups <- ROs_isoquant_groups

plot_heatmap <- function(dge_dir, quant_name, compare, tmm, groups, output_plots_dir) {
  # Read the DEXSeqSwitchList file
  
  if(compare == "FT_vs_RGC"){
    dge_file <- read_tsv(file.path(dge_dir, quant_name, compare,  "DGEs.tsv"))
    dge_file <- dge_file |> 
      filter(FDR < 0.05 & abs(logFC) > 0.5) |>
      arrange(FDR) |> select(gene_id, gene_name) |>
      head(n = 200)
  }
  
  if(compare == "ROs" & quant_name == "bambu"){ 
    dge_file <- bambu_top_RO_genes
  } else if(compare == "ROs" & quant_name == "isoquant"){
    dge_file <- isoquant_top_RO_genes
  }
  
  # remove version number from tmm rownames
  rownames(tmm) <- gsub("\\.\\d+$", "", rownames(tmm))
  # Keep only tmm rows that are in DTU_isoforms based on its rownames
  tmm <- as.data.frame(tmm)
  tmm$gene_id <- rownames(tmm)
  TMM_significant_genes <- tmm |>
    inner_join(dge_file, by = "gene_id")
  
  # Create a new column with gene name and isoform name
  TMM_significant_genes <- TMM_significant_genes |> 
    mutate(gene = paste(gene_name, gene_id, sep = "_")) |>
    column_to_rownames(var = "gene") |>
    select(-gene_id, -gene_name)
  
  
  # Convert to matrix and scale
  TMM_matrix <- as.matrix(TMM_significant_genes)
  TMM_matrix <- t(scale(t(TMM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Heatmap annotation
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen")
                            ))
  }
  
  
  # Create PDF output
  pdf(file.path(output_plots_dir, quant_name, paste0(compare, "/top_genes_heatmap.pdf")))
  
  # Generate and draw heatmap
  ht_list <- Heatmap(TMM_matrix, name = paste0("scaled TMM expression of top", str(nrow(TMM_matrix))," DE genes"), row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
}

plot_heatmap(dge_dir, "bambu", "ROs", ROs_bambu_tpm, ROs_bambu_groups, output_plots_dir)
plot_heatmap(dge_dir, "isoquant", "ROs", ROs_isoquant_tpm, ROs_isoquant_groups, output_plots_dir)
plot_heatmap(dge_dir, "bambu", "FT_vs_RGC", FT_vs_RGC_bambu_tpm, FT_vs_RGC_bambu_groups, output_plots_dir)
plot_heatmap(dge_dir, "isoquant", "FT_vs_RGC", FT_vs_RGC_isoquant_tpm, FT_vs_RGC_isoquant_groups, output_plots_dir)
