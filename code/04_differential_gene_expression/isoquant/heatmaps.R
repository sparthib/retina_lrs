#####clustering based on only top isoforms that show DTU. 
library(readr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library("DGEobj.utils")
library(grid)

## DexSeqSwitchList.tsv


#bambu
FT_vs_RGC_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/CPM_transcript.txt",
                                  sep = "\t", row.names = 1, header = TRUE)
ROs_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/CPM_transcript.txt",
                                   sep = "\t", row.names = 1, header = TRUE)
FT_vs_RGC_isoquant_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/OUT.transcript_model_grouped_tpm.tsv",
                                     sep = "\t", row.names = 1)
ROs_isoquant_tpm <-read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/OUT.transcript_model_grouped_tpm.tsv",
                                     sep = "\t", row.names = 1)
FT_vs_RGC_bambu_tpm <- FT_vs_RGC_bambu_tpm[,2:5]    
ROs_bambu_tpm <- ROs_bambu_tpm[,2:8]

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

#groups 
FT_vs_RGC_bambu_groups <- c("FT", "FT", "RGC", "RGC")
ROs_bambu_groups <- c("RO_D200", "RO_D100", "RO_D45", 
                             "RO_D100", "RO_D200", "RO_D100", "RO_D45")
FT_vs_RGC_isoquant_groups <- c("FT", "FT", "RGC", "RGC")
ROs_isoquant_groups <- c("RO_D200", "RO_D45", "RO_D100",
                                "RO_D200", "RO_D100", "RO_D45", "RO_D100")


RO_D100_vs_RO_D45, RO_D200_vs_RO_D45, FT_vs_RGC, RO_D100_vs_RO_D200, DEXSeqSwitchList.tsv
#### Hierarchical Clustering #####

output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/heatmaps"
dex_seq_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/"

#ROs isoform dexseq 
bambu_RO_D100_vs_RO_D45 <- read_tsv(file.path(dex_seq_dir, "bambu", "RO_D100_vs_RO_D45", "DEXSeqSwitchList.tsv"))
bambu_RO_D200_vs_RO_D45 <- read_tsv(file.path(dex_seq_dir, "bambu", "RO_D200_vs_RO_D45", "DEXSeqSwitchList.tsv"))
bambu_RO_D100_vs_RO_D200 <- read_tsv(file.path(dex_seq_dir, "bambu", "RO_D100_vs_RO_D200", "DEXSeqSwitchList.tsv"))

#get values for top 200 isoforms in each comparison
bambu_RO_D100_vs_RO_D45 <- bambu_RO_D100_vs_RO_D45 %>%
  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |> distinct() |> 
  arrange(isoform_switch_q_value) |> 
  head(n = 200)

bambu_RO_D200_vs_RO_D45 <- bambu_RO_D200_vs_RO_D45 %>%
  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |>  distinct() |> 
  arrange(isoform_switch_q_value) |> 
  head(n = 200)

bambu_RO_D100_vs_RO_D200 <- bambu_RO_D100_vs_RO_D200 %>%
  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |> distinct() |> 
  arrange(isoform_switch_q_value) |> 
  head(n = 200)

#concatenate all three comparisons
bambu_top_RO_isoforms <- rbind(bambu_RO_D100_vs_RO_D45, bambu_RO_D200_vs_RO_D45, bambu_RO_D100_vs_RO_D200) |> 
  arrange(isoform_switch_q_value) |> select(isoform_id, ensembl_gene_name) |> distinct() |> head(n = 200)

#what are the common isoforms in all three comparisons
intersect(bambu_RO_D100_vs_RO_D45$isoform_id, bambu_RO_D200_vs_RO_D45$isoform_id)
intersect(bambu_RO_D100_vs_RO_D45$isoform_id, bambu_RO_D100_vs_RO_D200$isoform_id)
intersect(bambu_RO_D200_vs_RO_D45$isoform_id, bambu_RO_D100_vs_RO_D200$isoform_id)


#ROs isoquant dexseq
# isoquant_RO_D100_vs_RO_D45 <- read_tsv(file.path(dex_seq_dir, "isoquant", "RO_D100_vs_RO_D45", "DEXSeqSwitchList.tsv"))
isoquant_RO_D200_vs_RO_D45 <- read_tsv(file.path(dex_seq_dir, "isoquant", "RO_D200_vs_RO_D45", "DEXSeqSwitchList.tsv"))
isoquant_RO_D100_vs_RO_D200 <- read_tsv(file.path(dex_seq_dir, "isoquant", "RO_D100_vs_RO_D200", "DEXSeqSwitchList.tsv"))

isoquant_top_RO_isoforms <- rbind(isoquant_RO_D200_vs_RO_D45, isoquant_RO_D100_vs_RO_D200) |> 
  arrange(isoform_switch_q_value) |> select(isoform_id, ensembl_gene_name) |> distinct() |> head(n = 200)

quant_name <- "bambu"
compare <- "ROs"

plot_heatmap <- function(dex_seq_dir, quant_name, compare, tmm, groups, output_plots_dir) {
  # Read the DEXSeqSwitchList file

  if(compare == "FT_vs_RGC"){
    dex_seq_file <- read_tsv(file.path(dex_seq_dir, quant_name, compare, "DEXSeqSwitchList.tsv"))
    dex_seq_file <- dex_seq_file %>%
      filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |>
      arrange(isoform_switch_q_value) |> 
      head(n = 200)
  }
  # Filter significant DTUs
  significant_DTUs <- dex_seq_file %>%
    filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |>
    arrange(isoform_switch_q_value) |> 
    select(isoform_id, ensembl_gene_name) |>
    distinct() |> 
    head(n = 200)
  
  if(compare == "ROs" & quant_name == "bambu"){ 
    significant_DTUs <- bambu_top_RO_isoforms
  } else if(compare == "ROs" & quant_name == "isoquant"){
    significant_DTUs <- isoquant_top_RO_isoforms
  }
  
  # Keep only tmm rows that are in DTU_isoforms based on its rownames
  tmm <- as.data.frame(tmm)
  tmm$isoform_id <- rownames(tmm)
  TMM_significant_isoforms <- tmm |>
    inner_join(significant_DTUs, by = "isoform_id")
  
  # Create a new column with gene name and isoform name
  TMM_significant_isoforms <- TMM_significant_isoforms |> 
    mutate(gene_isoform = paste(ensembl_gene_name, isoform_id, sep = "_")) |>
    column_to_rownames(var = "gene_isoform") |>
    select(-isoform_id, -ensembl_gene_name)
    
  
  # Convert to matrix and scale
  TMM_matrix <- as.matrix(TMM_significant_isoforms )
  TMM_matrix <- t(scale(t(TMM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Heatmap annotation
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  
  # Create PDF output
  pdf(file.path(output_plots_dir, quant_name, paste0(compare, "/top_isoforms_heatmap.pdf")))
  
  # Generate and draw heatmap
  ht_list <- Heatmap(TMM_matrix, name = "scaled TMM expression of top 176 DTU isoforms", row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
}

plot_heatmap(dex_seq_dir, "bambu", "FT_vs_RGC", FT_vs_RGC_bambu_tpm, FT_vs_RGC_bambu_groups, output_plots_dir)
plot_heatmap(dex_seq_dir, "isoquant", "FT_vs_RGC", FT_vs_RGC_isoquant_tpm, FT_vs_RGC_isoquant_groups, output_plots_dir)
plot_heatmap(dex_seq_dir, "bambu", "ROs", ROs_bambu_tpm, ROs_bambu_groups, output_plots_dir)
plot_heatmap(dex_seq_dir, "isoquant", "ROs", ROs_isoquant_tpm, ROs_isoquant_groups, output_plots_dir)









