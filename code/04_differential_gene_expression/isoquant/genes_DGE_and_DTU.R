library(readr)
library(dplyr)

dge_dir <- "/users/sparthib/retina_lrs/processed_data/dge/edgeR/isoquant/"
FT_vs_RGC_dges <- paste0(dge_dir,"FT_vs_RGC/", "DGEs.tsv")
D100_vs_D45_dges <- paste0(dge_dir,"ROs/", "D100_vs_D45_DGEs.tsv")
D100_vs_D200_dges <- paste0(dge_dir,"ROs/", "D200_vs_D100_DGEs.tsv")
D200_vs_D45_dges <- paste0(dge_dir,"ROs/", "D200_vs_D45_DGEs.tsv")

FT_vs_RGC_dges <- read_tsv(FT_vs_RGC_dges)
D100_vs_D45_dges <- read_tsv(D100_vs_D45_dges)
D200_vs_D100_dges <- read_tsv(D100_vs_D200_dges)
D200_vs_D45_dges <- read_tsv(D200_vs_D45_dges)


dtu_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/isoquant/"
FT_vs_RGC_dtus <- paste0(dtu_dir,"FT_vs_RGC/", "DEXSeqSwitchList.tsv")
D100_vs_D45_dtus <- paste0(dtu_dir,"RO_D100_vs_RO_D45/", "DEXSeqSwitchList.tsv")
D100_vs_D200_dtus <- paste0(dtu_dir,"RO_D100_vs_RO_D200/", "DEXSeqSwitchList.tsv")
D200_vs_D45_dtus <- paste0(dtu_dir,"RO_D200_vs_RO_D45/", "DEXSeqSwitchList.tsv")

FT_vs_RGC_dtus <- read_tsv(FT_vs_RGC_dtus)
D100_vs_D45_dtus <- read_tsv(D100_vs_D45_dtus)
D200_vs_D100_dtus <- read_tsv(D100_vs_D200_dtus)
D200_vs_D45_dtus <- read_tsv(D200_vs_D45_dtus)


merge_dge_with_dtu <- function(dge_df, dtu_df){
  dtu_df <- dtu_df |>  
    dplyr::filter(isoform_switch_q_value < 0.05) |> 
    dplyr::filter(abs(dIF) >= 0.1) |> distinct()
  dge_df <- dge_df |> 
    dplyr::filter(FDR < 0.05) |> 
    dplyr::filter(abs(logFC) >= 0.5) |> distinct()
  merged_df <- merge(dge_df, dtu_df, by = "gene_id") |>
    dplyr::select(gene_id, ensembl_gene_name, isoform_id,logFC, FDR, isoform_switch_q_value, dIF)
  return(merged_df)
  
}

FT_vs_RGC_merged <- merge_dge_with_dtu(FT_vs_RGC_dges, FT_vs_RGC_dtus)
nrow(FT_vs_RGC_merged)

D100_vs_D45_merged <- merge_dge_with_dtu(D100_vs_D45_dges, D100_vs_D45_dtus)
nrow(D100_vs_D45_merged)

D200_vs_D100_merged <- merge_dge_with_dtu(D200_vs_D100_dges, D200_vs_D100_dtus)
nrow(D200_vs_D100_merged)

D200_vs_D45_merged <- merge_dge_with_dtu(D200_vs_D45_dges, D200_vs_D45_dtus)
nrow(D200_vs_D45_merged)

write_tsv(FT_vs_RGC_merged, "/users/sparthib/retina_lrs/processed_data/dge_dtu/isoquant/FT_vs_RGC_dtu_dge.tsv")
# write_tsv(D100_vs_D45_merged,"/users/sparthib/retina_lrs/processed_data/dge_dtu/bambu/D100_vs_D45_dtu_dge.tsv")
write_tsv(D200_vs_D100_merged, "/users/sparthib/retina_lrs/processed_data/dge_dtu/isoquant/D200_vs_D45_dtu_dge.tsv")
write_tsv(D200_vs_D45_merged, "/users/sparthib/retina_lrs/processed_data/dge_dtu/isoquant/D200_vs_D100_dtu_dge.tsv")

