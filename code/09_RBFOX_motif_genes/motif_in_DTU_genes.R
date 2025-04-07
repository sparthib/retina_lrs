library(readr)
library(dplyr)


motif_locations <- readr::read_csv("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_gene_overlap.csv")

RO_RGC_DTE_DTU_DGE <- readr::read_tsv(file.path(RO_RGC_data_dir, "DGE_DTE_DTU.tsv"))
FT_RGC_DTE_DTU_DGE <- readr::read_tsv(file.path(FT_RGC_data_dir, "DGE_DTE_DTU.tsv"))
ROs_DTE_DTU_DGE <- readr::read_tsv(file.path(ROs_data_dir, "DGE_DTE_DTU.tsv"))


RO_RGC_DTE_DTU_DGE  <- RO_RGC_DTE_DTU_DGE |> filter(DTU == TRUE)
FT_RGC_DTE_DTU_DGE  <- FT_RGC_DTE_DTU_DGE |> filter(DTU == TRUE)
ROs_DTE_DTU_DGE  <- ROs_DTE_DTU_DGE |> filter(DTU == TRUE)

RO_RGC_genes_and_isoforms <- RO_RGC_DTE_DTU_DGE |> 
  select(gene_id, isoform_id) |> 
  distinct() 

FT_RGC_genes_and_isoforms <- FT_RGC_DTE_DTU_DGE |>
  select(gene_id, isoform_id) |> 
  distinct()

ROs_genes_and_isoforms <- ROs_DTE_DTU_DGE |>
  select(gene_id, isoform_id) |> 
  distinct()

# > nrow(RO_RGC_DTE_DTU_DGE)
# [1] 6718
# > nrow(FT_RGC_DTE_DTU_DGE)
# [1] 464
# > nrow(ROs_DTE_DTU_DGE)
# [1] 6664

colnames(motif_locations)
# [1] "seqnames"    "motif_start" "score"       "motif_end"   "gene_id"    
# [6] "gene_start"  "gene_end"    "gene_name"  


## remove gene id version from motif_locations
motif_locations$gene_id <- gsub("\\..*", "", motif_locations$gene_id)

RO_RGC_motif_locations <- motif_locations |> 
  filter(gene_id %in% RO_RGC_DTE_DTU_DGE$gene_id) 

# [1] 203793
nrow(RO_RGC_motif_locations)
FT_RGC_motif_locations <- motif_locations |>
  filter(gene_id %in% FT_RGC_DTE_DTU_DGE$gene_id)
nrow(FT_RGC_motif_locations)
# [1] 31499


ROs_motif_locations <- motif_locations |>
  filter(gene_id %in% ROs_DTE_DTU_DGE$gene_id)
nrow(ROs_motif_locations)
# [1] 191759

FT_RGC_motif_locations <- FT_RGC_motif_locations |>
  group_by(gene_id, gene_name) |>
  summarise(all_motif_starts = paste(motif_start, collapse = ","), 
            .groups = "drop")
saveRDS(FT_RGC_motif_locations, 
        file = "/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/protein_coding/rds/FT_RGC_motif_locations.rds")

RO_RGC_motif_locations <- RO_RGC_motif_locations |>
  group_by(gene_id, gene_name) |>
  summarise(all_motif_starts = paste(motif_start, collapse = ","), 
            .groups = "drop")
saveRDS(RO_RGC_motif_locations, 
        file = "/users/sparthib/retina_lrs/processed_data/dtu/bambu/RO_vs_RGC/protein_coding/rds/RO_RGC_motif_locations.rds")

ROs_motif_locations <- ROs_motif_locations |>
  group_by(gene_id, gene_name) |>
  summarise(all_motif_starts = paste(motif_start, collapse = ","), 
            .groups = "drop")
saveRDS(ROs_motif_locations, 
        file = "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/protein_coding/rds/ROs_motif_locations.rds")


#### Load AS dataframe ####
FT_RGC_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/protein_coding/"
RO_RGC_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/RO_vs_RGC/protein_coding/"
ROs_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/protein_coding/"
FT_RGC_motif_locations <- readRDS(file.path( FT_RGC_data_dir, "rds","FT_RGC_motif_locations.rds"))
RO_RGC_motif_locations <- readRDS(file.path( RO_RGC_data_dir, "rds", "RO_RGC_motif_locations.rds"))
ROs_motif_locations <- readRDS(file.path ( ROs_data_dir, "rds" , "ROs_motif_locations.rds"))

library(IsoformSwitchAnalyzeR)
library(dplyr)

FT_RGC_switchlist_part2_path = file.path(FT_RGC_data_dir, "rds", "SwitchList_part2.rds")
RO_RGC_switchlist_part2_path = file.path(RO_RGC_data_dir, "rds", "SwitchList_part2.rds")
ROs_switchlist_part2_path = file.path(ROs_data_dir, "rds", "SwitchList_part2.rds")

FT_RGC_switchlist_part2 <- readRDS(FT_RGC_switchlist_part2_path)
RO_RGC_switchlist_part2 <- readRDS(RO_RGC_switchlist_part2_path)
ROs_switchlist_part2 <- readRDS(ROs_switchlist_part2_path)

FT_RGC_Alternative_Splicing <- FT_RGC_switchlist_part2$AlternativeSplicingAnalysis
RO_RGC_Alternative_Splicing <- RO_RGC_switchlist_part2$AlternativeSplicingAnalysis
ROs_Alternative_Splicing <- ROs_switchlist_part2$AlternativeSplicingAnalysis


get_AS_data <- function(sample_AS_data, 
                        sample_genes_and_isoforms, 
                        sample_motif_locations, 
                        AS_type, 
                        output_path = NULL) {
  # Dynamically create column names
  start_col <- paste0(AS_type, "_genomic_start")
  end_col <- paste0(AS_type, "_genomic_end")
  
  # Select relevant AS columns and filter non-zero events
  filtered_AS <- sample_AS_data |>
    dplyr::select("isoform_id", dplyr::all_of(c(AS_type, start_col, end_col))) |>
    dplyr::filter(.data[[AS_type]] != 0)
  
  # Remove version number from isoform_id
  filtered_AS$isoform_id <- gsub("\\..*", "", filtered_AS$isoform_id)
  
  # Filter isoforms that exist in genes_and_isoforms
  filtered_AS <- dplyr::filter(filtered_AS, isoform_id %in% sample_genes_and_isoforms$isoform_id)
  
  # Join with gene metadata and motif locations
  filtered_AS <- dplyr::left_join(filtered_AS, sample_genes_and_isoforms, by = "isoform_id")
  filtered_AS <- dplyr::left_join(filtered_AS, sample_motif_locations, by = "gene_id")
  
  # Write to file if output_path is provided
  if (!is.null(output_path)) {
    readr::write_tsv(filtered_AS, output_path)
  }
  
  return(filtered_AS)
}

## IR

ROs_IR <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "IR",
                      file.path(ROs_data_dir, "motif_locations", "ROs_IR_motif_locations.tsv"))
RO_RGC_IR <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "IR",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_IR_motif_locations.tsv"))
FT_RGC_IR <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "IR",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_IR_motif_locations.tsv"))

nrow(ROs_IR)
nrow(RO_RGC_IR)
nrow(FT_RGC_IR)

## ES
ROs_ES <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "ES",
                      file.path(ROs_data_dir, "motif_locations", "ROs_ES_motif_locations.tsv"))
RO_RGC_ES <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "ES",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_ES_motif_locations.tsv"))
FT_RGC_ES <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "ES",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_ES_motif_locations.tsv"))

nrow(ROs_ES)
nrow(RO_RGC_ES)
nrow(FT_RGC_ES)

## A5
ROs_A5 <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "A5",
                      file.path(ROs_data_dir, "motif_locations", "ROs_A5_motif_locations.tsv"))
RO_RGC_A5 <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "A5",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_A5_motif_locations.tsv"))
FT_RGC_A5 <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "A5",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_A5_motif_locations.tsv"))

nrow(ROs_A5)
nrow(RO_RGC_A5)
nrow(FT_RGC_A5)

## A3
ROs_A3 <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "A3",
                      file.path(ROs_data_dir, "motif_locations", "ROs_A3_motif_locations.tsv"))
RO_RGC_A3 <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "A3",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_A3_motif_locations.tsv"))
FT_RGC_A3 <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "A3",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_A3_motif_locations.tsv"))
nrow(ROs_A3)
nrow(RO_RGC_A3)
nrow(FT_RGC_A3)


## ATSS
ROs_ATSS <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "ATSS",
                      file.path(ROs_data_dir, "motif_locations", "ROs_ATSS_motif_locations.tsv"))
RO_RGC_ATSS <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "ATSS",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_ATSS_motif_locations.tsv"))
FT_RGC_ATSS <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "ATSS",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_ATSS_motif_locations.tsv"))
                      
##ATTS
ROs_ATTS <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "ATTS",
                      file.path(ROs_data_dir, "motif_locations", "ROs_ATTS_motif_locations.tsv"))
RO_RGC_ATTS <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "ATTS",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_ATTS_motif_locations.tsv"))
FT_RGC_ATTS <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "ATTS",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_ATTS_motif_locations.tsv"))
nrow(ROs_ATTS)
nrow(RO_RGC_ATTS)
nrow(FT_RGC_ATTS)



## MEE
ROs_MEE <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "MEE",
                      file.path(ROs_data_dir, "motif_locations", "ROs_MEE_motif_locations.tsv"))
RO_RGC_MEE <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "MEE",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_MEE_motif_locations.tsv"))
FT_RGC_MEE <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "MEE",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_MEE_motif_locations.tsv"))
nrow(ROs_MEE)
nrow(RO_RGC_MEE)
nrow(FT_RGC_MEE)


## MES
ROs_MES <- get_AS_data(ROs_Alternative_Splicing, ROs_genes_and_isoforms, ROs_motif_locations, "MES",
                      file.path(ROs_data_dir, "motif_locations", "ROs_MES_motif_locations.tsv"))
RO_RGC_MES <- get_AS_data(RO_RGC_Alternative_Splicing, RO_RGC_genes_and_isoforms, RO_RGC_motif_locations, "MES",
                         file.path(RO_RGC_data_dir, "motif_locations", "RO_RGC_MES_motif_locations.tsv"))
FT_RGC_MES <- get_AS_data(FT_RGC_Alternative_Splicing, FT_RGC_genes_and_isoforms, FT_RGC_motif_locations, "MES",
                         file.path(FT_RGC_data_dir, "motif_locations", "FT_RGC_MES_motif_locations.tsv"))

nrow(ROs_MES)
nrow(RO_RGC_MES)
nrow(FT_RGC_MES)





