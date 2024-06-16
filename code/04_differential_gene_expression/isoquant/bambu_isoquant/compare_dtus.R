library(readr)
library(dplyr)

# "FT_vs_RGC" "RO_D100_vs_RO_D45" "RO_D200_vs_RO_D45" "RO_D100_vs_RO_D200"
compare_datasets <- function(compare){ 
  isoquant_dir <- paste0("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/isoquant/", compare)
  bambu_dir <- paste0("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/", compare)
  output_dir <- paste0("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu_isoquant/", compare)
  
  isoquant_tsv <- read_tsv(file.path(isoquant_dir, "DEXSeqSwitchList.tsv"))
  bambu_tsv <- read_tsv(file.path(bambu_dir, "DEXSeqSwitchList.tsv"))
  
  colnames(isoquant_tsv)
  
  isoquant_data <- isoquant_tsv |> select("gene_id", "isoform_id" ,
                                          "dIF", "isoform_switch_q_value",
                                          "ensembl_gene_name", "condition_1", "condition_2")
  bambu_data <- bambu_tsv |> select("gene_id", "isoform_id" ,
                                    "dIF", "isoform_switch_q_value",
                                    "ensembl_gene_name", "condition_1", "condition_2")   
  #filter out the isoforms with q value < 0.05 and abs(dIF) > 0.1
  isoquant_data <- isoquant_data |> filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)
  bambu_data <- bambu_data |> filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)
  
  colnames(isoquant_data) <- c("gene_id", "isoform_id", "isoquant_dIF", "isoquant_isoform_switch_q_value",
                               "ensembl_gene_name", "isoquant_condition_1", "isoquant_condition_2")
  
  colnames(bambu_data) <- c("gene_id", "isoform_id", "bambu_dIF", "bambu_isoform_switch_q_value",
                            "ensembl_gene_name", "bambu_condition_1", "bambu_condition_2")
  
  
  # inner join the two dataframes based on isoform id
  isoform_switch_data <- inner_join(isoquant_data, bambu_data, by = "isoform_id")
  #remove duplicate rows 
  isoform_switch_data <- isoform_switch_data[!duplicated(isoform_switch_data),]
  #write the data to a file
  write_tsv(isoform_switch_data, file.path(output_dir, "/common_isoforms.tsv"))

  
  }

compare_datasets("FT_vs_RGC")
compare_datasets("RO_D100_vs_RO_D45")
compare_datasets("RO_D200_vs_RO_D45")
compare_datasets("RO_D100_vs_RO_D200")


## load output datasets now
FT_vs_RGC <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu_isoquant/FT_vs_RGC/common_isoforms.tsv")
RO_D100_vs_RO_D45 <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu_isoquant/RO_D100_vs_RO_D45/common_isoforms.tsv")
RO_D200_vs_RO_D45 <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu_isoquant/RO_D200_vs_RO_D45/common_isoforms.tsv")
RO_D100_vs_RO_D200 <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu_isoquant/RO_D100_vs_RO_D200/common_isoforms.tsv")


nrow(RO_D200_vs_RO_D45)
RO_D200_vs_RO_D45$ensembl_gene_name.x

# [1] "GTF2IRD1" "CLK1"     "CLK1"     "PLEKHB1"  "GRAMD1B"  "GUCA1A"  
# [7] "ELN"      "LLGL2"    "DGKD"     "GNB1"     "LIPE"     "PAPOLA"  
# [13] "STRN4"    "SLC22A17" "MYL6"     "MED25"    "CRX"      "EZH2"    
# [19] "CDH23"    "NR2E1"    "OGG1"     "REEP6"    "RTN4"     "PROSER1" 
# [25] "CDKN2C"   "BCL2L12"  "POR"      "DUT"      "ZNF331"   "BIN1"    
# [31] "NUMA1"    "TPM1"     "FTO"      "UCK2"     "GUK1"     "GUK1"    
# [37] "TAGLN3"   "FADS1"    "TEX30"    "NCAPD3"   "FARP1"    "FARP1"   
# [43] "UBE2L6"   "UBE2L6"   "SYNJ1"    "EMC10"    "TRA2A"    "STX3"    
# [49] "STX3"     "TSC22D4"  "TSC22D4"  "LUZP1"    "MAP6"     "MAP6"    
# [55] "CAMTA1"   "BSG"      "BSG"      "RBM4"     "CLTB"     "CLTB"    
# [61] "SEMA4B"   "EVL"      "AP2A1"    "AP2A1"    "RPL12"    "NIFK-AS1"


RO_D100_vs_RO_D200$isoform_id
