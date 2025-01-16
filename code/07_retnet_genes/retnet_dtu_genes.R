library(readr)
library(tidyr)
library(dplyr)
library(readxl)
library(biomaRt)
library(stringr)
library(purrr)
raw_data_dir <- "/users/sparthib/retina_lrs/raw_data"


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
gene_vector <- c(
  "ADIPOR1", "ARL6", "BBIP1", "BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBS10", "BBS12", "C8orf37", "CEP19", "CEP290", "IFT172", "IFT27", "INPP5E", "LZTFL1", "MKKS", "MKS1", "NPHP1", "SDCCAG8", "TRIM32", "TTC8",
  "PRDM13", "RGR", "TEAD1",
  "AIPL1", "CRX", "GUCA1A", "GUCY2D", "PITPNM3", "PROM1", "PRPH2", "RIMS1", "SEMA4A", "UNC119",
  "ABCA4", "ADAM9", "ATF6", "C21orf2", "C8orf37", "CACNA2D4", "CDHR1", "CEP78", "CERKL", "CNGA3", "CNGB3", "CNNM4", "DYNC2I2", "GNAT2", "IFT81", "KCNV2", "PDE6C", "PDE6H", "POC1B", "RAB28", "RAX2", "RDH5", "RPGRIP1", "SLC4A7", "TTLL5",
  "CACNA1F", "RPGR",
  "GNAT1", "PDE6B", "RHO",
  "CABP4", "GNAT1", "GNB3", "GPR179", "GRK1", "GRM6", "LRIT3", "RDH5", "SAG", "SLC24A1", "TRPM1",
  "CACNA1F", "NYX",
  "ESPN", "WFS1",
  "CDH23", "CIB2", "ESPN", "MYO7A", "PCDH15", "PDZD7", "USH1C", "WHRN",
  "CRX", "IMPDH1", "OTX2",
  "AIPL1", "CABP4", "CCT2", "CEP290", "CLUAP1", "CRB1", "CRX", "DTHD1", "GDF6", "GUCY2D", "IFT140", "IQCB1", "KCNJ13", "LCA5", "LRAT", "NMNAT1", "PRPH2", "RD3", "RDH12", "RPE65", "RPGRIP1", "SPATA7", "TULP1",
  "BEST1", "C1QTNF5", "CTNNA1", "EFEMP1", "ELOVL4", "FSCN2", "GUCA1B", "HMCN1", "IMPG1", "LRRTM4", "OTX2", "PRDM13", "PROM1", "PRPH2", "RP1L1", "TIMP3",
  "ABCA4", "CFH", "DRAM2", "IMPG1", "MFSD8", "SLC37A3",
  "RPGR", "VCAN",
  "AFG3L2", "MFN2", "MIEF1", "NR2F1", "OPA1",
  "ACO2", "NBAS", "RTN4IP1", "TMEM126A", "TIMM8A",
  "ADIPOR1", "ARL3", "BEST1", "CA4", "CRX", "FSCN2", "GUCA1B", "HK1", "IMPDH1", "IMPG1", "KIF3B", "KLHL7", "NR2E3", "NRL", "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF31", "PRPH2", "RDH12", "RHO", "ROM1", "RP1", "RP9", "RPE65", "SAG", "SEMA4A", "SNRNP200", "SPP2", "TOPORS",
  "ABCA4", "ADGRA3", "AGBL5", "AHR", "ARHGEF18", "ARL6", "ARL2BP", "BBS1", "BBS2", "BEST1", "C8orf37", "CC2D2A", "CERKL", "CLCC1", "CLRN1", "CNGA1", "CNGB1", "COQ2", "COQ4", "COQ5", "CRB1", "CWC27", "CYP4V2", "DHDDS", "DHX38", "EMC1", "ENSA", "EYS", "FAM161A", "HGSNAT", "IDH3B", "IFT140", "IFT172", "IMPG2", "KIAA1549", "KIZ", "LRAT", "MAK", "MERTK", "MVK", "NEK2", "NEUROD1", "NR2E3", "NRL", "PCARE", "PDE6A", "PDE6B", "PDE6G", "PDSS1", "POMGNT1", "PRCD", "PROM1", "PROS1", "RAX2", "RBP3", "REEP6", "RGR", "RHO", "RLBP1", "RP1", "RP1L1", "RPE65", "SAG", "SAMD11", "SLC37A3", "SLC39A12", "SLC66A1", "SLC7A14", "SPATA7", "TRNT1", "TTC8", "TULP1", "USH2A", "ZNF408", "ZNF513",
  "OFD1", "RP2", "RPGR",
  "ABCC6", "AFG3L2", "ATXN7", "COL11A1", "COL2A1", "JAG1", "KCNJ13", "KIF11", "MFN2", "OPA3", "PAX2", "TREX1", "VCAN",
  "ABCC6", "ABHD12", "ACBD5", "ACO2", "ADAMTS18", "ADIPOR1", "AFG3L2", "AHI1", "ALMS1", "CC2D2A", "CEP164", "CEP290", "CLN3", "COL9A1", "CSPP1", "CWC27", "ELOVL4", "EXOSC2", "FLVCR1", "GNPTG", "HARS", "HGSNAT", "HMX1", "IFT140", "IFT81", "INPP5E", "INVS", "IQCB1", "LAMA1", "LRP5", "MKS1", "MTTP", "NPHP1", "NPHP3", "NPHP4", "OPA3", "PANK2", "PCYT1A", "PDSS1", "PEX1", "PEX2", "PEX7", "PHYH", "PLK4", "PNPLA6", "POC5", "POC1B", "PPT1", "PRPS1", "RDH11", "RIMS2", "RPGRIP1L", "SDCCAG8", "SLC25A46", "TMEM216", "TMEM237", "TRNT1", "TTPA", "TUB", "TUBGCP4", "TUBGCP6", "WDPCP", "WDR19", "WFS1", "ZNF423",
  "OFD1", "TIMM8A",
  "ABHD12", "ADGRV1", "ARSG", "CDH23", "CEP250", "CEP78", "CIB2", "CLRN1", "ESPN", "HARS", "MYO7A", "PCDH15", "USH1C", "USH1G", "USH2A", "WHRN",
  "BEST1", "CAPN5", "CRB1", "ELOVL1", "FZD4", "ITM2B", "KIF3B", "LRP5", "MAPKAPK3", "MIR204", "OPN1SW", "RB1", "RCBTB1", "RGR", "TSPAN12", "ZNF408",
  "ASRGL1", "BEST1", "C12orf65", "CDH3", "CNGA3", "CNGB3", "CNNM4", "COQ2", "CYP4V2", "DYNC2H1", "LRP5", "MFRP", "MVK", "NBAS", "NR2E3", "OAT", "PLA2G5", "PROM1", "RBP4", "RCBTB1", "RGS9", "RGS9BP", "RLBP1",
  "KSS", "LHON", "MT-ATP6", "MT-TH", "MT-TL1", "MT-TP", "MT-TS2",
  "CACNA1F", "CHM", "DMD", "NDP", "OPN1LW", "OPN1MW", "PGK1", "RS1"
)

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


#convert results_list to a dataframe
results_df <- data.frame(
  gene_name = names(results_list),
  disease_category = as.character(results_list),
  stringsAsFactors = FALSE
)



save_retnet_dtu_genes <- function(method, comparison) {
  # Define input and output directories
  input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison)
  plots_dir <- file.path(input_data_dir, "plots", "retnet")
  if (!dir.exists(plots_dir)){
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Read and filter data
  retn <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv")) |> 
    dplyr::filter(DTU == TRUE)
  
  # Merge with results_df
  retn <- dplyr::inner_join(retn, results_df, by = "gene_name")

  write_tsv(retn, file.path(plots_dir, paste0(comparison, "_retnet_DTU_genes.tsv")))

}

parameters <- list(
  list("bambu", "ROs"),
  list("bambu", "FT_vs_RGC"),
  list("Isoquant", "ROs"),
  list("Isoquant", "FT_vs_RGC")
)

# Execute save_retnet_dtu_genes for each parameter set
purrr::walk(parameters, ~ save_retnet_dtu_genes(.x[[1]], .x[[2]]))


#Prepare for heatmap 


# Read all the files
files <- list.files("/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots/retnet", full.names = TRUE)



plot_heatmap <- function(method, comparison){ 
  
  
  co

  
  }



load_gene_counts_matrix <- function(analysis_type, quant_method, counts_matrix_dir, splicing_factors_path,
                                    table_type = "DTE") {
  # Validate inputs
  if (!analysis_type %in% c("FT_vs_RGC", "ROs")) {
    stop("Invalid analysis_type. Choose 'FT_vs_RGC' or 'ROs'.")
  }
  if (!quant_method %in% c("bambu", "Isoquant")) {
    stop("Invalid quant_method. Choose 'bambu' or 'isoquant'.")
  }
  
  input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", quant_method, analysis_type)
  plots_dir <- file.path(input_data_dir, "plots", "retnet")
  DGE_DTE_DTU <- read_tsv(file.path(plots_dir, paste0(analysis_type, "_retnet_DTU_genes.tsv")))

  gene_file <- file.path(counts_matrix_dir, quant_method, analysis_type, "gene_cpm.RDS")
  
  # Load and process gene TPM
  gene_tpm <- readRDS(gene_file)
  rownames(gene_tpm) <- gsub("\\..*", "", rownames(gene_tpm))
  gene_tpm <- gene_tpm[rownames(gene_tpm) %in%  DGE_DTE_DTU$gene_id, ]
  
  input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", quant_method, analysis_type)
  DGE_DTE_DTU <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))
  
  gene_tpm <- gene_tpm[rownames(gene_tpm) %in% significant_genes, ]
  gene_tpm <- remove_zero_var_rows(gene_tpm)
  
  # Define groups and output directory
  groups <- if (analysis_type == "ROs") {
    c("RO_D45", "RO_D45", "RO_D100", "RO_D100", "RO_D100", "RO_D200", "RO_D200")
  } else {
    c("FT", "FT", "RGC", "RGC")
  }
  
  output_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu", 
                                quant_method, analysis_type, "plots", "retnet")
  return(list( gene_tpm = gene_tpm, 
               splicing_factors = splicing_factors, groups = groups, 
               output_plots_dir = output_plots_dir))
}


