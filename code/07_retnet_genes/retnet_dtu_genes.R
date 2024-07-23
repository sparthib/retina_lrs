library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(rtracklayer)
library(edgeR)
library(tidyr)
library(dplyr)
library(readxl)
library(here)
library(biomaRt)


RetNet_gene_list <- read_excel(here("raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")
genes_and_diseases <- read_excel(here("raw_data", "RetNet.xlsx"),
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
disease_gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = disease_genes$gene,
                 mart = mart)
# nrow(disease_genes)
#get rows with duplicate gene ids in gene_id
duplicated_gene_id <- disease_gene_id[duplicated(disease_gene_id$hgnc_symbol),]

# > duplicated_gene_id 
# ensembl_gene_id hgnc_symbol
# 3   ENSG00000091262       ABCC6
# 8   ENSG00000282230       ADAM9
# 42  ENSG00000151062    CACNA2D4
# 104 ENSG00000277399      GPR179
# 106 ENSG00000185974        GRK1
# 107 ENSG00000281988        GRK1
# 116 ENSG00000215612        HMX1
# 173 ENSG00000129535         NRL
# 214 ENSG00000274894      PRPF31
# 215 ENSG00000274651      PRPF31
# 216 ENSG00000276421      PRPF31
# 217 ENSG00000274144      PRPF31
# 218 ENSG00000275117      PRPF31
# 219 ENSG00000277953      PRPF31
# 220 ENSG00000277707      PRPF31
# 221 ENSG00000275885      PRPF31
# 222 ENSG00000105618      PRPF31
# 226 ENSG00000174231       PRPF8
# 250 ENSG00000183638       RP1L1
# 260 ENSG00000130561         SAG
# 263 ENSG00000054282     SDCCAG8
# 286 ENSG00000134160       TRPM1
  

DGE_DTU_DTE <- readr::read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE.tsv")

DGE_DTU_DTE <- readr::read_tsv("./processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")

DGE_DTU <- DGE_DTU_DTE |> dplyr::select(gene_id) |> distinct()
#get hgnc symbol
DGE_DTU_gene_id <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = DGE_DTU$gene_id,
                 mart = mart)
DGE_DTU_gene_id |> dplyr::filter(hgnc_symbol %in% duplicated_gene_id$hgnc_symbol)
# ensembl_gene_id hgnc_symbol
# 1 ENSG00000168615       ADAM9

# 7   ENSG00000168615       ADAM9
# 8   ENSG00000282230       ADAM9

##### for both ROs and FT_vs_RGC the only gene with duplicated gene id was ADAM9



#remove rows with duplicated gene ids from disease_gene_id
disease_gene_id <- disease_gene_id |> dplyr::filter(!hgnc_symbol %in% duplicated_gene_id$hgnc_symbol)
#add ENSG00000168615       ADAM9 to disease_gene_id
disease_gene_id <- rbind(disease_gene_id, data.frame(ensembl_gene_id = "ENSG00000168615", hgnc_symbol = "ADAM9"))


# for each gene in disease_gene_id$hgnc_symbol, get the corresponding disease category if that string is found in
# genes_and_diseases$mapped_and_identified_genes
results_list <- list()
for (gene in disease_gene_id$hgnc_symbol) {
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

#merge results_df with disease_gene_id
results_df <- merge(results_df, disease_gene_id, by.x = "gene_name", by.y = "hgnc_symbol", all.x = TRUE)
results_df <- results_df |> dplyr::select(-gene_name)


save_retnet_dtu_genes <- function(comparison, cond1, cond2) {
  
  if (comparison == "FT_vs_RGC") {
    retn <- read_tsv(here("processed_data", "dtu", "DTU_gandall",
                          "bambu", comparison, "DGE_DTU_DTE.tsv"))
    retn <- retn |> filter(DTU_qval < 0.05 & abs(DTU_dIF) > 0.1)
    # Merge retn with results_df if needed here
    
    retn <- merge(retn, results_df, by.x = "gene_id", by.y = "ensembl_gene_id", 
                  all.x = TRUE)
    retn <- retn |> filter(!is.na(gene_name)) |> 
      dplyr::select(gene_id, isoform_id, DTU_dIF, DTU_qval, gene_name, disease_category)
    
    retn <- retn |> filter(!is.na(disease_category))
    
    write_tsv(retn, here("processed_data", "dtu", "DTU_gandall",
                         "bambu",comparison, "retnet_DTU", paste0(cond1, "_vs_", cond2, ".tsv")))
    
    
  } else {
    retn <- read_tsv(here("processed_data", "dtu", "DTU_gandall",
                          "bambu", "ROs", "DGE_DTU_DTE.tsv"))
    retn <- retn |> filter(DTU_qval < 0.05 & abs(DTU_dIF) > 0.1) |>
      filter(condition_1 == cond1 & condition_2 == cond2)

    retn <- merge(retn, results_df, by.x = "gene_id", by.y = "ensembl_gene_id", 
                  all.x = TRUE)
    retn <- retn |> filter(!is.na(gene_name)) |> 
      dplyr::select(gene_id, isoform_id, DTU_dIF, DTU_qval, gene_name, disease_category)
    
    retn <- retn |> filter(!is.na(disease_category))
    
    write_tsv(retn, here("processed_data", "dtu", "DTU_gandall",
                         "bambu","ROs", "retnet_DTU", paste0(cond1, "_vs_", cond2, ".tsv")))
  }
  
}


save_retnet_dtu_genes("FT_vs_RGC", "FT", "RGC")
save_retnet_dtu_genes("RO_D100_vs_RO_D45", "RO_D100", "RO_D45")
save_retnet_dtu_genes("RO_D100_vs_RO_D200", "RO_D100", "RO_D200")
save_retnet_dtu_genes("RO_D200_vs_RO_D45", "RO_D200", "RO_D45")


