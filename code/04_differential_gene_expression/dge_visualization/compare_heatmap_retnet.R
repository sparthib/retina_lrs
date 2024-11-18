library(here)
library(readxl)
library(stringr)

FT_RGC_top_DGE <- c("CORO2B", "PDE4B", "NCAN", "SLIT1", "FGF14−IT1", "DNER", "CELF4",
                    "FBXL16", "CACNA1E", "KCNH7", "GDAP1", "GNG2",
                    "SCRT2", "DISP3", "TFAP2D", "FBLL1", "CREG2",
                    "MTSS1", "NACAD", "RBFOX3", "NDRG4", "POU4F2",
                    "SLC18A2", "EBF3", "CELF3", "CPLX2", "TRIM67",
                    "B4GALNT1", "CNR1", "CAMK2N1", "CELF5", "SLCO5A1",
                    "DACH2", "GAP43", "SH2D3C", "PRPH", "HECW1", "NRSN1", "RIPOR2",
                    "KIF5A", "NEFL", "P2RX3", "KCTD16", "FSTL4", 
                    "L1CAM", "ELAVL3", "SRRM3", "MIR124−2HG", "ONECUT2", "GNAO1", "GDAP1L1",
                    "GNG3", "DPP6", "RGS8", "DUSP4", "ATOH7", "NA", 
                    "C1QL1", "DRD2", "NSG2", "ASIC4", "RUNDC3A", 
                    "CHRNA4", "STMN4", "LINGO1", "NEFM", "PRMT8",
                    "NMNAT2", "RTN1", "TSPAN7", "MICAL1", "ADGRG1", 
                    "SEPTIN3", "PAK5", "PPP1R1A", "RAB6B", "SERPINI1", 
                    "AFF3", "MYO16", "MLLT11", "TUBB2A", "SNCG",
                    "SYT14", "POU3F1", "NUSAP1", "MTCYBP44", "REST", 
                    "NA", "SMOC1", "NA", "BIRC5", "CCN2", "PLK2", "SOX2", "MKI67", 
                    "NA", "ZFP36L1", "PCDH18", "LRATD2", "IFITM1", "ZIC1", "HES1",
                    "CYP1B1", "ADAMTS19", "ZFAND3−DT", "CLEC19A", "CXCL14", "FAT4",
                    "APOE", "P3H4", "TOX−DT", "CCND1", "SPARC", "VIM", "SOSTDC1", "RDH10",
                    "AURKB", "PCLAF", "FNDC3B", "WNT5A", "CCN1", "CRYAB", "MAF", "CGN",
                    "VSX2", "RAB13", "COL5A1", "GAS1", "TPX2", "LEFTY2", "NOTCH3", "FGF19",
                    "SEL1L3", "RREB1", "NPC2", "STK3", "ALDH1A1", "COL4A6", "COL2A1",
                    "FOXO1", "KIF11", "COL1A2", "GLI3", "TPD52L1", "ANLN", "FZD5",
                    "ATP6V1B1", "MFGE8", "P3H2", "NA", "YAP1", "KIF20A", "EFHD2",
                    "BNC2", "KIAA1217", "PTTG1", "CDC25C", "GPX8", "SFRP2", "SELENOP",
                    "MDFI", "AJUBA", "PDGFRA", "MYOCD", "PRSS35", "CDH6", "IGF2", "GXYLT2",
                    "AXL", "LAMB2", "ZFP36L2", "FZD2", "TSC22D1−AS1", "OAF", "SHISA2",
                    "IFITM3", "EVA1B", "CPAMD8", "MYL9", "GJA1", "FKBP10", "CFI", "KRT18",
                    "IFITM2", "FOXN4", "COLEC12", "LYPD1", "TAGLN", "MYO5C", "HSPA8P4", "PDPN",
                    "ADAMTS9−AS2", "UBE2C", "NA", "PLPP3", "CDH7", "NOTCH2", "SDC4", "SERPINH1",
                    "ZBTB20−AS4")

RetNet_gene_list <- read_excel(here("raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")



pattern <- paste(intersect(FT_RGC_top_DGE, RetNet_gene_list$Symbol), collapse = "|") 
# "SPARC"  "COL2A1" "KIF11" 

genes_and_diseases <- read_excel(here("raw_data", "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")

#return Disease Category where the gene is in "Mapped and Identified Genes" 
for (gene in intersect(FT_RGC_top_DGE, RetNet_gene_list$Symbol)) {
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  print(gene)
  print(result)
}


result <- genes_and_diseases |>
  filter(str_detect("Mapped and Identified Genes", pattern))
#if any of the values in a vector 

#### ROs ####

ro_genes <- c("PRRT2", "CUX2", "POLD1", "CERK", "ODC1", "DEPDC1B", "HMMR", 
              "GMNN", "MIS18BP1", "GAL", "PIMREG", "CIP2A", "HMGA2", "H2AC17",
              "DNAJC9", "SPHK1", "RAD50", "EMID1", "KIF4A", "RAB31", "SNCG",
              "FOXN4", "TRIM59", "PLPP3", "HAUS1", "ZBTB16", "TYMS", "IPO5",
              "PSRC1", "TMPO", "PGD", "UCP2", "CAP2", "CBX2", "MCM6", "C21ORF58",
              "H3C4", "CLVS1", "LINC01414", "NUF2", "CENPK", "POU2F2-AS2", "ATOH7",
              "MCM5", "EFHD2", "NCAPH", "LSM4", "IGF2BP1", "MCM2", "NA", "SUV39H2",
              "HMGA2-AS1", "H1-4", "H1-3", "DHFR", "H2AX", "SORCS1", "CYP1B1", "H2AC20",
              "GULP1", "KRT8", "TMPO-AS1", "STK26", "NHSL1-AS1", "HMGA1", "NA", "VSX2",
              "PAUPAR", "KCNMB4", "SYT7", "TGIF1", "MIR17HG", "MYCN", "GPX8", "ELL2",
              "MN1", "DANCR", "ACVR2B", "ISYNA1", "IGFBPL1", "MIR5688", "FGF19",
              "C14ORF132", "CD24", "GAP43", "MYBL1", "TAGLN3", "FMC1", "CNTN2", "PPIAP41",
              "NA", "NA", "DHRS3", "RXRG", "SLC38A5", "IGSF21", "PFKP", "PDE6H", "SNAP25",
              "PDC", "GNB3", "ALDOC", "ESPN", "TAFA3", "MPP4", "EOGT", "SLC17A7", "SH3GL2", "MAOA", "PYGL", "NAV2", "SYP", "EGFLAM",
              "LMOD1", "KCNH6", "DIRAS2", "SEPTIN4", "CHST9", "CRX", "CLIP4", "ARFGEF3", 
              "CYP26B1", "LARGE1", "AANAT", "NA", "USP2", "FAM107A", "NCOA7", "EML5", "CDKL5",
              "TULP1", "SLC38A3", "EPB41L2", "ARAP2", "NTM", "SYNE1", "UNC13C", "EGF", 
              "GABRG3-AS1", "LHFPL3", "SLC1A3", "APOE", "NA", "UNC80", "CDKN2B-AS1", 
              "YPEL2", "NRXN3", "FLT1", "STAC2", "C8ORF34", "EPAS1", "SVOP", "DLG2", 
              "PLEKHB1", "RD3", "PACSIN1", "MYO3B", "DOCK8", "GUCA1ANB-GUCA1A", "ABCA5",
              "ACSL6", "GNGT2", "RPGRIP1", "SNCB", "ABCA8", "UNC119", "MCF2L2", "FRMPD4",
              "RS1", "GUCA1B", "ABCA4", "FSTL5", "CNGB3", "AIPL1", "CDHR1", "FRMPD1", "SV2B",
              "ARHGAP42", "ABCA7", "ZNF385B", "RAX2", "COBLL1", "KIFC3", "ARR3", "DSCAML1",
              "NA", "USH2A", "NA", "ANO2", "ROM1", "YWHAQP4", "KCNB1", "SLC1A2", "IMPDH1",
              "GNAT1", "PRPH2", "RP1", "CNGB1", "RDH12")
intersect(ro_genes, RetNet_gene_list$Symbol)

# Initialize an empty list to store results
results_list <- list()

# Loop through each gene
for (gene in intersect(ro_genes, RetNet_gene_list$Symbol)) {
  # Filter the data
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  
  # Append results to the list
  results_list[[gene]] <- result$disease_category
}

# Convert the list to a dataframe
results_df <- data.frame(
  Gene = names(results_list),
  Disease_Category = I(results_list),
  stringsAsFactors = FALSE
)

# If you want to view the dataframe, you can use:
print(results_df)


## all rows with Retinitis Pigmentosa in the Disease Category column
 