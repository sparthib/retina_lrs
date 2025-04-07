library(dplyr)
library(readr)
library(IsoformSwitchAnalyzeR)
method <- "bambu"
comparison <- "RO_vs_RGC"


switchlist_part2_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                  method, comparison, "protein_coding", "rds", "SwitchList_part2.rds")

SwitchList_part2 <- readRDS(switchlist_part2_path)


isoformFeatures <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                      method, comparison, "protein_coding",  "isoformFeatures.tsv"))

isoformFeatures <- SwitchList_part2$isoformFeatures
isoformFeatures |> filter(isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1) |> nrow()

isoformFeatures <- isoformFeatures |> distinct()


colnames(isoformFeatures)

nrow(isoformFeatures)
# 101101

isoformFeatures <- isoformFeatures |> dplyr::select( 
  gene_id, isoform_id, condition_1, condition_2, 
  gene_name, gene_biotype, iso_biotype,gene_value_1, gene_value_2, 
  iso_value_1, iso_value_2, 
  dIF, IF1, IF2, isoform_switch_q_value
  )

colnames(isoformFeatures) <- c("gene_id", "isoform_id", "condition_1", "condition_2", 
                               "gene_name", "gene_biotype", "isoform_biotype", 
                               "gene_value_1", "gene_value_2", 
                            "iso_value_1", "iso_value_2",  
                            "dIF", "IF1", "IF2", "isoform_switch_q_value")

colnames(isoformFeatures)
isoformFeatures$gene_id <- gsub("\\..*", "", isoformFeatures$gene_id)
isoformFeatures$isoform_id <- ifelse(
  grepl("^ENST", isoformFeatures$isoform_id),  # Check if isoform_id starts with "ENST"
  gsub("\\..*", "", isoformFeatures$isoform_id),  # Remove everything after the first dot
  isoformFeatures$isoform_id  # Keep other isoform_id values unchanged
)


DTE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu",
                                 method, comparison, "protein_coding","DTE_table.tsv"))
DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu",
                                 method, comparison,"protein_coding","DGE_table.tsv"))

DTE_table$isoform_id <- ifelse(
  grepl("^ENST", DTE_table$isoform_id),  # Check if isoform_id starts with "ENST"
  gsub("\\..*", "", DTE_table$isoform_id),  # Remove everything after the first dot
  DTE_table$isoform_id  # Keep other isoform_id values unchanged
)

colnames(DTE_table)
DTE_table <- DTE_table |> dplyr::select(c( "isoform_id", "logFC", "logCPM", "F", "PValue", "FDR",  "condition_1", "condition_2"))
colnames(DTE_table) <- c("isoform_id", "DTE_log2FC", "DTE_logCPM", "DTE_F", "DTE_pval", "DTE_qval", "condition_1", "condition_2")


DGE_table$gene_id <- gsub("\\..*", "", DGE_table$gene_id)
colnames(DGE_table)
DGE_table <- DGE_table |> dplyr::select(c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR", "condition_1", "condition_2"))
colnames(DGE_table) <- c("gene_id", "DGE_log2FC", "DGE_logCPM", "DGE_F", "DGE_pval", "DGE_qval", "condition_1", "condition_2")


DTE_table <- DTE_table |> distinct()


new_DGE_DTE_DTU <- left_join(isoformFeatures, DTE_table, 
            by = c("isoform_id", "condition_1", "condition_2"))

nrow(new_DGE_DTE_DTU)




# diff_isoforms <- setdiff( s2_s3_isoforms$isoform_id, s1_s3_isoforms$isoform_id) 

## found in s2 vs s3 not in other comparisons 
# [1] "BambuTx22"       "ENST00000375108" "ENST00000409110" "ENST00000471884"
# [5] "ENST00000476500" "ENST00000393233" "ENST00000589681" "ENST00000458671"
# [9] "ENST00000589272" "ENST00000447342" "ENST00000350427" "ENST00000393963"
# [13] "ENST00000261047" "ENST00000370657" "BambuTx287"      "ENST00000590311"
# [17] "ENST00000284425" "ENST00000511880" "ENST00000286719" "ENST00000546816"
# [21] "ENST00000546386" "ENST00000428432" "ENST00000335556" "ENST00000703267"
# [25] "ENST00000703266"

# new_DGE_DTE_DTU |> filter(isoform_id %in% diff_isoforms) |> select(gene_name) |> unique()
# 


new_DGE_DTE_DTU <- left_join(new_DGE_DTE_DTU, DGE_table, 
            by = c("gene_id", "condition_1", "condition_2"))

nrow(new_DGE_DTE_DTU)

# Add DTU column based on isoform_switch_q_value and dIF
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTU = isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1)

# Add DTE column based on DTE_qval and DTE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTE = DTE_qval < 0.05 & abs(DTE_log2FC) >= 1)

# Add DGE column based on DGE_qval and DGE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DGE = DGE_qval < 0.05 & abs(DGE_log2FC) >= 1)


#convert NA to FALSE
new_DGE_DTE_DTU$DTU[is.na(new_DGE_DTE_DTU$DTU)] <- FALSE
new_DGE_DTE_DTU$DGE[is.na(new_DGE_DTE_DTU$DGE)] <- FALSE
new_DGE_DTE_DTU$DTE[is.na(new_DGE_DTE_DTU$DTE)] <- FALSE

new_DGE_DTE_DTU$DTU_qval <- new_DGE_DTE_DTU$isoform_switch_q_value


# Add DTE column based on DTE_qval and DTE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTE_0.5 = DTE_qval < 0.05 & abs(DTE_log2FC) >= 0.5) 

# Add DGE column based on DGE_qval and DGE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DGE_0.5 = DGE_qval < 0.05 & abs(DGE_log2FC) >= 0.5)

table(new_DGE_DTE_DTU$DGE_0.5, useNA = "always")
table(new_DGE_DTE_DTU$DTE_0.5, useNA = "always")

colnames(new_DGE_DTE_DTU)

library(biomaRt)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_transcript_id",
                "external_gene_name",
                "gene_biotype", "transcript_biotype"),
  filter="ensembl_transcript_id",
  values=new_DGE_DTE_DTU$isoform_id,
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("isoform_id", "gene_name","gene_biotype",
                           "isoform_biotype")

annotLookup <- annotLookup |> distinct()

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> dplyr::select(-c(gene_name, gene_biotype,
                                                       isoform_biotype)) |> 
  left_join(annotLookup, by = "isoform_id")
nrow(new_DGE_DTE_DTU)
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "protein_coding" )
write_tsv(new_DGE_DTE_DTU, file.path(input_data_dir, "DGE_DTE_DTU.tsv"))

library(dplyr)
df <- readr::read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))





