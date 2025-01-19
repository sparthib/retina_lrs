library(dplyr)
library(readr)

method <- "bambu"
comparison <- "ROs"

if( comparison == "ROs") {
  groups  = c("Stage_1", "Stage_1", "Stage_2","Stage_2",
              "Stage_2", "Stage_3","Stage_3" )
} 


isoformFeatures <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                      method, comparison, "isoformFeatures.tsv"))

isoformFeatures <- isoformFeatures |> distinct()

colnames(isoformFeatures)

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

DTE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                method, comparison, "DTE_table.tsv"))
DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                method, comparison, "DGE_table.tsv"))

DTE_table$isoform_id <- ifelse(
  grepl("^ENST", DTE_table$isoform_id),  # Check if isoform_id starts with "ENST"
  gsub("\\..*", "", DTE_table$isoform_id),  # Remove everything after the first dot
  DTE_table$isoform_id  # Keep other isoform_id values unchanged
)

colnames(DTE_table)
if (comparison == "FT_vs_RGC") {
  DTE_table$condition_1 <- "FT"
  DTE_table$condition_2 <- "RGC"
}
DTE_table <- DTE_table |> dplyr::select(c( "isoform_id", "logFC", "logCPM", "F", "PValue", "FDR",  "condition_1", "condition_2"))
colnames(DTE_table) <- c("isoform_id", "DTE_log2FC", "DTE_logCPM", "DTE_F", "DTE_pval", "DTE_qval", "condition_1", "condition_2")


DGE_table$gene_id <- gsub("\\..*", "", DGE_table$gene_id)
colnames(DGE_table)
if (comparison == "FT_vs_RGC") {
  DGE_table$condition_1 <- "FT"
  DGE_table$condition_2 <- "RGC"
}
DGE_table <- DGE_table |> dplyr::select(c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR", "condition_1", "condition_2"))
colnames(DGE_table) <- c("gene_id", "DGE_log2FC", "DGE_logCPM", "DGE_F", "DGE_pval", "DGE_qval", "condition_1", "condition_2")


DTE_table <- DTE_table |> distinct()
new_DGE_DTE_DTU <- left_join(isoformFeatures, DTE_table, 
            by = c("isoform_id", "condition_1", "condition_2"))

new_DGE_DTE_DTU <- left_join(new_DGE_DTE_DTU, DGE_table, 
            by = c("gene_id", "condition_1", "condition_2"))


# Add DTU column based on isoform_switch_q_value and dIF
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTU = isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1)

# Add DTE column based on DTE_qval and DTE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTE = DTE_qval < 0.05 & abs(DTE_log2FC) >= 1)

# Add DGE column based on DGE_qval and DGE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DGE = DGE_qval < 0.05 & abs(DGE_log2FC) >= 1)

table(new_DGE_DTE_DTU$DTU, useNA = "always")
table(new_DGE_DTE_DTU$DGE, useNA = "always")
table(new_DGE_DTE_DTU$DTE, useNA = "always")

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

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)

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

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> dplyr::select(-c(gene_name, gene_biotype, isoform_biotype)) |> 
  left_join(annotLookup, by = "isoform_id")
new_DGE_DTE_DTU
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison)
write_tsv(new_DGE_DTE_DTU, file.path(input_data_dir, "DGE_DTE_DTU.tsv"))







