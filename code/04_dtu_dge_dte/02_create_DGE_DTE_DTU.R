library(dplyr)
library(readr)
library(IsoformSwitchAnalyzeR)
library(biomaRt)

method <- "bambu"
comparison <- "RO_vs_RGC"


data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

isoformFeatures <- read_tsv(file.path(code_dir, "processed_data/dtu/",
                                      method, comparison, "protein_coding", "isoformFeatures_part2.tsv"))

matrix_dir <-  file.path(data_dir, "/06_quantification/counts_matrices/",
                        method, comparison, "filtered_by_counts_and_biotype")
counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)
nrow(counts)


isoformFeatures |> filter(isoform_switch_q_value < 0.05 & 
                            abs(dIF) >= 0.1) |> nrow()

isoformFeatures <- isoformFeatures |> distinct()

range(isoformFeatures$isoform_switch_q_value, na.rm = TRUE)

colnames(isoformFeatures)

nrow(isoformFeatures)/3

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

nrow(isoformFeatures)/3

DTE_table <- read_tsv(file.path(code_dir, "processed_data/dtu",
                                 method, comparison, "protein_coding","DTE_table.tsv"))
DGE_table <- read_tsv(file.path(code_dir, "processed_data/dtu",
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
nrow(DTE_table)

new_DGE_DTE_DTU <- left_join(isoformFeatures, DTE_table, 
            by = c("isoform_id", "condition_1", "condition_2"))

nrow(new_DGE_DTE_DTU)


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

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

#### get isoform biotypes ####
annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_transcript_id", "transcript_biotype"),
  filter="ensembl_transcript_id",
  values=new_DGE_DTE_DTU$isoform_id,
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("isoform_id","isoform_biotype")

annotLookup <- annotLookup |> distinct()

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> dplyr::select(-c( isoform_biotype)) |> 
  left_join(annotLookup, by = "isoform_id")

nrow(new_DGE_DTE_DTU)

#### get gene names and biotypes ####
annotLookup <- getBM(
  mart=mart,
  attributes=c( 
                "external_gene_name",
                "gene_biotype", "ensembl_gene_id"),
  filter="ensembl_gene_id",
  values=new_DGE_DTE_DTU$gene_id,
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id")
annotLookup <- annotLookup |> distinct()

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> dplyr::select(-c(gene_biotype, gene_name)) |> 
  left_join(annotLookup, by = "gene_id")

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> distinct()

input_data_dir <- file.path(code_dir, "processed_data/dtu/", method,
                            comparison, "protein_coding" )
write_tsv(new_DGE_DTE_DTU, file.path(input_data_dir, "DGE_DTE_DTU.tsv"))

