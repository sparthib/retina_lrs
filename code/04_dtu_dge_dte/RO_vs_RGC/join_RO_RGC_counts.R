library(edgeR)
library(readr)

ROs_input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype/genes_counts.RDS"
FT_vs_RGC_input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype/genes_counts.RDS"
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype/"
dir.create(output_dir, recursive = TRUE)

ROs_mat <- readRDS(ROs_input_dir)
FT_vs_RGC_mat <- readRDS(FT_vs_RGC_input_dir)


nrow(ROs_mat)
nrow(FT_vs_RGC_mat)

genes_not_in_FT_vs_RGC <- setdiff(rownames(ROs_mat), rownames(FT_vs_RGC_mat))
# 1950
genes_not_in_ROs <- setdiff(rownames(FT_vs_RGC_mat), rownames(ROs_mat))
#80

#add the missing genes to the matrix

# create 1950 rows of zeros
zeros_FT_vs_RGC <- matrix(0, nrow = length(genes_not_in_FT_vs_RGC), ncol = ncol(FT_vs_RGC_mat))
rownames(zeros_FT_vs_RGC) <- genes_not_in_FT_vs_RGC
colnames(zeros_FT_vs_RGC) <- colnames(FT_vs_RGC_mat)
FT_vs_RGC_mat <- rbind(FT_vs_RGC_mat, zeros_FT_vs_RGC)

# create 80 rows of zeros
zeros_ROs <- matrix(0, nrow = length(genes_not_in_ROs), ncol = ncol(ROs_mat))
rownames(zeros_ROs) <- genes_not_in_ROs
colnames(zeros_ROs) <- colnames(ROs_mat)
ROs_mat <- rbind(ROs_mat, zeros_ROs)

# check that the number of rows is the same
nrow(ROs_mat)
nrow(FT_vs_RGC_mat)

# cbind after aligning by the same row names
common_rows <- intersect(rownames(ROs_mat), rownames(FT_vs_RGC_mat))
FT_vs_RGC_mat <- FT_vs_RGC_mat[common_rows, ]
ROs_mat <- ROs_mat[common_rows, ]
RO_vs_RGC_gene_counts <- cbind(ROs_mat, FT_vs_RGC_mat)

##### Isoform counts ######
ROs_isoform_input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype/isoform_counts.RDS"
FT_vs_RGC_isoform_input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype/isoform_counts.RDS"


ROs_isoform_mat <- readRDS(ROs_isoform_input_dir)
FT_vs_RGC_isoform_mat <- readRDS(FT_vs_RGC_isoform_input_dir)

nrow(ROs_isoform_mat)
nrow(FT_vs_RGC_isoform_mat)

isoforms_not_in_FT_vs_RGC <- setdiff(rownames(ROs_isoform_mat), rownames(FT_vs_RGC_isoform_mat))
#15169
isoforms_not_in_ROs_isoform <- setdiff(rownames(FT_vs_RGC_isoform_mat), rownames(ROs_isoform_mat))
#1042

#add missing isoforms to the matrix
# create 15169 rows of zeros
zeros_FT_vs_RGC_isoform <- matrix(0, nrow = length(isoforms_not_in_FT_vs_RGC), ncol = ncol(FT_vs_RGC_isoform_mat))
rownames(zeros_FT_vs_RGC_isoform) <- isoforms_not_in_FT_vs_RGC
colnames(zeros_FT_vs_RGC_isoform) <- colnames(FT_vs_RGC_isoform_mat)
FT_vs_RGC_isoform_mat <- rbind(FT_vs_RGC_isoform_mat, zeros_FT_vs_RGC_isoform)

# create 1042 rows of zeros
zeros_ROs_isoform <- matrix(0, nrow = length(isoforms_not_in_ROs_isoform), ncol = ncol(ROs_isoform_mat))
rownames(zeros_ROs_isoform) <- isoforms_not_in_ROs_isoform
colnames(zeros_ROs_isoform) <- colnames(ROs_isoform_mat)
ROs_isoform_mat <- rbind(ROs_isoform_mat, zeros_ROs_isoform)

# check that the number of rows is the same
nrow(ROs_isoform_mat)
nrow(FT_vs_RGC_isoform_mat)

# cbind after aligning by the same row names
common_rows_isoform <- intersect(rownames(ROs_isoform_mat), rownames(FT_vs_RGC_isoform_mat))
FT_vs_RGC_isoform_mat <- FT_vs_RGC_isoform_mat[common_rows_isoform, ]
ROs_isoform_mat <- ROs_isoform_mat[common_rows_isoform, ]
RO_vs_RGC_isoform_counts <- cbind(ROs_isoform_mat, FT_vs_RGC_isoform_mat)

# save the counts matrices
saveRDS(RO_vs_RGC_gene_counts, file = paste0(output_dir, "RO_vs_RGC_gene_counts.RDS"))
saveRDS(RO_vs_RGC_isoform_counts, file = paste0(output_dir, "RO_vs_RGC_isoform_counts.RDS"))

# Define directories
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype/"

# Desired column order
new_order <- c("H9_CRX_ROs_D45", "EP1_WT_ROs_D45", "EP1_WT_hRO_2", 
               "H9_BRN3B_hRO_2", "H9_CRX_hRO_2", "EP1_BRN3B_RO", 
               "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2","H9_FT_1", "H9_FT_2")

# Function to process data tables
process_table <- function(filepath, file_prefix, remove_gene_id = FALSE) {
  data <- read.table(filepath, header = TRUE, row.names = 1, sep = "\t", comment.char = "")
  colnames(data) <- gsub("_primary_over_30_(chr_only_sorted|sorted)", "", colnames(data))
  colnames(data) <- gsub("\\.", "_", colnames(data))
  if (remove_gene_id) {
    data <- data[, -1]
  }
  data <- data[, new_order]
  rownames(data) <- ifelse(
    grepl("^ENST", rownames(data)),  # Check if isoform_id starts with "ENST"
    gsub("\\..*", "", rownames(data)),  # Remove everything after the first dot
    rownames(data))
  # Split into subsets
  RO_vs_RGC <- data[, 1:9]
  # Save subsets
  dir.create(file.path(output_dir, "RO_vs_RGC"), showWarnings = FALSE)
  saveRDS(RO_vs_RGC, file.path(output_dir, "RO_vs_RGC", paste0(file_prefix, "_RO_vs_RGC.RDS")))
}

process_table(file.path(bambu_dir, "counts_transcript.txt"),
              "isoform_counts", remove_gene_id = TRUE)



##### filter by common novel isoforms between bambu and isoquant ######

common_isoforms <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)

common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                                      !grepl("^BambuGene", isoquant_gene_id))


colnames(common_isoforms)

isoform_counts <- readRDS(file.path(output_dir, "RO_vs_RGC", "isoform_counts_RO_vs_RGC.RDS"))
print(paste("pre filter", nrow(isoform_counts)))

isoform_counts <- isoform_counts[
  (grepl("^Bambu", rownames(isoform_counts)) & rownames(isoform_counts) %in% common_isoforms$ref_id) |
    (grepl("^transcript", rownames(isoform_counts)) & rownames(isoform_counts) %in% common_isoforms$isoquant_isoform_id) |
    grepl("^ENST", rownames(isoform_counts)),
]
print(paste("post filter", nrow(isoform_counts)))


gene_counts <- readRDS(file.path(output_dir, "RO_vs_RGC", "gene_counts_RO_vs_RGC.RDS"))
gene_counts <- gene_counts[grepl("^ENSG", rownames(gene_counts)),]  


###### filter by counts and biotype ######

# Function to remove rows with zero variance
remove_zero_variance <- function(mat) {
  mat[rowSums(mat != mat[, 1]) > 0, , drop = FALSE]  # Keep rows with at least one differing value
}


gene_counts <- remove_zero_variance(gene_counts)
nrow(gene_counts)
isoform_counts <- remove_zero_variance(isoform_counts)
nrow(isoform_counts)


##### only keep protein coding genes ##### 
library(biomaRt)

# Function to load and filter gene counts or CPM data
filter_genes <- function(counts, mart) {
  genes <- counts
  genes$gene_nums <- gsub("\\..*", "", rownames(genes))  # Remove Ensembl version numbers
  
  gene_annotLookup <- getBM(
    mart = mart,
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    filter = "ensembl_gene_id",
    values = genes$gene_nums,
    uniqueRows = TRUE
  )
  
  # Keep only protein-coding genes
  genes <- genes[genes$gene_nums %in% gene_annotLookup$ensembl_gene_id[gene_annotLookup$gene_biotype == "protein_coding"], ]
  
  #remove genes$gene_nums column
  genes$gene_nums <- NULL

  return(genes)
}

# Function to load and filter isoform counts or CPM data
filter_isoforms <- function(counts, mart) {
  isoforms <- counts
  isoforms$isoform_nums <- gsub("\\..*", "", rownames(isoforms))  # Remove Ensembl version numbers
  
  isoform_annotLookup <- getBM(
    mart = mart,
    attributes = c("ensembl_transcript_id", "external_gene_name", "gene_biotype", "transcript_biotype"),
    filter = "ensembl_transcript_id",
    values =  isoforms$isoform_nums,
    uniqueRows = TRUE
  )
  
  # Identify isoforms belonging to protein-coding genes
  protein_coding_isoforms <- isoform_annotLookup$ensembl_transcript_id[isoform_annotLookup$gene_biotype == "protein_coding"]
  
  # Identify isoforms with "Bambu" in their rownames
  bambu_isoforms <- grepl("Bambu", rownames(isoforms))
  
  # Keep only isoforms that belong to protein-coding genes or contain "Bambu"
  isoforms <- isoforms[ isoforms$isoform_nums %in% protein_coding_isoforms | bambu_isoforms, ]
  
  # Remove isoforms$isoform_nums column
  isoforms$isoform_nums <- NULL
  return(isoforms)
}


us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)  


genes <- filter_genes(gene_counts, mart)
isoforms <- filter_isoforms(isoform_counts, mart)

groups <- c("Stage_1", "Stage_1", "Stage_2", "Stage_2",
            "Stage_2", "Stage_3", "Stage_3", "RGC1", "RGC2")

convert_to_cpm <- function(counts_matrix, group) {
  # Create a DGEList object
  dge <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, ]
  
  # Calculate normalization factors using TMM (Trimmed Mean of M-values)
  dge <- calcNormFactors(dge)
  
  # Convert to CPM (Counts Per Million)
  cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
  
  return(list(filtered_counts = counts_matrix[keep, ], cpm = cpm_matrix))
}

# Convert gene counts to CPM
gene_counts_result <- convert_to_cpm(genes, groups)
gene_counts <- gene_counts_result$filtered_counts
nrow(gene_counts)
# 18016
gene_cpm <- gene_counts_result$cpm
nrow(gene_cpm)
# 53668

#Convert isoform counts to CPM
isoform_counts_result <- convert_to_cpm(isoforms, groups)
isoform_counts <- isoform_counts_result$filtered_counts
isoform_cpm <- isoform_counts_result$cpm

# Save the filtered counts and CPM matrices
saveRDS(gene_counts, file.path(output_dir,  "filtered_gene_counts.RDS"))
saveRDS(isoform_counts, file.path(output_dir, "filtered_isoform_counts.RDS"))
saveRDS(gene_cpm, file.path(output_dir, "filtered_gene_cpm.RDS"))
saveRDS(isoform_cpm, file.path(output_dir, "filtered_isoform_cpm.RDS"))


