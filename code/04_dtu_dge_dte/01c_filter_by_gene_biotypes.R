library(biomaRt)

# Initialize Ensembl biomart connection
initialize_mart <- function() {
  us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
  useDataset("hsapiens_gene_ensembl", us_mart)
}

# Function to remove rows with zero variance
remove_zero_variance <- function(mat) {
  mat[rowSums(mat != mat[, 1]) > 0, , drop = FALSE]  # Keep rows with at least one differing value
}

# Function to load and filter gene counts or CPM data
filter_genes <- function(counts_file, mart) {
  genes <- readRDS(counts_file)
  rownames(genes) <- gsub("\\..*", "", rownames(genes))  # Remove Ensembl version numbers
  
  gene_annotLookup <- getBM(
    mart = mart,
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    filter = "ensembl_gene_id",
    values = rownames(genes),
    uniqueRows = TRUE
  )
  
  # Keep only protein-coding genes
  genes <- genes[rownames(genes) %in% gene_annotLookup$ensembl_gene_id[gene_annotLookup$gene_biotype == "protein_coding"], ]
  
  # Remove rows with zero variance
  genes <- remove_zero_variance(genes)
  
  return(genes)
}

# Function to load and filter isoform counts or CPM data
filter_isoforms <- function(counts_file, mart) {
  isoforms <- readRDS(counts_file)
  rownames(isoforms) <- gsub("\\..*", "", rownames(isoforms))  # Remove Ensembl version numbers
  
  isoform_annotLookup <- getBM(
    mart = mart,
    attributes = c("ensembl_transcript_id", "external_gene_name", "gene_biotype", "transcript_biotype"),
    filter = "ensembl_transcript_id",
    values = rownames(isoforms),
    uniqueRows = TRUE
  )
  
  # Identify isoforms belonging to protein-coding genes
  protein_coding_isoforms <- isoform_annotLookup$ensembl_transcript_id[isoform_annotLookup$gene_biotype == "protein_coding"]
  
  # Identify isoforms with "Bambu" in their rownames
  bambu_isoforms <- grepl("Bambu", rownames(isoforms))
  
  # Keep only isoforms that belong to protein-coding genes or contain "Bambu"
  isoforms <- isoforms[rownames(isoforms) %in% protein_coding_isoforms | bambu_isoforms, ]
  
  # Remove rows with zero variance
  isoforms <- remove_zero_variance(isoforms)
  
  return(isoforms)
}

# Function to process gene and isoform counts or CPM data for a given dataset
process_counts_or_cpm <- function(mat_dir, data_type, mart) {
  # Define file paths
  genes_file <- file.path(mat_dir, paste0("gene_", data_type, ".RDS"))
  isoforms_file <- file.path(mat_dir, paste0("isoform_", data_type, ".RDS"))
  
  # Process data
  genes <- filter_genes(genes_file, mart)
  isoforms <- filter_isoforms(isoforms_file, mart)
  
  return(list(genes = genes, isoforms = isoforms))
}

# Initialize BioMart
mart <- initialize_mart()

# Dataset directories
ROs_mat_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered"
FT_vs_RGC_mat_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered"

# Process both counts and CPM for ROs dataset
ROs_counts <- process_counts_or_cpm(ROs_mat_dir, "counts", mart)

# Process both counts and CPM for FT_vs_RGC dataset
FT_vs_RGC_counts <- process_counts_or_cpm(FT_vs_RGC_mat_dir, "counts", mart)

# Check number of rows
cat("ROs Counts - Genes:", nrow(ROs_counts$genes), " | Isoforms:", nrow(ROs_counts$isoforms), "\n")
cat("FT_vs_RGC Counts - Genes:", nrow(FT_vs_RGC_counts$genes), " | Isoforms:", nrow(FT_vs_RGC_counts$isoforms), "\n")


# 
# > cat("ROs Counts - Genes:", nrow(ROs_counts$genes), " | Isoforms:", nrow(ROs_counts$isoforms), "\n")
# ROs Counts - Genes: 19837  | Isoforms: 119690 
# > cat("FT_vs_RGC Counts - Genes:", nrow(FT_vs_RGC_counts$genes), " | Isoforms:", nrow(FT_vs_RGC_counts$isoforms), "\n")
# FT_vs_RGC Counts - Genes: 19411  | Isoforms: 94311 

####### Calculate size factors and convert to CPM in edgeR ########

library(edgeR)

# Function to compute size factors and convert counts to CPM
convert_to_cpm <- function(counts_matrix, group) {
  # Create a DGEList object
  dge <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, ]
  
  # Calculate normalization factors using TMM (Trimmed Mean of M-values)
  dge <- calcNormFactors(dge)
  
  # Convert to CPM (Counts Per Million)
  cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
  
  #return counts and CPM
  return(list(filtered_counts = counts_matrix[keep, ], cpm = cpm_matrix))
  
}

RO_group <- c("Stage_1", "Stage_1", "Stage_2", "Stage_2","Stage_2", "Stage_3", "Stage_3")
FT_vs_RGC_group <- c("FT", "FT", "RGC", "RGC")


# Convert gene and isoform counts to CPM
ROs_genes_result <- convert_to_cpm(ROs_counts$genes, group = RO_group)
ROs_isoforms_result <- convert_to_cpm(ROs_counts$isoforms, group = RO_group)

FT_vs_RGC_genes_result <- convert_to_cpm(FT_vs_RGC_counts$genes, group = FT_vs_RGC_group)
FT_vs_RGC_isoforms_result <- convert_to_cpm(FT_vs_RGC_counts$isoforms, group = FT_vs_RGC_group)

# Access filtered counts and CPM matrices
ROs_genes_filtered <- ROs_genes_result$filtered_counts
ROs_genes_cpm <- ROs_genes_result$cpm

FT_vs_RGC_genes_filtered <- FT_vs_RGC_genes_result$filtered_counts
FT_vs_RGC_genes_cpm <- FT_vs_RGC_genes_result$cpm

ROs_isoforms_filtered <- ROs_isoforms_result$filtered_counts
ROs_isoforms_cpm <- ROs_isoforms_result$cpm

FT_vs_RGC_isoforms_filtered <- FT_vs_RGC_isoforms_result$filtered_counts
FT_vs_RGC_isoforms_cpm <- FT_vs_RGC_isoforms_result$cpm




# Check dimensions
cat("Filtered ROs Genes:", dim(ROs_genes_filtered), "\n")
cat("CPM ROs Genes:", dim(ROs_genes_cpm), "\n")
cat("Filtered ROs Isoforms:", dim(ROs_isoforms_filtered), "\n")
cat("CPM ROs Isoforms:", dim(ROs_isoforms_cpm), "\n")

cat("Filtered FT_vs_RGC Genes:", dim(FT_vs_RGC_genes_filtered), "\n")
cat("CPM FT_vs_RGC Genes:", dim(FT_vs_RGC_genes_cpm), "\n")
cat("Filtered FT_vs_RGC Isoforms:", dim(FT_vs_RGC_isoforms_filtered), "\n")
cat("CPM FT_vs_RGC Isoforms:", dim(FT_vs_RGC_isoforms_cpm), "\n")

# 
# > cat("Filtered ROs Genes:", dim(ROs_genes_filtered), "\n")
# Filtered ROs Genes: 16971 7 
# > cat("CPM ROs Genes:", dim(ROs_genes_cpm), "\n")
# CPM ROs Genes: 16971 7 
# > cat("Filtered ROs Isoforms:", dim(ROs_isoforms_filtered), "\n")
# Filtered ROs Isoforms: 41290 7 
# > cat("CPM ROs Isoforms:", dim(ROs_isoforms_cpm), "\n")
# CPM ROs Isoforms: 41290 7 
# > cat("Filtered FT_vs_RGC Genes:", dim(FT_vs_RGC_genes_filtered), "\n")
# Filtered FT_vs_RGC Genes: 15101 4 
# > cat("CPM FT_vs_RGC Genes:", dim(FT_vs_RGC_genes_cpm), "\n")
# CPM FT_vs_RGC Genes: 15101 4 
# > cat("Filtered FT_vs_RGC Isoforms:", dim(FT_vs_RGC_isoforms_filtered), "\n")
# Filtered FT_vs_RGC Isoforms: 27163 4 
# > cat("CPM FT_vs_RGC Isoforms:", dim(FT_vs_RGC_isoforms_cpm), "\n")
# CPM FT_vs_RGC Isoforms: 27163 4 

### save counts and CPM

ROs_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype"
FT_vs_RGC_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype"

# Create output directories if they do not exist
dir.create(ROs_output_dir, showWarnings = FALSE)
dir.create(FT_vs_RGC_output_dir, showWarnings = FALSE)

# Save filtered counts and CPM matrices
saveRDS(ROs_genes_filtered, file.path(ROs_output_dir, "genes_counts.RDS"))
saveRDS(ROs_genes_cpm, file.path(ROs_output_dir, "genes_cpm.RDS"))

saveRDS(ROs_isoforms_filtered, file.path(ROs_output_dir, "isoform_counts.RDS"))
saveRDS(ROs_isoforms_cpm, file.path(ROs_output_dir, "isoform_cpm.RDS"))




