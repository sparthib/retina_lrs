library(biomaRt)

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)  

# Function to remove rows with zero variance
remove_zero_variance <- function(mat) {
  mat[rowSums(mat != mat[, 1]) > 0, , drop = FALSE]  # Keep rows with at least one differing value
}

# Function to load and filter gene counts or CPM data
filter_genes <- function(counts_file, mart) {
  genes <- readRDS(counts_file)
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
  # Remove rows with zero variance
  genes <- remove_zero_variance(genes)
  
  return(genes)
}

# Function to load and filter isoform counts or CPM data
filter_isoforms <- function(counts_file, mart) {
  isoforms <- readRDS(counts_file)
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

# Calculate size factors and convert to CPM in edgeR
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
  
  return(list(filtered_counts = counts_matrix[keep, ], cpm = cpm_matrix))
}

RO_group <- c("Stage_1", "Stage_1", "Stage_2", "Stage_2","Stage_2", "Stage_3", "Stage_3")
FT_vs_RGC_group <- c("FT", "FT", "RGC", "RGC")

# Convert gene and isoform counts to CPM
ROs_genes_result <- convert_to_cpm(ROs_counts$genes, group = RO_group)
ROs_isoforms_result <- convert_to_cpm(ROs_counts$isoforms, group = RO_group)

FT_vs_RGC_genes_result <- convert_to_cpm(FT_vs_RGC_counts$genes, group = FT_vs_RGC_group)
FT_vs_RGC_isoforms_result <- convert_to_cpm(FT_vs_RGC_counts$isoforms, group = FT_vs_RGC_group)

# Save filtered counts and CPM matrices
dirs <- list(
  ROs = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype",
  FT_vs_RGC = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype"
)

lapply(dirs, dir.create, showWarnings = FALSE)

saveRDS(ROs_genes_result$filtered_counts, file.path(dirs$ROs, "genes_counts.RDS"))
saveRDS(ROs_genes_result$cpm, file.path(dirs$ROs, "genes_cpm.RDS"))

saveRDS(ROs_isoforms_result$filtered_counts, file.path(dirs$ROs, "isoform_counts.RDS"))
saveRDS(ROs_isoforms_result$cpm, file.path(dirs$ROs, "isoform_cpm.RDS"))

saveRDS(FT_vs_RGC_genes_result$filtered_counts, file.path(dirs$FT_vs_RGC, "genes_counts.RDS"))
saveRDS(FT_vs_RGC_genes_result$cpm, file.path(dirs$FT_vs_RGC, "genes_cpm.RDS"))

saveRDS(FT_vs_RGC_isoforms_result$filtered_counts, file.path(dirs$FT_vs_RGC, "isoform_counts.RDS"))
saveRDS(FT_vs_RGC_isoforms_result$cpm, file.path(dirs$FT_vs_RGC, "isoform_cpm.RDS"))

###### 

library(GenomicRanges)
library(rtracklayer)

# Function to filter GTF based on transcripts in the counts matrix
filter_gtf_by_counts <- function(gtf_file, counts_matrix, output_gtf) {
  # Load GTF
  gtf <- import(gtf_file)
  
  # Extract transcript IDs from the GTF file
  gtf_transcripts <- gtf$transcript_id
  
  # Get transcript IDs from counts matrix row names
  counts_transcripts <- rownames(counts_matrix)
  
  # Filter GTF to keep only matching transcript IDs
  gtf_filtered <- gtf[gtf_transcripts %in% counts_transcripts]
  
  # Save the filtered GTF file
  export(gtf_filtered, output_gtf)
  
  return(gtf_filtered)
}

bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"

# Define file paths
gtf_file <- paste0(bambu_dir, "/extended_annotations.gtf")
ROs_output_gtf <- paste0(bambu_dir, "/ROs_protein_coding_annotations.gtf")
FT_vs_RGC_output_gtf <- paste0(bambu_dir, "/FT_vs_RGC_protein_coding_annotations.gtf")

# Apply function to your isoform counts matrix
ROs_filtered_gtf <- filter_gtf_by_counts(gtf_file, ROs_isoforms_result$filtered_counts, ROs_output_gtf)
FT_vs_RGC_filtered_gtf <- filter_gtf_by_counts(gtf_file, FT_vs_RGC_isoforms_result$filtered_counts, FT_vs_RGC_output_gtf)



cat("Filtered GTF saved to:", output_gtf, "\n")


