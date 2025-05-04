library(biomaRt)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
# Define file paths
gtf_file <- paste0(bambu_dir, "/extended_annotations.gtf")
gtf <- import(gtf_file)
#match isoforms to genes
isoforms_and_genes <- as.data.frame(gtf) |> 
  dplyr::select(transcript_id, gene_id) 
isoforms_and_genes$transcript_id <- gsub("\\..*", "", isoforms_and_genes$transcript_id)
isoforms_and_genes$gene_id <- gsub("\\..*", "", isoforms_and_genes$gene_id)


# Define directories
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype/"

# Desired column order
new_order <- c("H9_CRX_ROs_D45", "EP1_WT_ROs_D45", "EP1_WT_hRO_2", 
               "H9_BRN3B_hRO_2", "H9_CRX_hRO_2", "EP1_BRN3B_RO", 
               "H9_BRN3B_RO", "H9_hRGC_1", "H9_hRGC_2","H9_FT_1", "H9_FT_2")

# 1. Function to process data tables (column, row renaming and reordering columns )
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
process_table(file.path(bambu_dir, "counts_gene.txt"),
              "gene_counts", remove_gene_id = FALSE)

# 2. load processed data and remove isoforms only appearing in bambu 
isoform_counts <- readRDS(file.path(output_dir, 
                                    "RO_vs_RGC", "isoform_counts_RO_vs_RGC.RDS"))

gene_counts <- readRDS(file.path(output_dir, "RO_vs_RGC", "gene_counts_RO_vs_RGC.RDS"))
## remove version id from gene_counts
rownames(gene_counts) <- gsub("\\..*", "", rownames(gene_counts))

##### filter by common novel isoforms between bambu and isoquant ######
## common_isoforms only has the known gene related novel isoforms. 

common_isoforms <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)

common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                                      !grepl("^BambuGene", isoquant_gene_id))

colnames(common_isoforms)
common_isoforms$ref_gene <- gsub("\\..*", "", common_isoforms$ref_gene)
# 497


print(paste("pre filter", nrow(isoform_counts)))
# "pre filter 278023"

isoform_counts <- isoform_counts[
  (grepl("^Bambu", rownames(isoform_counts)) & rownames(isoform_counts) %in% common_isoforms$ref_id) |
    grepl("^ENST", rownames(isoform_counts)),
]


print(paste("post filter", nrow(isoform_counts)))
# "post filter 277401",  622 isoforms were labeled novel in Bambu but not found in isoquant. 


## only keep annotated genes in gene counts 
nrow(gene_counts)
#70341
gene_counts <- gene_counts[grepl("^ENSG", rownames(gene_counts)),]  
nrow(gene_counts)
# 70116

### 3. filter by biotype ###

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

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)  

genes <- filter_genes(gene_counts, mart)
nrow(genes)

novel_isoforms_from_protein_coding_genes <- common_isoforms$ref_id[
  common_isoforms$ref_gene %in% rownames(genes) 
]


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
  
  # Identify isoforms with "Bambu" in their rownames from PTC genes
  bambu_isoforms = novel_isoforms_from_protein_coding_genes
  
  all_isoforms <- c(protein_coding_isoforms, bambu_isoforms)
  
  isoforms <- isoforms[ isoforms$isoform_nums %in% all_isoforms, ]
  
  # Remove isoforms$isoform_nums column
  isoforms$isoform_nums <- NULL
  return(isoforms)
}

isoforms <- filter_isoforms(isoform_counts, mart)

nrow(genes)
# [1] 23216
nrow(isoforms)
# [1] 188091

groups <- c("Stage_1", "Stage_1", "Stage_2", "Stage_2",
            "Stage_2", "Stage_3", "Stage_3", "RGC", "RGC")


# Calculate size factors and convert to CPM in edgeR, and filter based on counts
filter_gene_counts <- function(counts_matrix, group ){ 
  min_counts <- 10
  dge <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(dge, group = group, min.count = min_counts)
  dge <- dge[keep, ]
  dge <- calcNormFactors(dge)
  # Convert to CPM (Counts Per Million)
  cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
  return(list(filtered_counts = counts_matrix[keep, ], cpm = cpm_matrix))
  
}

filter_isoform_counts <- function(isoform_gene_df, gene_counts, isoform_counts, group) { 
  
  genes_in_comparison <- rownames(gene_counts)
  #remove version numbers
  genes_in_comparison <- gsub("\\..*", "", genes_in_comparison)
  
  isoforms_in_comparison <- isoform_gene_df$transcript_id[ 
    isoform_gene_df$gene_id %in% genes_in_comparison
  ]
  rownames(isoform_counts) <- gsub("\\..*", "", rownames(isoform_counts))
  isoform_counts <- isoform_counts[rownames(isoform_counts) %in% isoforms_in_comparison, ]
  nrow(isoform_counts)
  
  dge <- DGEList(counts = isoform_counts)
  keep <- filterByExpr(dge, group = group, min.count = 2) #since we already filtered
  dge <- dge[keep, ]
  # Calculate normalization factors using TMM (Trimmed Mean of M-values)
  dge <- calcNormFactors(dge)
  # Convert to CPM (Counts Per Million)
  cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
  return(list(filtered_counts = isoform_counts[keep, ], cpm = cpm_matrix))
  
  
}

RO_RGC_genes_result <- filter_gene_counts(genes, groups)
RO_RGC_gene_counts <- RO_RGC_genes_result$filtered_counts
RO_RGC_genes_cpm <- RO_RGC_genes_result$cpm

nrow(RO_RGC_gene_counts)
#17112

ROs_isoforms_result <- filter_isoform_counts(isoforms_and_genes, 
                                             RO_RGC_gene_counts, 
                                             isoforms, 
                                             groups)

RO_RGC_isoform_counts <- ROs_isoforms_result$filtered_counts
RO_RGC_isoform_cpm <- ROs_isoforms_result$cpm
nrow(RO_RGC_isoform_counts)


# Save the filtered counts and CPM matrices
saveRDS(RO_RGC_gene_counts, file.path(output_dir,  "filtered_gene_counts.RDS"))
saveRDS(RO_RGC_isoform_counts, file.path(output_dir, "filtered_isoform_counts.RDS"))
saveRDS(RO_RGC_genes_cpm, file.path(output_dir, "filtered_gene_cpm.RDS"))
saveRDS(RO_RGC_isoform_cpm, file.path(output_dir, "filtered_isoform_cpm.RDS"))




