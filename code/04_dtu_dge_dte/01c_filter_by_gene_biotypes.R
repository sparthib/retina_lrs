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

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)  

common_isoforms <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)

common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                                      !grepl("^BambuGene", isoquant_gene_id))


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
  
  return(genes)
}


# Function to load and filter isoform counts or CPM data
filter_isoforms <- function(counts_path, mart) {
  
  isoforms <- readRDS(counts_path)
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


RO_group <- c("Stage_1", "Stage_1", "Stage_2", "Stage_2","Stage_2", "Stage_3", "Stage_3")
FT_vs_RGC_group <- c("FT", "FT", "RGC", "RGC")



##### ROs filter genes and isoforms based on PTC #####

counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"

method <- "bambu"
comparison <- "ROs"
counts_dir <- file.path(counts_dir, method, comparison, "filtered")
isoform_counts_path <- file.path(counts_dir, "isoform_counts.RDS")
gene_counts_path <- file.path(counts_dir, "gene_counts.RDS")

### don't remove zero variance rows in genes or isoforms. 
genes <- filter_genes(gene_counts_path, mart)
nrow(genes)

novel_isoforms_from_protein_coding_genes <- common_isoforms$ref_id[
  common_isoforms$ref_gene %in% rownames(genes) 
]

isoforms <- filter_isoforms(isoform_counts_path, mart)

nrow(isoforms)

#### final gene and isoform counts ####
ROs_genes_result <- filter_gene_counts(genes, RO_group)
ROs_gene_counts <- ROs_genes_result$filtered_counts
ROs_genes_cpm <- ROs_genes_result$cpm

ROs_isoforms_result <- filter_isoform_counts(isoforms_and_genes, 
                                             ROs_gene_counts, 
                                             isoforms, 
                                             RO_group)


ROs_isoform_counts <- ROs_isoforms_result$filtered_counts
ROs_isoforms_cpm <- ROs_isoforms_result$cpm


nrow(ROs_isoforms_cpm)
nrow(ROs_genes_cpm)
nrow(ROs_isoform_counts)
nrow(ROs_gene_counts)

ROs_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype/"

saveRDS(ROs_isoforms_cpm, file.path(ROs_output_dir, "filtered_isoform_cpm.RDS"))
saveRDS(ROs_genes_cpm, file.path(ROs_output_dir, "filtered_gene_counts_cpm.RDS"))
saveRDS(ROs_isoform_counts, file.path(ROs_output_dir, "filtered_isoform_counts.RDS"))
saveRDS(ROs_gene_counts, file.path(ROs_output_dir, "filtered_gene_counts.RDS"))

##### FT vs RGC filter genes and isoforms based on PTC #####

method <- "bambu"
comparison <- "FT_vs_RGC"
counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"
counts_dir <- file.path(counts_dir, method, comparison, "filtered")
isoform_counts_path <- file.path(counts_dir, "isoform_counts.RDS")
gene_counts_path <- file.path(counts_dir, "gene_counts.RDS")

genes <- filter_genes(gene_counts_path, mart)
nrow(genes)

novel_isoforms_from_protein_coding_genes <- common_isoforms$ref_id[
  common_isoforms$ref_gene %in% rownames(genes) 
]

isoforms <- filter_isoforms(isoform_counts_path, mart)

nrow(isoforms)

FT_vs_RGC_gene_result <- filter_gene_counts(genes, FT_vs_RGC_group)
FT_vs_RGC_gene_counts <- FT_vs_RGC_gene_result$filtered_counts
FT_vs_RGC_gene_cpm <- FT_vs_RGC_gene_result$cpm

FT_vs_RGC_isoform_result <- filter_isoform_counts(isoforms_and_genes, 
                                                  FT_vs_RGC_gene_counts, 
                                                  isoforms, 
                                                  FT_vs_RGC_group)

FT_vs_RGC_isoform_counts <- FT_vs_RGC_isoform_result$filtered_counts
FT_vs_RGC_isoform_cpm <- FT_vs_RGC_isoform_result$cpm


nrow(FT_vs_RGC_isoform_counts)
nrow(FT_vs_RGC_gene_counts)

FT_vs_RGC_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype/"


saveRDS(FT_vs_RGC_isoform_counts, file.path(FT_vs_RGC_output_dir, "filtered_isoform_counts.RDS"))
saveRDS(FT_vs_RGC_gene_counts, file.path(FT_vs_RGC_output_dir, "filtered_gene_counts.RDS"))

saveRDS(FT_vs_RGC_isoform_cpm, file.path(FT_vs_RGC_output_dir, "filtered_isoform_cpm.RDS"))
saveRDS(FT_vs_RGC_gene_cpm, file.path(FT_vs_RGC_output_dir, "filtered_gene_counts_cpm.RDS"))




###### filter gtf based on counts ######
# Function to filter GTF based on transcripts in the counts matrix
filter_gtf_by_counts <- function(gtf_file, counts_matrix, output_gtf) {
  # Load GTF
  gtf <- import(gtf_file)
  #remove version numbers
  gtf$transcript_id <- gsub("\\..*", "", gtf$transcript_id)
  
  # Get transcript IDs from counts matrix row names
  counts_transcripts <- rownames(counts_matrix)
  counts_transcripts <- gsub("\\..*", "", counts_transcripts)  
  
  # Filter GTF to keep only matching transcript IDs
  gtf_filtered <- gtf[gtf$transcript_id %in% counts_transcripts]

  # Save the filtered GTF file
  export(gtf_filtered, output_gtf)
  
  return(gtf_filtered)
}


# Apply function to your isoform counts matrix
ROs_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype"
ROs_counts <- readRDS(file.path(ROs_dir, "filtered_isoform_counts.RDS"))
nrow(ROs_counts)

ROs_output_gtf <- paste0(bambu_dir,
                          "/ROs_protein_coding_annotations.gtf")
ROs_filtered_gtf <- filter_gtf_by_counts(gtf_file, 
                                          ROs_counts,
                                          ROs_output_gtf)

FT_vs_RGC_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype"
FT_vs_RGC_counts <- readRDS(file.path(FT_vs_RGC_dir, "filtered_isoform_counts.RDS"))
nrow(FT_vs_RGC_counts)

FT_vs_RGC_output_gtf <- paste0(bambu_dir,
                                "/FT_vs_RGC_protein_coding_annotations.gtf")
FT_vs_RGC_filtered_gtf <- filter_gtf_by_counts(gtf_file, 
                                                FT_vs_RGC_counts,
                                                FT_vs_RGC_output_gtf)

RO_vs_RGC_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype"
RO_vs_RGC_counts <- readRDS(file.path(RO_vs_RGC_dir, "filtered_isoform_counts.RDS"))
nrow(RO_vs_RGC_counts)

RO_vs_RGC_output_gtf <- paste0(bambu_dir,
                               "/RO_vs_RGC_protein_coding_annotations.gtf")

RO_vs_RGC_filtered_gtf <- filter_gtf_by_counts(gtf_file, 
                                               RO_vs_RGC_counts,
                                               RO_vs_RGC_output_gtf)



