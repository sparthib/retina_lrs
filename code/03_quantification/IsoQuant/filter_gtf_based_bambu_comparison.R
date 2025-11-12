
# Load necessary libraries
library(rtracklayer)
library(dplyr)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

# Define file paths
isoquant_GTF <- file.path(data_dir ,"/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.extended_annotation.gtf")


method <- "Isoquant"
comparison <- "ROs"
matrix_dir <- file.path(data_dir,"06_quantification/counts_matrices/",
                        method, comparison, "filtered")
counts <- file.path(matrix_dir, "isoform_counts.RDS") 

# Load counts matrix
counts <- readRDS(counts)

# Read GTF file
gtf <- import(isoquant_GTF)
gtf <- as.data.frame(gtf)
gtf |>
  filter(!str_starts(gene_id, "ENSG")) |>
  distinct(gene_id) |> nrow()
# Extract isoform IDs from the counts matrix
isoform_ids <- rownames(counts)

# Filter GTF to retain only features corresponding to isoform IDs
filtered_gtf <- subset(gtf, gtf$transcript_id %in% isoform_ids)

# Save filtered GTF to a new file
filtered_gtf_file <- file.path(data_dir,"06_quantification/isoquant/high_quality/all_samples/OUT", "ROs_filtered_isoforms.gtf")
export(filtered_gtf, con = filtered_gtf_file)

# Message to indicate success
cat("Filtered GTF file saved at:", filtered_gtf_file, "\n")

