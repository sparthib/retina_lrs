
# Load necessary libraries
library(rtracklayer)
library(dplyr)

# Define file paths
isoquant_GTF <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT",
                          "OUT.transcript_models.gtf")

method <- "Isoquant"
comparison <- "ROs"
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered")
counts <- file.path(matrix_dir, "isoform_counts.RDS") 

# Load counts matrix
counts <- readRDS(counts)

# Read GTF file
gtf <- import(isoquant_GTF)

# Extract isoform IDs from the counts matrix
isoform_ids <- rownames(counts)

# Filter GTF to retain only features corresponding to isoform IDs
filtered_gtf <- subset(gtf, gtf$transcript_id %in% isoform_ids)

# Save filtered GTF to a new file
filtered_gtf_file <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT", "ROs_filtered_isoforms.gtf")
export(filtered_gtf, con = filtered_gtf_file)

# Message to indicate success
cat("Filtered GTF file saved at:", filtered_gtf_file, "\n")

