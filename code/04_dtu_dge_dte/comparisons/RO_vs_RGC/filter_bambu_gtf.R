library(stringr)
library(dplyr)
library(rtracklayer)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")

input_gtf_dir <- file.path(data_dir,"06_quantification/bambu/all_samples_extended_annotation_track_reads")
input_gtf <- file.path(input_gtf_dir, "extended_annotations.gtf")

#load gtf
gtf <- rtracklayer::import(input_gtf)
head(gtf)


# 3. Strip version numbers from transcript_id (keep only the part before the dot)
transcript_ids <- mcols(gtf)$transcript_id
transcript_ids_base <- str_replace(transcript_ids, "\\..*", "")  # Remove version




matrix_dir <- file.path(data_dir, "06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype")

counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)

rownames(counts)

# Match base transcript IDs to rownames of counts matrix
keep <- transcript_ids_base %in% rownames(counts)

# Subset the GTF (but keep versioned transcript_id in result)
gtf_filtered <- gtf[keep]


export(gtf_filtered, file.path(input_gtf_dir, "RO_RGC_protein_coding.gtf"), format = "gtf")

