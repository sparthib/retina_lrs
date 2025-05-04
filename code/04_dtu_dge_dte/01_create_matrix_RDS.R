# Define directories
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"

# Desired column order
new_order <- c("H9_CRX_ROs_D45", "EP1_WT_ROs_D45", "EP1_WT_hRO_2", 
               "H9_BRN3B_hRO_2", "H9_CRX_hRO_2", "EP1_BRN3B_RO", 
               "H9_BRN3B_RO", "H9_FT_1", "H9_FT_2", 
               "H9_hRGC_1", "H9_hRGC_2")

# Function to process data tables
process_table <- function(filepath, output_subdir, file_prefix, remove_gene_id = FALSE) {
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
  ROs <- data[, 1:7]
  FT_vs_RGC <- data[, 8:11]
  
  # Save subsets
  saveRDS(ROs, file.path(output_dir, output_subdir, "ROs", paste0(file_prefix, "_ROs.RDS")))
  saveRDS(FT_vs_RGC, file.path(output_dir, output_subdir, "FT_vs_RGC", paste0(file_prefix, "_FT_vs_RGC.RDS")))
}

# Bambu processing
process_table(file.path(bambu_dir, "counts_gene.txt"), "bambu", "gene_counts")
process_table(file.path(bambu_dir, "counts_transcript.txt"), "bambu",
              "isoform_counts", remove_gene_id = TRUE)


### load and check nrow 
# Load the RDS files
bambu_RO_genes <- readRDS(file.path(output_dir, "bambu", "ROs", "gene_counts_ROs.RDS"))
bambu_FT_vs_RGC_genes <- readRDS(file.path(output_dir, "bambu", "FT_vs_RGC", "gene_counts_FT_vs_RGC.RDS"))

bambu_RO_isoforms <- readRDS(file.path(output_dir, "bambu", "ROs", "isoform_counts_ROs.RDS"))
bambu_FT_vs_RGC_isoforms <- readRDS(file.path(output_dir, "bambu", "FT_vs_RGC", "isoform_counts_FT_vs_RGC.RDS"))

nrow(bambu_RO_isoforms)
nrow(bambu_FT_vs_RGC_isoforms)





