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
  rownames(data) <- sub("\\..*", "", rownames(data))
  
  # Split into subsets
  ROs <- data[, 1:7]
  FT_vs_RGC <- data[, 8:11]
  
  # Save subsets
  saveRDS(ROs, file.path(output_dir, output_subdir, "ROs", paste0(file_prefix, "_ROs.RDS")))
  saveRDS(FT_vs_RGC, file.path(output_dir, output_subdir, "FT_vs_RGC", paste0(file_prefix, "_FT_vs_RGC.RDS")))
}

# Bambu processing
process_table(file.path(bambu_dir, "counts_gene.txt"), "bambu", "gene_counts")
process_table(file.path(bambu_dir, "counts_transcript.txt"), "bambu", "isoform_counts", remove_gene_id = TRUE)
process_table(file.path(bambu_dir, "CPM_transcript.txt"), "bambu", "isoform_cpm", remove_gene_id = TRUE)

# Isoquant processing
process_table(file.path(isoquant_dir, "OUT.gene_grouped_counts.tsv"), "isoquant", "gene_counts")
process_table(file.path(isoquant_dir, "OUT.gene_grouped_tpm.tsv"), "isoquant", "gene_cpm")




isoquant_isoform_counts <- read.table(file.path(isoquant_dir,
                                                "OUT.transcript_model_grouped_counts.tsv"),
                                      header = TRUE,
                                      row.names = 1, sep = "\t", comment.char = "")
colnames(isoquant_isoform_counts) <- gsub("_primary_over_30_sorted", "", 
                                          colnames(isoquant_isoform_counts))
colnames(isoquant_isoform_counts) <- gsub("\\.", "_", colnames(isoquant_isoform_counts))

# Rearrange the columns in the data frame
isoquant_isoform_counts <- isoquant_isoform_counts[, new_order]

##keep only ROs
RO_isoquant_isoform_counts <- isoquant_isoform_counts[,1:7]
FT_vs_RGC_isoquant_isoform_counts <- isoquant_isoform_counts[,8:11]

saveRDS(RO_isoquant_isoform_counts, 
        file.path(output_dir, "isoquant", "ROs", "isoform_counts.RDS"))
saveRDS(FT_vs_RGC_isoquant_isoform_counts,
        file.path(output_dir, "isoquant", "FT_vs_RGC", "isoform_counts.RDS"))


###isoquant isoform cpm 

isoquant_isoform_cpm <- read.table(file.path(isoquant_dir, "OUT.transcript_model_grouped_tpm.tsv"), header = TRUE,
                                   row.names = 1, sep = "\t", comment.char = "")
colnames(isoquant_isoform_cpm) <- gsub("_primary_over_30_sorted", "", 
                                       colnames(isoquant_isoform_cpm))
colnames(isoquant_isoform_cpm) <- gsub("\\.", "_", colnames(isoquant_isoform_cpm))

# Rearrange the columns in the data frame
isoquant_isoform_cpm <- isoquant_isoform_cpm[, new_order]

##keep only ROs
RO_isoquant_isoform_cpm <- isoquant_isoform_cpm[,1:7]
FT_vs_RGC_isoquant_isoform_cpm <- isoquant_isoform_cpm[,8:11]

saveRDS(RO_isoquant_isoform_cpm, 
        file.path(output_dir, "isoquant", "ROs", "isoform_cpm.RDS"))
saveRDS(FT_vs_RGC_isoquant_isoform_cpm, 
        file.path(output_dir, "isoquant", "FT_vs_RGC", "isoform_cpm.RDS"))







