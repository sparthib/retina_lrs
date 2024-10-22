library(biomaRt)
library(GenomicAlignments)
#get transcript_coverage function from transcriptome_coverage.R
source("/users/sparthib/retina_lrs/code/coverage/transcriptome_coverage.R")

# Define BAM file directory and sample name
bam_dir <- "/dcs04/hicks/data/sparthib/sg_nex_data/transcriptome_bam_files"
#get list of files in bam dir 
bam_files <- list.files(bam_dir, pattern = ".bam$", full.names = TRUE)
samples <- gsub(".bam", "", basename(bam_files))



# Read alignments from BAM file
output_dir <- "/dcs04/hicks/data/sparthib/sg_nex_data/coverage/alignment_files"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (sample_name in samples) {
  if (!file.exists(file.path(output_dir, paste0(sample_name, ".rds")))) {
    bam <- file.path(bam_dir, paste0(sample_name, ".bam"))
    aln <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
    
    # Save the alignments as an RDS file
    saveRDS(aln, file = file.path(output_dir, paste0(sample_name, ".rds")))
  } else {
    aln <- readRDS(file.path(output_dir, paste0(sample_name, ".rds")))
  }
}



#remove version num in seqnames
seqlevels(aln) <- gsub("\\..*", "", seqlevels(aln))
seqnames(aln) <- gsub("\\..*", "", seqnames(aln))

#get GC content of each transcript

process_sample_alignment <- function(aln, sample_name, output_dir, mart, annotation_cache = list(), length_bins = c(0, 1, 2, 5, 10, Inf)) {
  
  # Path to save the processed transcript info
  aln_rds_path <- file.path(output_dir, paste0(sample_name, ".rds"))
  
  # Check if transcript info has already been processed and saved
  if (file.exists(aln_rds_path)) {
    message(paste("Transcript info for", sample_name, "already exists, loading from file."))
    return(readRDS(aln_rds_path))
  }
  
  # Modify seqlevels/seqnames of alignment
  seqlevels(aln) <- gsub("\\..*", "", seqlevels(aln))
  seqnames(aln) <- gsub("\\..*", "", seqnames(aln))
  
  # Get sequence levels (transcript IDs)
  seq_ids <- seqlevels(aln)
  
  # Check if the seqlevels are already cached, and if not, fetch missing ones
  if (!all(seq_ids %in% names(annotation_cache))) {
    missing_ids <- setdiff(seq_ids, names(annotation_cache))
    
    annotLookup <- getBM(
      mart = mart,
      attributes = c("ensembl_transcript_id", "hgnc_symbol", 
                     "percentage_gene_gc_content", "transcript_biotype", 
                     "transcript_length"),
      filters = "ensembl_transcript_id",
      values = missing_ids,
      uniqueRows = TRUE
    )
    
    # Cache new transcript annotations
    annotation_cache <- c(annotation_cache, split(annotLookup, annotLookup$ensembl_transcript_id))
  }
  
  # Combine cached annotations into a single data frame for current seqlevels
  annotLookup <- do.call(rbind, annotation_cache[seq_ids])
  
  # Get isoforms with more than 10 occurrences
  isoform <- names(table(seqnames(aln)))[table(seqnames(aln)) > 10]
  
  # Get coverage info for transcripts
  transcript_info <- transcript_coverage(aln, isoform, length_bins)
  
  # Add isoform info and merge with annotations
  transcript_info$isoform <- rownames(transcript_info)
  transcript_info <- merge(transcript_info, annotLookup, by.x = "isoform", 
                           by.y = "ensembl_transcript_id", all.x = TRUE)
  
  # Save the processed transcript info
  saveRDS(transcript_info, file = aln_rds_path)
}


library(httr)
library(biomaRt)

# Set a longer global timeout (e.g., 20 seconds)
set_config(config(timeout = 20))

# Establish the Ensembl connection
us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)
# mart <- useDataset("hsapiens_gene_ensembl", us_mart)

annotation_cache <- list()

# Directory to save the processed RDS files
output_dir <- "/dcs04/hicks/data/sparthib/sg_nex_data/coverage/transcript_info"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Loop through samples and process each alignment
for (sample_name in samples) {
  input_dir <- "/dcs04/hicks/data/sparthib/sg_nex_data/coverage/alignment_files"
  aln <- readRDS(file.path(input_dir, paste0(sample_name, ".rds")))  # Load alignment
  output_dir <- "/dcs04/hicks/data/sparthib/sg_nex_data/coverage/transcript_info"
  transcript_info <- process_sample_alignment(aln, sample_name, output_dir, mart, annotation_cache)
  
  # Optionally: Do something with `transcript_info` here if needed
}
