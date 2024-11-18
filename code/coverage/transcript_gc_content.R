#compute transcript GC content from a transcriptome fasta file 
library(Biostrings)

fa_dir <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
transcript_fa <- "release_46_all_transcripts_short_header.fa" 
transcript_fa <- file.path(fa_dir, transcript_fa)

transcripts <- readDNAStringSet(transcript_fa)

length(transcripts)
# Calculate GC content for each transcript
gc_content <- letterFrequency(transcripts, letters = c("G", "C"), as.prob = TRUE)
gc_percent <- rowSums(gc_content) * 100  # Convert to percentage

# get length of transcript 
transcript_length <- width(transcripts)

# Create a data frame with transcript IDs and GC content
transcript_gc_df <- data.frame(
  isoform_id = names(transcripts),
  GC_percent = gc_percent,
  transcript_length = transcript_length
)

#remove version number from transcript id
transcript_gc_df$isoform_id <- gsub("\\..*", "", transcript_gc_df$isoform_id)

# Display the result
saveRDS(transcript_gc_df, file = paste0(fa_dir, "/transcript_meta_info.rds"))
