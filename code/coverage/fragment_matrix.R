# Load necessary libraries
library(GenomicAlignments)
library(Rsamtools)
library(sessioninfo)
library(ggplot2)
library(Biostrings)  # For efficient sequence manipulation
library(dplyr)

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered"
# sample <- commandArgs(trailingOnly = TRUE)[1]
sample <- "H9-FT_1"
bam <- file.path(bam_dir, paste0(sample, ".bam"))

param <- ScanBamParam(what = c("seq"), 
                      mapqFilter = 5)

aln <- readGAlignments(bam, param = param)

# Extract the read sequences
sequences <- mcols(aln)$seq

# Calculate GC content percentage for each sequence
gc_content <- letterFrequency(sequences, letters = c("G", "C"), as.prob = TRUE)
gc_percent <- rowSums(gc_content) * 100  # Convert to percentage

# Add GC content as a new column to the GAlignments metadata
mcols(aln)$GC_percent <- gc_percent

# Check the result
head(mcols(aln))




# output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/GAlignments/"
# aln <- readRDS(file.path(output_dir, paste0(sample, ".rds")))

# Check the GAlignments object


seqlevels(aln) <- gsub("\\..*", "", seqlevels(aln))
seqnames(aln) <- gsub("\\..*", "", seqnames(aln))


# get seqnames, width, GC_percent
df <- data.frame(seqnames = seqnames(aln), 
                 width = width(aln), 
                 GC_percent = mcols(aln)$GC_percent)

max(df$width)

# bin width into 10bp bins
df$read_length_bin <- cut(df$width, breaks = seq(0, max(df$width), 10), include.lowest = TRUE)

#round GC percent to nearest integer
df$GC_percent <- round(df$GC_percent)

#group by seqnames,  GC_percent and read_length_bin and get count for each group 
fragment_table <- df |>  group_by(seqnames, GC_percent, read_length_bin) |>
  summarise(count = n())



