# Load necessary libraries
library(GenomicAlignments)
library(Rsamtools)
library(sessioninfo)
library(ggplot2)

# Define BAM file directory and sample name
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered"
sample <- commandArgs(trailingOnly = TRUE)[1]
bam <- file.path(bam_dir, paste0(sample, ".bam"))

# Read alignments from BAM file
aln <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))

# Get isoforms with more than 10 occurrences
isoform <- names(table(seqnames(aln)))[table(seqnames(aln)) > 10]
length(isoform)


length_bins = c(0, 1, 2, 5, 10, Inf)




# Define the transcript coverage function
transcript_coverage <- function(bam, isoform, length_bins, weight_fn = "read_counts") {
  if (!is.null(isoform)) {
    bam <- bam[GenomicAlignments::seqnames(bam) %in% isoform]
  }
  
  read_counts <- table(GenomicAlignments::seqnames(bam))
  transcript_names <- names(read_counts)
  
  transcript_info <- data.frame(
    tr_length = GenomeInfoDb::seqlengths(bam)[transcript_names],
    read_counts = as.data.frame(read_counts[transcript_names])$Freq
  )
  
  # Remove isoforms with zero read counts
  transcript_info <- transcript_info[transcript_info$read_counts != 0, ]
  transcript_info$length_bin <- cut(transcript_info$tr_length / 1000, length_bins)

  cover <- bam |>
    GenomicRanges::granges() |>
    GenomicRanges::coverage() |>
    sapply(function(x) {
      x[round(seq(1, length(x), length.out = 100), 0)] |>
        as.integer()
    }) |>
    subset(select = rownames(transcript_info)) |>
    t() |>
    as.data.frame()
  
  colnames(cover) <- paste0("coverage_", 1:100)
  # Define weight function
  weight_fn <- function(mat, read_counts) { read_counts }
  
  # Scale coverage by read counts
  cover <- cover / transcript_info$read_counts
  transcript_info <- cbind(transcript_info, cover)
  transcript_info$read_counts <- weight_fn(mat, transcript_info$read_counts)
  
  return(transcript_info)
}


transcript_info <- transcript_coverage(aln, isoform, length_bins)
nrow(transcript_info)


# Define output directory for plots
output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/"

# Define the plot coverage function
plot_coverage <- function(output_dir, transcript_info, sample, length_bin) {
  # Filter based on length_bin if provided
  if (!is.null(length_bin)) {
    transcript_info <- transcript_info[transcript_info$length_bin == length_bin, ]
  }
  
  # Create the output PDF
  pdf(paste0(output_dir, sample, "_", length_bin, "_coverage.pdf"), width = 10, height = 5)
  
  # Iterate through each row (isoform) of transcript_info
  for (i in 1:nrow(transcript_info)) {
    isoform <- rownames(transcript_info)[i]
    coverage <- rev(as.numeric(unlist(transcript_info[i, 4:103])))
    
    # Create data frame for plotting
    coverage_df <- data.frame(
      position = 1:100,
      coverage = coverage
    )
    
    # Generate the plot
    p <- ggplot(coverage_df, aes(x = position, y = coverage)) + 
      geom_line() + 
      labs(
        title = paste("Coverage of", isoform),
        x = "3' to 5' position",
        y = "Coverage Value"
      ) + 
      theme_minimal()
    
    # Print the plot to the PDF
    print(p)
  }
  
  # Close the PDF device
  dev.off()
}

# Generate coverage plots for each length bin
length_bins_to_plot <- c("(10,Inf]", "(5,10]", "(2,5]", "(1,2]", "(0,1]")
for (bin in length_bins_to_plot) {
  plot_coverage(output_dir, transcript_info, sample, bin)
}

# Print session info
sessionInfo()



