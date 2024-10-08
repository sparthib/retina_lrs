# Load necessary libraries
library(GenomicAlignments)
library(Rsamtools)
library(sessioninfo)
library(ggplot2)

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

# Define plot_coverage function
plot_coverage <- function(output_dir, transcript_info, sample, length_bin = NULL, gc_content = NULL, transcript_biotype = NULL) {
  
  # Filter based on length_bin if provided
  if (!is.null(length_bin)) {
    transcript_info <- transcript_info[transcript_info$length_bin == length_bin, ]
    length_bin_str <- as.character(length_bin)
    length_bin_str <- gsub("[\\(\\),\\[\\]]", "_", length_bin_str)  # Clean the string
  } else {
    length_bin_str <- "all_lengths"
  }
  
  # Filter based on gc_content if provided
  if (!is.null(gc_content)) {
    transcript_info <- transcript_info[transcript_info$gc_content_bin == gc_content, ]
    gc_content_str <- as.character(gc_content)
    gc_content_str <- gsub("[\\(\\),\\[\\]]", "_", gc_content_str)  # Clean the string
  } else {
    gc_content_str <- "all_gc_content"
  }
  
  # Filter based on transcript_biotype if provided
  if (!is.null(transcript_biotype)) {
    transcript_info <- transcript_info[transcript_info$transcript_biotype == transcript_biotype, ]
    biotype_str <- transcript_biotype
  } else {
    biotype_str <- "all_biotypes"
  }

  # Create filename for the PDF output
  pdf_file <- paste0(output_dir, sample, "_", length_bin_str, "_", gc_content_str, "_", biotype_str, "_coverage.pdf")
  
  # Ensure there are rows to plot
  if (nrow(transcript_info) == 0) {
    stop("No data left to plot after filtering.")
  }
  
  # Open PDF device for plotting
  pdf(pdf_file, width = 10, height = 5)
  
  # Ensure dev.off() is called even if there is an error
  on.exit(dev.off(), add = TRUE)
  
  # Iterate through each row (isoform) of transcript_info
  for (i in 1:nrow(transcript_info)) {
    isoform <- transcript_info$isoform[i]
    
    # Safeguard for potential data issues
    if (ncol(transcript_info) < 104) {
      stop("transcript_info does not have enough columns (expected at least 104 columns)")
    }
    
    # Extract and reverse coverage data
    coverage <- rev(as.numeric(unlist(transcript_info[i, 5:104])))
    
    # Check for invalid coverage values
    if (any(is.na(coverage))) {
      warning("NA values found in coverage data for isoform ", isoform)
      next  # Skip to the next isoform
    }
    
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
}
