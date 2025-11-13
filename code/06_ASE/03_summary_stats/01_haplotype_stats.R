# gets the percentage of haplotagged reads after whatshap haplotagging
code_dir <- Sys.getenv("retina_lrs_code")

txt <- file.path(code_dir, "/refs/heads/main/code/08_ASE/logs/haplotag_subset.txt")
file_content <- readLines(txt)

# Initialize empty vectors to store data
sample_names <- c()
total_alignments <- c()
haplotype1_reads <- c()
haplotype2_reads <- c()

# Parse the file
current_sample <- NULL

for (line in file_content) {
  # Extract sample name
  if (grepl("\\*\\*\\*\\* Haplotagging sample:", line)) {
    current_sample <- gsub(".*Haplotagging sample: (.*) \\*\\*\\*\\*", "\\1", line)
    current_sample <- trimws(current_sample)
  }
  
  # Extract total alignments processed
  if (grepl("Total alignments processed:", line)) {
    total_align <- gsub(".*Total alignments processed:\\s+(\\d+).*", "\\1", line)
    total_alignments <- c(total_alignments, as.numeric(total_align))
  }
  
  # Extract haplotype 1 reads
  if (grepl('Number of output reads haplotype 1:', line)) {
    hap1 <- gsub('.*Number of output reads haplotype 1:\\s+(\\d+).*', "\\1", line)
    haplotype1_reads <- c(haplotype1_reads, as.numeric(hap1))
  }
  
  # Extract haplotype 2 reads
  if (grepl('Number of output reads haplotype 2:', line)) {
    hap2 <- gsub('.*Number of output reads haplotype 2:\\s+(\\d+).*', "\\1", line)
    haplotype2_reads <- c(haplotype2_reads, as.numeric(hap2))
    sample_names <- c(sample_names, current_sample)
  }
}

# Create dataframe
df <- data.frame(
  sample_name = sample_names,
  total_alignments_processed = total_alignments,
  haplotype1_reads = haplotype1_reads,
  haplotype2_reads = haplotype2_reads
)

# View the result
print(df)

df$percentage_hap1 <- (df$haplotype1_reads / df$total_alignments_processed) * 100
df$percentage_hap2 <- (df$haplotype2_reads / df$total_alignments_processed) * 100
df$percent_hap <- df$percentage_hap1 + df$percentage_hap2
print(df)

# Save the dataframe to a CSV file
readr::write_tsv(df, file.path(code_dir,
                               "processed_data/ASE/bam_stats/haplotagging_summary_stats.tsv"))

