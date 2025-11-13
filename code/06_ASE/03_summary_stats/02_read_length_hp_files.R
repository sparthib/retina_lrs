# get average read length for haplotagged sample files from whatshap
code_dir <- Sys.getenv("retina_lrs_code")
data_dir <- Sys.getenv("retina_lrs_dir")
library(Rsamtools)
library(dplyr)

samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "EP1-WT_ROs_D45", "EP1-BRN3B-RO", "EP1-WT_hRO_2", "H9-FT_1", 
             "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")

whatshap_out_dir <- file.path(data_dir,
                              "09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1")

results <- data.frame()

for (sample in samples) {
  input_h1_bam <- file.path(whatshap_out_dir, paste0(sample,"_h1.bam"))
  input_h2_bam <- file.path(whatshap_out_dir, paste0(sample,"_h2.bam"))
  
  # get read lengths
  h1_reads <- scanBam(input_h1_bam, param=ScanBamParam(what="qwidth"))[[1]]$qwidth
  h2_reads <- scanBam(input_h2_bam, param=ScanBamParam(what="qwidth"))[[1]]$qwidth
  
  # calculate means
  h1_mean <- mean(h1_reads, na.rm=TRUE)
  h2_mean <- mean(h2_reads, na.rm=TRUE)
  
  # add to results
  results <- rbind(results, 
                   data.frame(sample=sample, haplotype="h1", mean_read_length=h1_mean),
                   data.frame(sample=sample, haplotype="h2", mean_read_length=h2_mean))
}

print(results)

# Save results to a TSV file
output_dir <- file.path(code_dir, "processed_data/ASE/bam_stats")
file_name <- "mean_read_lengths_per_sample_haplotype.tsv"

output_path <- file.path(output_dir, file_name)
write.table(results, file=output_path, sep="\t", row.names=FALSE, quote=FALSE)
