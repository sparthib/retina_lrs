# Load necessary libraries
library(plyranges)
library(Rsamtools)
library(dplyr)
library(tidyr)
library(stringr)

### 1. Load protein-coding multi-exonic genes
sample <- commandArgs(trailingOnly = TRUE)[1]q
print(sample)
alignment_dir <- '/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only'
alignment <- paste0(alignment_dir, '/', sample, '_primary_over_30_chr_only_sorted.bam')


### 2. Define the function to compute the number of exon-exon junctions covered by one read based on CIGAR string
compute_num_junction_per_read <- function(cigar_string) {
  m <- str_count(cigar_string, "N")
  return(m)
}

bamfile <- scanBam(BamFile(alignment))
# names(bamfile[[1]]) to print colnames in bamfile

nums <- c()  

### 3. Compute number of junctions in read
for (cig_string in bamfile[[1]]$cigar) {
  num <- compute_num_junction_per_read(cig_string)
  nums <- c(nums, num)
}

table(nums)
session_info::sessionInfo()

