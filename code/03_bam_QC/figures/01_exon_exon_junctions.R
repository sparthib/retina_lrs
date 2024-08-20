# Load necessary libraries
library(plyranges)
library(Rsamtools)
library(rtracklayer)
library(dplyr)
library(tidyr)

### 1. Load protein-coding multi-exonic genes
gtf_fn <- '/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf'
gtf <- readGFF(gtf_fn) %>% as_tibble()

# Filter genes of interest (multi-exonic protein-coding genes)
# curr_gene <- ''
# goi <- c()
# curr_num_exon <- 0
# 
# for (i in 1:nrow(gtf)) {
#   if (gtf[i, "type"] == "gene" & gtf[i, "gene_type"] == "protein_coding") {
#     if (curr_num_exon > 1) {
#       goi <- c(goi, curr_gene)
#     }
#     curr_gene <- gtf[i, "gene_id"]
#     curr_num_exon <- 0
#   }
#   
#   if (gtf[i, "type"] == "exon" & gtf[i, "gene_id"] == curr_gene) {
#     curr_num_exon <- curr_num_exon + 1
#   }
# }
# 
# genic_gtf <- gtf %>% filter(type == "gene" & gene_id %in% goi)

# saveRDS(genic_gtf, '/dcs04/hicks/data/sparthib/references/genome/GENCODE/genic_gtf.rds')
# genic_gtf <- readRDS('/dcs04/hicks/data/sparthib/references/genome/GENCODE/genic_gtf.rds')
# 
# #get all rows that have chr prefix in seqid 
# #chrM is missing as we're only looking at protein-coding
# genic_gtf <- genic_gtf[grep("chr", genic_gtf$seqid), ]
# 
# unique(genic_gtf$seqid)
### 2. Get BAM file name
sample <- commandArgs(trailingOnly = TRUE)[1]
sample <- "H9-FT_1"
print(sample)
alignment_dir <- '/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only'
alignment <- paste0(alignment_dir, '/', sample, '_primary_over_30_chr_only_sorted.bam')


### 3. Define the function to compute the number of exon-exon junctions covered by one read based on CIGAR string
compute_num_junction_per_read <- function(cigar_string) {
  m <- str_count(cigar_string, "N")
  return(m)
}


bamfile <- scanBam(BamFile(alignment))
colnames(bamfile)

(bamfile[[1]]$cigar)[1] 

nums <- c()  

for (cig_string in bamfile[[1]]$cigar) {
  num <- compute_num_junction_per_read(cig_string)
  nums <- c(nums, num)
}

table(nums)



session_info::sessionInfo()

