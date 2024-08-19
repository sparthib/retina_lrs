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
genic_gtf <- readRDS('/dcs04/hicks/data/sparthib/references/genome/GENCODE/genic_gtf.rds')
### 2. Get BAM file names
sample <- commandArgs(trailingOnly = TRUE)[1]
alignment_dir <- '/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only'
alignment <- paste0(alignment_dir, '/', sample, '_primary_over_30_chr_only_sorted.bam')


### 3. Define the function to compute the number of exon-exon junctions covered by one read based on CIGAR string
compute_num_junction_per_read <- function(cigar_string) {
  m <- str_count(cigar_string, "N")
  return(m)
}

### 4. Calculate junction coverage for each BAM file
max_junction <- 11 
df_list <- list()

nums <- c()  
for (i in 1:nrow(genic_gtf)) {
  chr <- as.character(genic_gtf[i, "seqid"])
  start <- as.integer(genic_gtf[i, "start"])
  end <- as.integer(genic_gtf[i, "end"])
  
  reads <- scanBam(BamFile(alignment), param = ScanBamParam(which = GRanges(chr, IRanges(start, end))))
  
  for (read in reads[[1]]$cigar) {
    num <- compute_num_junction_per_read(read)
    nums <- c(nums, num)
  }
}
tmp_counter <- table(factor(nums, levels = 0:(max_junction - 1)))
df <- as.numeric(tmp_counter)
df 

per_df <- sweep(df, 1, rowSums(df), FUN = "/") 

per_df 


### 7. Save summary to CSV
output_dir <- '/users/sparthib/retina_lrs/processed_data/bam_qc/'
write.csv(per_df, paste0(output_dir, sample, "_exon-exon_junctions_perdf.csv"),
          row.names = FALSE)
write.csv(df, paste0(output_dir, sample, "_exon-exon_junctions_df.csv"),
          row.names = FALSE)


session_info::sessionInfo()


### 5. Calculate summary statistics for two types of sequencing
# ont_per_df <- per_df[1:9, ] %>% t() %>% as.data.frame()
# ont_per_df$mean <- rowMeans(ont_per_df)
# ont_per_df$sem <- apply(ont_per_df, 1, sd) / sqrt(ncol(ont_per_df))
# 
# rna_per_df <- per_df[10:nrow(per_df), ] %>% t() %>% as.data.frame()
# rna_per_df$mean <- rowMeans(rna_per_df)
# rna_per_df$sem <- apply(rna_per_df, 1, sd) / sqrt(ncol(rna_per_df))
# 
# ### 6. Prepare summary DataFrame
# summary_df <- data.frame(
#   group = rep(c("long_read", "short_read"), each = max_junction),
#   num = rep(0:(max_junction - 1), 2),
#   mean = c(ont_per_df$mean, rna_per_df$mean),
#   sem = c(ont_per_df$sem, rna_per_df$sem)
# )

