library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

sample <- commandArgs(trailingOnly = TRUE)

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/"
bam_dir <- paste0(bam_dir, sample, "_chromosome_level/")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"




for (chr in c(1:22, "X", "Y", "M")){ 
  se_output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/", sample, "/")
  print(sample)
  print(chr)
  se_quant_chr_sample <- bambu(reads = paste0(bam_dir, sample, "_", chr, ".bam"),
                     annotations = annotation,
                     genome = fa.file,
                     quant = FALSE, NDR = 1)
  
  writeBambuOutput(se_quant_chr_sample, 
                   path = se_output_dir,
                   prefix = sample)
  }

sessioninfo::session_info()
                          
