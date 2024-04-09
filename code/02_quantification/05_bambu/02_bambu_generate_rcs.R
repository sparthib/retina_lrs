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
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"


output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/", sample, "/")
print(sample)
se_quant_sample <- bambu(reads = paste0(bam_dir, sample, "_sorted.bam"),
                               annotations = annotation,
                               genome = fa.file,
                          discovery = FALSE, 
                         quant = FALSE, 
                         rcOutDir = output_dir)

sessioninfo::session_info()

