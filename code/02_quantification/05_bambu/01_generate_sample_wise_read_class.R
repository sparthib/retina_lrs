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
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/"
bam_dir <- paste0(bam_dir, sample, "_primary_over_30_chr_only_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"

output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/track_reads/", sample, "/")
print(sample)

se_read_class <- bambu(reads = paste0(bam_dir),
                     annotations = annotation,
                     genome = fa.file,
                     rcOutDir = output_dir,
                     lowMemory = TRUE,
                     trackReads = TRUE)


sessioninfo::session_info()
