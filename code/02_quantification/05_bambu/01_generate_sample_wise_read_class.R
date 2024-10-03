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
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality"
bam_dir <- paste0(bam_dir, sample, "_primary_over_30_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"

output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/", sample, "/")

if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
print(sample)

se_read_class <- bambu(reads = paste0(bam_dir),
                     annotations = annotation,
                     genome = fa.file,
                     rcOutDir = output_dir,
                     lowMemory = TRUE,
                     trackReads = TRUE)


sessioninfo::session_info()


# class(metadata(se)$readToTranscriptMaps)
# colnames(metadata(se)$readToTranscriptMaps[[1]])
# 
# metadata(se)$readToTranscriptMaps[1][1]
# length(metadata(se)$readToTranscriptMaps[[7]])

# temp <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/track_reads/EP1-WT_hRO_2/25dd8655e094e7_25dd8655e094e7.rds")
