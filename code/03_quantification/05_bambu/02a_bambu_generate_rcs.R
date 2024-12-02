library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)


annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality"
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"

root_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/rc_output"
files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
rds_files <- files[grep("\\.rds$", files)]

#select RO samples 
rds_files <- rds_files[c(2,3,4,5,6,7,8,9,10,11,12)]

se   <- bambu(reads = rds_files,
              annotations = annotation,
              genome = fa.file,
              trackReads = TRUE,
              quant = TRUE,
              discovery = TRUE)

dir.create("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/", showWarnings = FALSE)
writeBambuOutput(se, path = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/")
saveRDS(se, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/se.rds")
sessioninfo::session_info()
