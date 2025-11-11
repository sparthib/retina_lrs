library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

annotation <- readRDS(file.path(data_dir, "/06_quantification/bambu/annotations.rds"))
bam_dir <- file.path(data_dir, "05_bams/genome/primary_assembly/high_quality")
fa.file <- file.path(ref_dir, "genome/GENCODE/primary_assembly/release_46_primary_genome.fa")

root_dir <- file.path(data_dir, "06_quantification/bambu/rc_output")
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

dir.create(file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/"), showWarnings = FALSE)
writeBambuOutput(se, path = file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/"))
saveRDS(se, file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/se.rds"))
sessioninfo::session_info()
