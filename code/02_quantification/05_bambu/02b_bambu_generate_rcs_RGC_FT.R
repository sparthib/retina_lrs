library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)


annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/"
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"

root_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/rc_output"
files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
rds_files <- files[grep("\\.rds$", files)]

#select H9 RGC FT samples 
rds_files <- rds_files[c(9,10,11,12)]

se <- bambu(reads = rds_files,
              annotations = annotation,
              genome = fa.file, 
            trackReads = TRUE)

dir.create("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation_track_reads/", showWarnings = FALSE)
writeBambuOutput(se, path = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation_track_reads/")
saveRDS(se, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation_track_reads/se.rds")
sessioninfo::session_info()



