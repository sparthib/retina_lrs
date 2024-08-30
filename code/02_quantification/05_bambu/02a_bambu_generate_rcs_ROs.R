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

#select RO samples 
rds_files <- rds_files[c(2,3,4,5,6,7,8)]

se   <- bambu(reads = rds_files,
              annotations = annotation,
              genome = fa.file,
              trackReads = TRUE)

dir.create("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation_track_reads/", showWarnings = FALSE)
writeBambuOutput(se, path = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation_track_reads/")

sessioninfo::session_info()
# library(rtracklayer)
# discovery_only <- readGFF("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/multisample_output.gtf")
# 
# plotBambu(discovery_only, type = "heatmap")
# 



