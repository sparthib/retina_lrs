library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/"
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
rds <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/track_reads/H9-FT_2/25d6d66a6f8e25_25d6d66a6f8e25.rds"

Sys.time()
se <- bambu(reads = rds,
            annotations = annotation,
            genome = fa.file,
            trackReads = TRUE)
Sys.time()

metadata(se)$readToTranscriptMaps[[1]]$equalMatches

#get entries where metadata(se)$readToTranscriptMaps[[1]]$equalMatches is not NULL
metadata(se)$readToTranscriptMaps[[1]]$equalMatches |>
  purrr::keep(~!is.null(.x))





