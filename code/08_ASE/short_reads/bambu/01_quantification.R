library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"
outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/H9_EP1_bambu"

sample <- commandArgs(trailingOnly = TRUE)

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1"
h1_bam <- file.path(bam_dir, paste0(sample, "_hp1.bam"))
h2_bam <- file.path(bam_dir, paste0(sample, "_hp2.bam"))

se.quantOnly_h1 <- bambu(reads = h1_bam, annotations = gtf.file, genome = fa.file, discovery = FALSE)
se.quantOnly_h2 <- bambu(reads = h2_bam, annotations = gtf.file, genome = fa.file, discovery = FALSE)
writeBambuOutput(se.quantOnly_h1, path = file.path(outdir, sample, "hp1"))
writeBambuOutput(se.quantOnly_h2, path = file.path(outdir, sample, "hp2"))

sessioninfo::session_info()
