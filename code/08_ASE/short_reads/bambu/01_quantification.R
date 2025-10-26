library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/extended_annotations.gtf"
# annotations <- prepareAnnotations(gtf_file)
# saveRDS(annotations, file.path(outdir, "annotations.rds"))

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/annotations.rds")

fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"

sample <- commandArgs(trailingOnly = TRUE)

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1"
h1_bam <- file.path(bam_dir, paste0(sample, "_h1.bam"))
h2_bam <- file.path(bam_dir, paste0(sample, "_h2.bam"))

se.quantOnly_h1 <- bambu(reads = h1_bam, annotations = annotation, genome = fa.file, discovery = FALSE)
se.quantOnly_h2 <- bambu(reads = h2_bam, annotations = annotation, genome = fa.file, discovery = FALSE)

outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/H9_EP1_bambu"
writeBambuOutput(se.quantOnly_h1, path = file.path(outdir, sample, "hp1"))
writeBambuOutput(se.quantOnly_h2, path = file.path(outdir, sample, "hp2"))

sessioninfo::session_info()
