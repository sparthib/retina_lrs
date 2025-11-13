library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

code_dir <- Sys.getenv("retina_lrs_code")
data_dir <- Sys.getenv("retina_lrs_dir")
ref_dir <- Sys.getenv("references_dir")

# gtf_file <- file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/extended_annotations.gtf")
# annotations <- prepareAnnotations(gtf_file)
# saveRDS(annotations, file.path(outdir, "annotations.rds"))

annotation <- readRDS(file.path(data_dir, "06_quantification/bambu/primary_assembly/annotations.rds"))

fa.file <- file.path(ref_dir,"genome/GENCODE/primary_assembly/release_46_primary_genome.fa")

sample <- commandArgs(trailingOnly = TRUE)

bam_dir <- file.path(data_dir,"09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1")
h1_bam <- file.path(bam_dir, paste0(sample, "_h1.bam"))
h2_bam <- file.path(bam_dir, paste0(sample, "_h2.bam"))

se.quantOnly_h1 <- bambu(reads = h1_bam, annotations = annotation, genome = fa.file, discovery = FALSE)
se.quantOnly_h2 <- bambu(reads = h2_bam, annotations = annotation, genome = fa.file, discovery = FALSE)

outdir <- file.path(data_dir,"09_ASE/H9_DNA_Seq_data/H9_EP1_bambu")
writeBambuOutput(se.quantOnly_h1, path = file.path(outdir, sample, "hp1"))
writeBambuOutput(se.quantOnly_h2, path = file.path(outdir, sample, "hp2"))

sessioninfo::session_info()
