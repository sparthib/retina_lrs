library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

# gtf.file <- file.path(ref_dir, "genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf")
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, file.path(data_dir, 
#   06_quantification/bambu/primary_assembly/annotations.rds"))

sample <- commandArgs(trailingOnly = TRUE)

annotation <- readRDS(file.path(data_dir, "06_quantification/bambu/primary_assembly/annotations.rds"))
bam_dir <- readRDS(file.path(data_dir, "05_bams/genome/primary_assembly/high_quality/"))
bam_file <- paste0(bam_dir, sample, "_primary_over_30_chr_only_sorted.bam")
fa.file <- file.path(ref_dir, "genome/GENCODE/primary_assembly/release_46_primary_genome.fa")

output_dir <- file.path(data_dir, "06_quantification/bambu/rc_output/",sample)

#check if bam file exists
if (!file.exists(bam_file)){
  print(bam_file)
  stop("Bam file does not exist")
}
#create output directory
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
print(sample)
read_class <- bambu(reads = paste0(bam_file),
                     annotations = annotation,
                     genome = fa.file,
                     rcOutDir = output_dir,
                     lowMemory = TRUE, 
                     discovery = FALSE, 
                     quant = FALSE)

sessioninfo::session_info()
