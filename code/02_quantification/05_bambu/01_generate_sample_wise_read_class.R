library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/annotations.rds")

sample <- commandArgs(trailingOnly = TRUE)

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/annotations.rds")
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality"
bam_dir <- paste0(bam_dir, sample, "_primary_over_30_chr_only_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"

output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/primary_assembly/", sample, "/")

#check if bam file exists
if (!file.exists(bam_dir)){
  stop("Bam file does not exist")
}

#create output directory
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
print(sample)

read_class <- bambu(reads = paste0(bam_dir),
                     annotations = annotation,
                     genome = fa.file,
                     lowMemory = TRUE, 
                     discovery = FALSE, 
                     quant = FALSE)

saveRDS(read_class, paste0(output_dir, "read_class.rds"))

sessioninfo::session_info()
