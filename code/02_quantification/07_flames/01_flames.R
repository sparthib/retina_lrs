library("FLAMES")
library("sessioninfo")

sample <- commandArgs(trailingOnly = TRUE)
gtf <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
genome_fa  <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"
outdir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/", sample)
genome_bam <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/H9-hRGC_1_sorted.bam"
minimap_path <- file.path("/users/sparthib/minimap2")
fastq <- paste0(output_dir, "/matched_reads.fastq")

#config_file <- FLAMES::create_config(outdir, do_barcode_demultiplex = FALSE)
config_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json"
config <- jsonlite::fromJSON("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json")

se <- bulk_long_pipeline(
  annot = annot, fastq = fastq, outdir = outdir,
  genome_fa = genome_fa, config_file = config_file
)
saveRDS(se, file = paste0(outdir, "/se.rds"))


sessioninfo::session_info()
