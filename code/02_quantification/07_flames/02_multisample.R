library("FLAMES")
library("sessioninfo")

gtf <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
genome_fa  <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"
outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/ROs"
minimap_path <- file.path("/users/sparthib/minimap2")

#config_file <- FLAMES::create_config(outdir, do_barcode_demultiplex = FALSE)
config_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json"
config <- jsonlite::fromJSON("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json")

RO_fastq_path = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/ROs/fastqs"
RO_fastqs <- list.files(RO_fastq_path, full.names = TRUE)

#name filepaths in list 
names(RO_fastqs) <- c("EP1-BRN3B-RO", "EP1-WT_hRO_2", "EP1-WT_ROs_D45",
                      "H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45")



se <- bulk_long_pipeline(
  annotation = gtf, 
  fastq= RO_fastq_path,
  outdir = outdir,
  genome_fa = genome_fa, config_file = config_file
)
saveRDS(se, file = paste0(outdir, "/se.rds"))

FT_vs_RGC_fastq_path = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC/fastqs/"
FT_vs_RGC_fastqs <- list.files(FT_vs_RGC_fastq_path, full.names = TRUE)
names(FT_vs_RGC_fastqs) <- c("H9-FT_1","H9-FT_2" , "H9-hRGC_1" ,"H9-hRGC_2")
 
outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC"
se <- bulk_long_pipeline(
  annotation = gtf,
  fastq = FT_vs_RGC_fastq_path,
  outdir = outdir,
  genome_fa = genome_fa, config_file = config_file
)
saveRDS(se, file = paste0(outdir, "/se.rds"))

sessioninfo::session_info()
