library("FLAMES")
library("sessioninfo")

gtf <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
genome_fa  <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa"
outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/ROs"
minimap_path <- file.path("/users/sparthib/minimap2")

#config_file <- FLAMES::create_config(outdir, do_barcode_demultiplex = FALSE)
config_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json"
config <- jsonlite::fromJSON("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/bulk.json")

se <- bulk_long_pipeline(
  annotation = gtf, outdir = outdir,
  genome_fa = genome_fa, config_file = config_file
)
saveRDS(se, file = paste0(outdir, "/se.rds"))


outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/FT_vs_RGC"
se <- bulk_long_pipeline(
  annotation = gtf, outdir = outdir,
  genome_fa = genome_fa, config_file = config_file
)
saveRDS(se, file = paste0(outdir, "/se.rds"))

sessioninfo::session_info()
