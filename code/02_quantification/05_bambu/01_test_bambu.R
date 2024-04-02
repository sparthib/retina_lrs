library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

# sample_config <- read_tsv(file ="/users/sparthib/retina_lrs/config.tsv")

samples <- c("EP1-BRN3B-RO", "H9-BRN3B-RO", "hRGC", "DG-WT-hRGC",
             "H9-CRX_ROs_D45", "EP1-WT_ROs_D45", "YZ-3T_hRGC",
             "YZ-15T_hRGC", "H9-FT_1","H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/"



fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
se_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/"


se_quantOnly <- bambu(reads = paste0(bam_dir, samples[task_id], "_sorted.bam"),
                                  annotations = annotation,
                                  genome = fa.file,
                                  quant = TRUE)

writeBambuOutput(se_quantOnly, 
                 path = se_output_dir,
                 prefix = samples[task_id])

sessioninfo::session_info()
                          
