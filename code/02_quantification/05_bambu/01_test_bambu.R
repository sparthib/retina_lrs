library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/"
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
se_output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/", sample, "/")

samples <- c("EP1-BRN3B-RO", "H9-BRN3B-RO", "hRGC", "DG-WT-hRGC",
                          "H9-CRX_ROs_D45", "EP1-WT_ROs_D45", "YZ-3T_hRGC",
                          "YZ-15T_hRGC", "H9-FT_1","H9-FT_2", "H9-hRGC_2")

for (sample in samples){ 
  print(sample)
  novelOnly <- bambu(reads = paste0(bam_dir, sample, "_sorted.bam"),
                     annotations = annotation,
                     genome = fa.file,
                     quant = FALSE, NDR = 1)
  
  writeBambuOutput(novelOnly, 
                   path = se_output_dir,
                   prefix = sample)
  }

sessioninfo::session_info()
                          
