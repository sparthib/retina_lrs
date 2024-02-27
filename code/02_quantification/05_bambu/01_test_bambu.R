library(bambu)
library(BiocFileCache)

sample  <- commandArgs(trailingOnly = TRUE)
test.bam <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/MAPQ_FILTERED/", sample, "_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"

se_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/EP1-BRN3B-RO"


bambuAnnotations <- prepareAnnotations(gtf.file)
# se <- bambu(reads = test.bam, 
#             annotations = bambuAnnotations,
#             genome = fa.file,
#             ncore = 17,
#             rcOutDir = se_output_dir) 



bfc <- BiocFileCache(se_output_dir, ask = FALSE)
info <- bfcinfo(bfc)


# se <- bambu(reads = c("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/EP1-BRN3B-RO/32836422b0c117_32836422b0c117.rds"),
#             annotations = bambuAnnotations, genome = fa.file)
# Reads should either be: a vector of paths to .bam files, a vector of paths to Bambu RCfile .rds files, 
# or a list of loaded Bambu RCfiles
# pdf("/users/sparthib/retina_lrs/plots/bambu/ENSG00000116670.pdf")
# p <- plotBambu(se, type = "annotation", gene_id = "ENSG00000116670.16")
# print(p)
# dev.off()

# writeBambuOutput(se, 
#                  path = "/users/sparthib/retina_lrs/processed_data/bambu",
#                  prefix = "EP1-BRN3B-RO")

se.discoveryOnly <- bambu(reads = test.bam,
                          annotations = bambuAnnotations,
                          genome = fa.file,
                          quant = FALSE,
                          NDR = 1,
                          ncore = 17,
                          rcOutDir = se_output_dir)
                          
                          
writeBambuOutput(se.discoveryOnly, 
                 path = "/users/sparthib/retina_lrs/processed_data/bambu",
                 prefix = "EP1-BRN3B-RO_discovery_only")






