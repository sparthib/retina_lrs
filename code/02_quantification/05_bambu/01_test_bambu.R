library(bambu)


sample  <- commandArgs(trailingOnly = TRUE)
test.bam <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/MAPQ_FILTERED/", sample, "_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"

se_output_dir <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/", sample)

bambuAnnotations <- prepareAnnotations(gtf.file)
se <- bambu(reads = test.bam, 
            annotations = bambuAnnotations,
            genome = fa.file,
            ncore = 17,
            rcOutDir = se_output_dir) 


# quant = FALSE ## transcript discovery only
# discovery = FALSE ## quantification only 

