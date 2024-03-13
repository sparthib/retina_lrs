library(bambu)
library(BiocFileCache)

sample  <- commandArgs(trailingOnly = TRUE)
bam <- paste0("/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/", "H9-hRGC_1", "_sorted.bam")
fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"

se_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/H9-hRGC_1"


bambuAnnotations <- prepareAnnotations(gtf.file)

bfc <- BiocFileCache(se_output_dir, ask = FALSE)
info <- bfcinfo(bfc)

se.discoveryOnly <- bambu(reads = bam,
                          annotations = bambuAnnotations,
                          genome = fa.file,
                          quant = FALSE,
                          NDR = 1)
                          

writeBambuOutput(se.discoveryOnly, 
                 path = "/users/sparthib/retina_lrs/processed_data/bambu",
                 prefix = "test_mar_13")

annotations.filtered <- se.discoveryOnly[(!is.na(mcols(newAnnotations)$NDR) & mcols(newAnnotations)$NDR <
                                          0.1) | is.na(mcols(newAnnotations)$NDR)]



# $warnings
# $warnings[[1]]
# [1] "not all chromosomes present in reference annotations, annotations might be incomplete. Please compare objects on the same reference"
# [2] "25701 reads are mapped outside the provided genomic regions. These reads will be dropped. Check you are using the same genome used for the alignment"
# [3] "No aligned spliced reads detected!Bambu expects spliced reads. If this is intended, see Documentation on how to handle single-exon transcripts"


