library(bambu)
library(BiocFileCache)
library(readr)
library(sessioninfo)
library(dplyr)

# gtf.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
# bambuAnnotations <- prepareAnnotations(gtf.file)
# saveRDS(bambuAnnotations, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

annotation <- readRDS("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/annotations.rds")

sample_config <- read_tsv(file ="/users/sparthib/retina_lrs/config.tsv")

sample_config <- sample_config |> dplyr::select(-c(1))
colnames(sample_config) <- c("sample_name", "fastq_path", "summary_stats_path")
## create list of bam file paths for all samples 
samples <- sample_config$sample_name
samples <- c(samples, "H9-FT_1","H9-FT_2","H9-hRGC_1","H9-hRGC_2" )
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/"


input_names <- c()
 for(i in samples){ 
  input_names <- c(input_names, paste0(bam_dir, i, "_sorted.bam"))
  }

fa.file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
se_output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/"


se_discoveryOnly_multisample <- bambu(reads = input_names,
                          annotations = annotation,
                          genome = fa.file,
                          quant = FALSE,
                          NDR = 1)


writeBambuOutput(se_discoveryOnly_multisample, 
                 path = se_output_dir,
                 prefix = "multi_sample_mar_29")

sessioninfo::session_info()
                          

# annotations.filtered <- se.discoveryOnly[(!is.na(mcols(se.discoveryOnly)$NDR) & mcols(se.discoveryOnly)$NDR <
#                                           0.1) | is.na(mcols(se.discoveryOnly)$NDR)]
# 
# 
# se.NDR_1 <- bambu(reads = bam, annotations = annotations.filtered, genome = fa.file,
#                   NDR = 1, discovery = FALSE)





# se.discoveryOnly <- bambu(reads = c("test_rc.rds"),
#                           annotations = bambuAnnotations, genome = fa.file)
# 
# annotations.filtered <- se.discoveryOnly[(!is.na(mcols(newAnnotations)$NDR) & mcols(newAnnotations)$NDR <
#                                           0.1) | is.na(mcols(newAnnotations)$NDR)]
# 
# 

# $warnings
# $warnings[[1]]
# [1] "not all chromosomes present in reference annotations, annotations might be incomplete. Please compare objects on the same reference"
# [2] "25701 reads are mapped outside the provided genomic regions. These reads will be dropped. Check you are using the same genome used for the alignment"
# [3] "No aligned spliced reads detected!Bambu expects spliced reads. If this is intended, see Documentation on how to handle single-exon transcripts"


