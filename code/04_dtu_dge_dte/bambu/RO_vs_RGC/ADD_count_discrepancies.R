

library(rtracklayer)
### bambu annotation 
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
gtf <- file.path(bambu_dir, "extended_annotations.gtf")

### load gtf
gtf <- import(gtf)

### keep only type "transcript"
gtf <- gtf[gtf$type == "transcript"]

## get rows where gene_id has the pattern "ENSG00000087274"
transcript_ids <- gtf[gtf$gene_id == "ENSG00000087274.19"]$transcript_id



###### get the counts for all these transcripts from bambu output 

# [1] "ENST00000683351.1"  "ENST00000264758.11" "ENST00000398125.5" 
# [4] "ENST00000508684.5"  "ENST00000510101.5"  "BambuTx210"        
# [7] "ENST00000511797.5"  "ENST00000539108.5"  "ENST00000513328.6" 
# [10] "ENST00000509039.5"  "ENST00000508277.5"  "ENST00000503455.6" 
# [13] "ENST00000540541.1"  "ENST00000355842.7"  "ENST00000651918.1" 
# [16] "ENST00000374281.2"  "ENST00000503169.5"  "ENST00000398123.6" 
# [19] "ENST00000398129.5"  "ENST00000534870.5"  "ENST00000506157.1" 
# [22] "ENST00000514940.5"  "ENST00000536078.1"  "ENST00000541051.5" 
# [25] "ENST00000513762.2"  "ENST00000536424.2"  "ENST00000539149.1" 
# [28] "ENST00000541843.2"  "ENST00000503062.1"  "ENST00000538860.1"

add_isoforms <- c("ENST00000683351.1", "ENST00000264758.11", "ENST00000398125.5", 
                  "ENST00000508684.5", "ENST00000510101.5", "BambuTx210", 
                  "ENST00000511797.5", "ENST00000539108.5", "ENST00000513328.6", 
                  "ENST00000509039.5", "ENST00000508277.5", "ENST00000503455.6", 
                  "ENST00000540541.1", "ENST00000355842.7", "ENST00000651918.1", 
                  "ENST00000374281.2", "ENST00000503169.5", "ENST00000398123.6", 
                  "ENST00000398129.5", "ENST00000534870.5", "ENST00000506157.1", 
                  "ENST00000514940.5", "ENST00000536078.1", "ENST00000541051.5", 
                  "ENST00000513762.2", "ENST00000536424.2", "ENST00000539149.1", 
                  "ENST00000541843.2", "ENST00000503062.1" ,  "ENST00000538860.1")


isoform_counts <- read.table(file.path(bambu_dir, "counts_transcript.txt"), 
                   header = TRUE, row.names = 1, sep = "\t",
                   comment.char = "")


isoform_counts <- isoform_counts[add_isoforms, ]

#remove first column
isoform_counts <- isoform_counts[, -1]

isoform_counts <- colSums(isoform_counts) 

gene_counts <- read.table(file.path(bambu_dir, "counts_gene.txt"), 
                   header = TRUE, row.names = 1, sep = "\t",
                   comment.char = "")

gene_counts <- gene_counts[rownames(gene_counts) == "ENSG00000087274.19", ]

# https://github.com/GoekeLab/bambu/issues/440
# These are read counts that could not be associated to any transcript belonging to the gene, 
# but would still be included in the gene count.




