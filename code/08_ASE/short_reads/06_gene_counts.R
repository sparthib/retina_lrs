library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)

# Load BAM file
array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
            "H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2") 
sample <- samples[array_id]
print(paste("Processing sample:", sample))

whatshap_out_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output"
input_h1_bam <- file.path(whatshap_out_dir,
                          paste0(sample,"_h1.bam"))
input_h2_bam <- file.path(whatshap_out_dir,
                          paste0(sample,"_h2.bam"))

# Import GTF/GFF annotation
gtf_file <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly",
                      "release_46_primary_assembly.gtf")
gtf <- import(gtf_file)

print("Loaded GTF file")

all_gene_gtf <- gtf[gtf$type == "gene"]
genes <- gtf[gtf$type == "gene"]

output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gene_counts"

fc_1 <- featureCounts(
  files         = input_h1_bam,
  annot.ext     = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  nthreads      = 19,
  isLongRead     = TRUE,
  allowMultiOverlap = TRUE
)

print("Finished featureCounts for h1")

write.table(fc_1$counts, file = file.path(output_dir, paste0(sample, "_h1_counts.txt")),
            sep = "\t", quote = FALSE, col.names = NA)


fc_2 <- featureCounts(
  files         = input_h2_bam,
  annot.ext     = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType     = "exon",
  GTF.attrType        = "gene_id",
  nthreads      = 19,
  isLongRead     = TRUE,
  allowMultiOverlap = TRUE
)

print("Finished featureCounts for h2")

write.table(fc_2$counts, file = file.path(output_dir, paste0(sample, "_h2_counts.txt")),
            sep = "\t", quote = FALSE, col.names = NA)





