library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)
library(readr)

# get number of variants per read in each of the haplotagged bam files

# Read in BAM and VCF
code_dir <- Sys.getenv("retina_lrs_code")
data_dir <- Sys.getenv("retina_lrs_dir")
vcf <- readVcf(file.path(data_dir, "09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf"),
               "hg38")
variants <- rowRanges(vcf)

samples <- c(
  "H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
  "H9-FT_1", "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2",
  "EP1-BRN3B-RO", "EP1-WT_hRO_2", "EP1-WT_ROs_D45"
)

array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Define BAM directory
genome_bam_dir <- file.path(data_dir,"09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1")

# Create full BAM file paths
h1_bam_files <- file.path(genome_bam_dir, paste0(samples, "_h1.bam"))
h2_bam_files <- file.path(genome_bam_dir,
                          paste0(samples, "_h2.bam"))


# Check result
print(h1_bam_files)
print(h2_bam_files)



# Extract reads and their positions
h1_bam <- BamFile(h1_bam_files[array_id])
h1_reads <- readGAlignments(h1_bam)
#total number of aligned reads
print(paste("Total aligned reads in h1:", length(h1_reads)))

# Find overlaps
h1_hits <- findOverlaps(h1_reads, variants)

# Count how many variants per read
h1_variant_counts <- table(queryHits(h1_hits))


# Convert to numeric
h1_counts_per_read <- as.integer(h1_variant_counts)


# Extract reads and their positions
h2_bam <- BamFile(h2_bam_files[array_id])
h2_reads <- readGAlignments(h2_bam)
#total number of aligned reads
print(paste("Total aligned reads in h2:", length(h2_reads)))

# Find overlaps
h2_hits <- findOverlaps(h2_reads, variants)

# Count how many variants per read
h2_variant_counts <- table(queryHits(h2_hits))


# Convert to numeric
h2_counts_per_read <- as.integer(h2_variant_counts)


plot_output_dir <- file.path(code_dir,"processed_data/ASE/vcf_stats/H9_EP1/variants_per_read")

write_tsv(data.frame(
  read_id = names(table(h1_counts_per_read)),
  variants_overlapped = table(h1_counts_per_read)
), file.path(plot_output_dir, paste0(sample, "_h1_variants_per_read.tsv")))

write_tsv(data.frame(
  read_id = names(table(h2_counts_per_read)),
  variants_overlapped = table(h2_counts_per_read)
), file.path(plot_output_dir, paste0(sample, "_h2_variants_per_read.tsv")))






