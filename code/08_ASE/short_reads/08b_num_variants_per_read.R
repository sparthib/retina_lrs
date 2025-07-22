
library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)

# Read in BAM and VCF

vcf <- readVcf("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf",
               "hg38")
variants <- rowRanges(vcf)

# Define samples
array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

samples <- c(
  "H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
  "H9-FT_1", "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2",
  "EP1-BRN3B-RO", "EP1-WT_hRO_2", "EP1-WT_ROs_D45"
)

sample <- samples[array_id]
# Define BAM directory
genome_bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality"

# Create full BAM file paths
bam_files <- file.path(genome_bam_dir, paste0(samples, "_primary_over_30_chr_only_sorted.bam"))

# Check result
print(bam_files)

# Extract reads and their positions
bam <- BamFile(bam_files[array_id])
reads <- readGAlignments(bam)

# Find overlaps
hits <- findOverlaps(reads, variants)

# Count how many variants per read
variant_counts <- table(queryHits(hits))


# Convert to numeric
counts_per_read <- as.integer(variant_counts)

plot_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/vcf_stats/H9_EP1/variants_per_read"
dir.create(plot_output_dir, showWarnings = FALSE)
# Plot histogram
pdf(file.path(plot_output_dir, paste0(sample,
                                      "variants_per_read_histogram.pdf")))
hist(counts_per_read,
     breaks = 30,
     main = "Number of Variants per Read",
     xlab = "Variants Overlapped",
     ylab = "Number of Reads",
     col = "steelblue",
     border = "white")

dev.off()





