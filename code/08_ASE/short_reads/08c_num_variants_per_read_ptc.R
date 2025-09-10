library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)
library(readr)
library(rtracklayer)
# Read in BAM and VCF

vcf <- readVcf("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf",
               "hg38")
variants <- rowRanges(vcf)

samples <- c(
  "H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
  "H9-FT_1", "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2",
  "EP1-BRN3B-RO", "EP1-WT_hRO_2", "EP1-WT_ROs_D45"
)

array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_id <- 1

# Define BAM directory
genome_bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1"

gtf_dir <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
gtf <- import(gtf_dir)
# Filter for protein-coding genes
gtf <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

# Find overlaps
hits <- findOverlaps(variants, gtf)

# Subset variants that overlap with protein-coding genes
variants_ptc <- variants[queryHits(hits)]


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
h1_hits <- findOverlaps(h1_reads, variants_ptc)

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
h2_hits <- findOverlaps(h2_reads, variants_ptc)

# Count how many variants per read
h2_variant_counts <- table(queryHits(h2_hits))


# Convert to numeric
h2_counts_per_read <- as.integer(h2_variant_counts)


plot_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/vcf_stats/H9_EP1/variants_per_read_ptc"
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

write_tsv(data.frame(
  read_id = names(table(h1_counts_per_read)),
  variants_overlapped = table(h1_counts_per_read)
), file.path(plot_output_dir, paste0(sample, "_h1_variants_per_read_ptc.tsv")))

write_tsv(data.frame(
  read_id = names(table(h2_counts_per_read)),
  variants_overlapped = table(h2_counts_per_read)
), file.path(plot_output_dir, paste0(sample, "_h2_variants_per_read_ptc.tsv")))



