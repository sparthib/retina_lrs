library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)
library(readr)
library(rtracklayer)
library(biomaRt)
library(readxl)
library(dplyr)
library(stringr)
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

# Define BAM directory
genome_bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output_phased_on_H9_and_EP1"

gtf_dir <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
gtf <- import(gtf_dir)
# Filter for protein-coding genes
gtf <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

## get number of variants per ptc gene 
# ptc_hits <- findOverlaps(gtf, variants)
# variant_counts <- tabulate(queryHits(ptc_hits), nbins = length(gtf))
# mcols(gtf)$n_variants <- variant_counts

# get gene_id and n_variants dataframe
# gene_variant_counts <- data.frame(
#   gene_id = gtf$gene_id,
#   gene_name = gtf$gene_name,
#   start = start(gtf),
#   end = end(gtf),
#   n_variants = mcols(gtf)$n_variants
# )

# get gened 

# readr::write_tsv(gene_variant_counts, file = "/users/sparthib/retina_lrs/processed_data/ASE/vcf_stats/number_of_variants_per_ptc_gene.tsv")


#### IRD Genes ######

#from RetNet database 
# disease_gene_vector <- c("ADIPOR1", "ARL6", "BBIP1", "BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBS10", "BBS12", "C8orf37", "CEP19", "CEP290", "IFT172", "IFT27", "INPP5E", "LZTFL1", "MKKS", "MKS1", "NPHP1", "SDCCAG8", "TRIM32", "TTC8",
#                  "PRDM13", "RGR", "TEAD1",
#                  "AIPL1", "CRX", "GUCA1A", "GUCY2D", "PITPNM3", "PROM1", "PRPH2", "RIMS1", "SEMA4A", "UNC119",
#                  "ABCA4", "ADAM9", "ATF6", "C21orf2", "C8orf37", "CACNA2D4", "CDHR1", "CEP78", "CERKL", "CNGA3", "CNGB3", "CNNM4", "DYNC2I2", "GNAT2", "IFT81", "KCNV2", "PDE6C", "PDE6H", "POC1B", "RAB28", "RAX2", "RDH5", "RPGRIP1", "SLC4A7", "TTLL5",
#                  "CACNA1F", "RPGR",
#                  "GNAT1", "PDE6B", "RHO",
#                  "CABP4", "GNAT1", "GNB3", "GPR179", "GRK1", "GRM6", "LRIT3", "RDH5", "SAG", "SLC24A1", "TRPM1",
#                  "CACNA1F", "NYX",
#                  "ESPN", "WFS1",
#                  "CDH23", "CIB2", "ESPN", "MYO7A", "PCDH15", "PDZD7", "USH1C", "WHRN",
#                  "CRX", "IMPDH1", "OTX2",
#                  "AIPL1", "CABP4", "CCT2", "CEP290", "CLUAP1", "CRB1", "CRX", "DTHD1", "GDF6", "GUCY2D", "IFT140", "IQCB1", "KCNJ13", "LCA5", "LRAT", "NMNAT1", "PRPH2", "RD3", "RDH12", "RPE65", "RPGRIP1", "SPATA7", "TULP1",
#                  "BEST1", "C1QTNF5", "CTNNA1", "EFEMP1", "ELOVL4", "FSCN2", "GUCA1B", "HMCN1", "IMPG1", "LRRTM4", "OTX2", "PRDM13", "PROM1", "PRPH2", "RP1L1", "TIMP3",
#                  "ABCA4", "CFH", "DRAM2", "IMPG1", "MFSD8", "SLC37A3",
#                  "RPGR", "VCAN",
#                  "AFG3L2", "MFN2", "MIEF1", "NR2F1", "OPA1",
#                  "ACO2", "NBAS", "RTN4IP1", "TMEM126A", "TIMM8A",
#                  "ADIPOR1", "ARL3", "BEST1", "CA4", "CRX", "FSCN2", "GUCA1B", "HK1", "IMPDH1", "IMPG1", "KIF3B", "KLHL7", "NR2E3", "NRL", "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF31", "PRPH2", "RDH12", "RHO", "ROM1", "RP1", "RP9", "RPE65", "SAG", "SEMA4A", "SNRNP200", "SPP2", "TOPORS",
#                  "ABCA4", "ADGRA3", "AGBL5", "AHR", "ARHGEF18", "ARL6", "ARL2BP", "BBS1", "BBS2", "BEST1", "C8orf37", "CC2D2A", "CERKL", "CLCC1", "CLRN1", "CNGA1", "CNGB1", "COQ2", "COQ4", "COQ5", "CRB1", "CWC27", "CYP4V2", "DHDDS", "DHX38", "EMC1", "ENSA", "EYS", "FAM161A", "HGSNAT", "IDH3B", "IFT140", "IFT172", "IMPG2", "KIAA1549", "KIZ", "LRAT", "MAK", "MERTK", "MVK", "NEK2", "NEUROD1", "NR2E3", "NRL", "PCARE", "PDE6A", "PDE6B", "PDE6G", "PDSS1", "POMGNT1", "PRCD", "PROM1", "PROS1", "RAX2", "RBP3", "REEP6", "RGR", "RHO", "RLBP1", "RP1", "RP1L1", "RPE65", "SAG", "SAMD11", "SLC37A3", "SLC39A12", "SLC66A1", "SLC7A14", "SPATA7", "TRNT1", "TTC8", "TULP1", "USH2A", "ZNF408", "ZNF513",
#                  "OFD1", "RP2", "RPGR",
#                  "ABCC6", "AFG3L2", "ATXN7", "COL11A1", "COL2A1", "JAG1", "KCNJ13", "KIF11", "MFN2", "OPA3", "PAX2", "TREX1", "VCAN",
#                  "ABCC6", "ABHD12", "ACBD5", "ACO2", "ADAMTS18", "ADIPOR1", "AFG3L2", "AHI1", "ALMS1", "CC2D2A", "CEP164", "CEP290", "CLN3", "COL9A1", "CSPP1", "CWC27", "ELOVL4", "EXOSC2", "FLVCR1", "GNPTG", "HARS", "HGSNAT", "HMX1", "IFT140", "IFT81", "INPP5E", "INVS", "IQCB1", "LAMA1", "LRP5", "MKS1", "MTTP", "NPHP1", "NPHP3", "NPHP4", "OPA3", "PANK2", "PCYT1A", "PDSS1", "PEX1", "PEX2", "PEX7", "PHYH", "PLK4", "PNPLA6", "POC5", "POC1B", "PPT1", "PRPS1", "RDH11", "RIMS2", "RPGRIP1L", "SDCCAG8", "SLC25A46", "TMEM216", "TMEM237", "TRNT1", "TTPA", "TUB", "TUBGCP4", "TUBGCP6", "WDPCP", "WDR19", "WFS1", "ZNF423",
#                  "OFD1", "TIMM8A",
#                  "ABHD12", "ADGRV1", "ARSG", "CDH23", "CEP250", "CEP78", "CIB2", "CLRN1", "ESPN", "HARS", "MYO7A", "PCDH15", "USH1C", "USH1G", "USH2A", "WHRN",
#                  "BEST1", "CAPN5", "CRB1", "ELOVL1", "FZD4", "ITM2B", "KIF3B", "LRP5", "MAPKAPK3", "MIR204", "OPN1SW", "RB1", "RCBTB1", "RGR", "TSPAN12", "ZNF408",
#                  "ASRGL1", "BEST1", "C12orf65", "CDH3", "CNGA3", "CNGB3", "CNNM4", "COQ2", "CYP4V2", "DYNC2H1", "LRP5", "MFRP", "MVK", "NBAS", "NR2E3", "OAT", "PLA2G5", "PROM1", "RBP4", "RCBTB1", "RGS9", "RGS9BP", "RLBP1",
#                  "KSS", "LHON", "MT-ATP6", "MT-TH", "MT-TL1", "MT-TP", "MT-TS2",
#                  "CACNA1F", "CHM", "DMD", "NDP", "OPN1LW", "OPN1MW", "PGK1", "RS1"
# )
# 
# disease_gene_vector <- disease_gene_vector |> unique()
# length(disease_gene_vector)
# 
# # Join with gene_variant_counts to get number of variants per disease gene
# disease_gene_variant_counts <- gene_variant_counts |>
#   filter(gene_name %in% disease_gene_vector)
# #get if there are duplicate gene names in gene_variant_counts
# sum(duplicated(disease_gene_variant_counts$gene_name))
# #which are the duplicated gene names
# disease_gene_variant_counts[duplicated(disease_gene_variant_counts$gene_name), ]
# 
# 
# nrow(disease_gene_variant_counts)
# # 303
# 
# readr::write_tsv(disease_gene_variant_counts, file = "/users/sparthib/retina_lrs/processed_data/ASE/vcf_stats/number_of_variants_per_ird_gene.tsv")


# Create full BAM file paths
h1_bam_files <- file.path(genome_bam_dir, paste0(samples, "_h1.bam"))
h2_bam_files <- file.path(genome_bam_dir,
                          paste0(samples, "_h2.bam"))


# Check result
print(h1_bam_files)
print(h2_bam_files)




# Open BAM
h1_bam <- BamFile(h1_bam_files[array_id])
bam1_seqlevels <- seqlevels(h1_bam)

# 1. Filter gtf rows to keep only seqnames in BAM
ptc_gtf1 <- gtf[seqnames(gtf) %in% bam1_seqlevels]

# 2. Drop unused seqlevels so param only has chromosomes in BAM
ptc_gtf1 <- keepSeqlevels(ptc_gtf1, intersect(seqlevels(ptc_gtf1), bam1_seqlevels),
                          pruning.mode = "coarse")

# 3. Now create ScanBamParam
param <- ScanBamParam(which = ptc_gtf1)

# 4. Read alignments
indexBam(h1_bam_files[array_id])

h1_reads <- readGAlignments(h2_bam, param = param, index = paste0(h1_bam_files[array_id], ".bai"))


#total number of aligned reads
print(paste("Total aligned reads in h1:", length(h1_reads)))

# Find overlaps
var1_hits <- findOverlaps(variants, ptc_gtf1)
variants_ptc1 <- variants[queryHits(var1_hits)]
h1_hits <- findOverlaps(h1_reads,variants_ptc1)

# Count how many variants per read
h1_variant_counts <- table(queryHits(h1_hits))


# Convert to numeric
h1_counts_per_read <- as.integer(h1_variant_counts)


# Extract reads and their positions
h2_bam <- BamFile(h2_bam_files[array_id])
bam2_seqlevels <- seqlevels(h2_bam)

# 2. Filter gtf rows to keep only seqnames in BAM
ptc_gtf2 <- gtf[seqnames(gtf) %in% bam2_seqlevels]

# 2. Drop unused seqlevels so param only has chromosomes in BAM
ptc_gtf2 <- keepSeqlevels(ptc_gtf2, intersect(seqlevels(ptc_gtf2), bam2_seqlevels),
                          pruning.mode = "coarse")

# 3. Now create ScanBamParam
param <- ScanBamParam(which = ptc_gtf2)

# 4. Read alignments
indexBam(h2_bam_files[array_id])

h2_reads <- readGAlignments(h2_bam, param = param, index = paste0(h2_bam_files[array_id], ".bai"))

#total number of aligned reads
print(paste("Total aligned reads in h2:", length(h2_reads)))

# Find overlaps
var2_hits <- findOverlaps(variants, ptc_gtf2)
variants_ptc2 <- variants[queryHits(var2_hits)]
h2_hits <- findOverlaps(h2_reads,variants_ptc2)


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



