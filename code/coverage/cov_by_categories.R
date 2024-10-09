library(biomaRt)

#get transcript_coverage function from transcriptome_coverage.R
source("/users/sparthib/retina_lrs/code/coverage/transcriptome_coverage.R")

# Define BAM file directory and sample name
bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/transcriptome/GENCODE/supplementary_filtered"
sample <- commandArgs(trailingOnly = TRUE)[1]
bam <- file.path(bam_dir, paste0(sample, ".bam"))

# Read alignments from BAM file

output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/GAlignments/"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!file.exists(file.path(output_dir, paste0(sample, ".rds")))){
  aln <- GenomicAlignments::readGAlignments(bam, param = Rsamtools::ScanBamParam(mapqFilter = 5))
  # Save the alignments as an RDS file
  saveRDS(aln, file = file.path(output_dir, paste0(sample, ".rds")))
} else {
  aln <- readRDS(file.path(output_dir, paste0(sample, ".rds")))
}

#remove version num in seqnames
seqlevels(aln) <- gsub("\\..*", "", seqlevels(aln))
seqnames(aln) <- gsub("\\..*", "", seqnames(aln))

#get GC content of each transcript
us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_transcript_id",
                "hgnc_symbol", "percentage_gene_gc_content",
                "transcript_biotype", "transcript_length"),
  filter="ensembl_transcript_id",
  values=seqlevels(aln),
  uniqueRows=TRUE)
nrow(annotLookup)

#missing mart info on 221 transcripts 

# Get isoforms with more than 10 occurrences
isoform <- names(table(seqnames(aln)))[table(seqnames(aln)) > 10]
length(isoform)

length_bins = c(0, 1, 2, 5, 10, Inf)
#get coverage info of transcripts from aln
if (exists("transcript_coverage")) {
  print("The function 'transcript_coverage' exists!")
} else {
  print("The function 'transcript_coverage' does not exist.")
}

transcript_info <- transcript_coverage(aln, isoform, length_bins)

#merge transcript info with annotLookup
transcript_info$isoform <- rownames(transcript_info)

transcript_info <- merge(transcript_info, annotLookup, by.x = "isoform", 
                         by.y = "ensembl_transcript_id", all.x = TRUE)

# unique(transcript_info$isoform) |> head()
#check if
#mutate percentage_gene_gc_content to bins of length 10

transcript_info <- transcript_info |> 
  dplyr::mutate(gc_content_bin = cut(percentage_gene_gc_content, breaks = seq(0, 100, 10))) 

saveRDS(transcript_info, file = file.path(output_dir, paste0(sample_transcript_info, ".rds")))
# 
# # Generate coverage plots for each gc content bin
# output_plot_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/gc_content"
# if (!dir.exists(output_plot_dir)) {
#   dir.create(output_plot_dir, recursive = TRUE)
# }
# 
# for (bin in unique(transcript_info$gc_content_bin)) {
#   plot_coverage(output_plot_dir, transcript_info, sample, gc_content = bin)
# }
# 
# output_plot_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/transcript_biotype"
# if (!dir.exists(output_plot_dir)) {
#   dir.create(output_plot_dir, recursive = TRUE)
# }
# 
# for (bin in unique(transcript_info$transcript_biotype)) {
#   plot_coverage(output_plot_dir, transcript_info, sample, transcript_biotype = bin)
# }
# 
# output_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/length_bins"
# 
# # Ensure the directory exists
# if (!dir.exists(output_dir)) {
#   stop("Output directory does not exist: ", output_dir)
# }
# 
# # Generate coverage plots for each length bin
# length_bins_to_plot <- c("(10,Inf]", "(5,10]", "(2,5]", "(1,2]", "(0,1]")
# for (bin in length_bins_to_plot) {
#   plot_coverage(output_dir, transcript_info, sample, bin)
# }


sessioninfo::session_info()



