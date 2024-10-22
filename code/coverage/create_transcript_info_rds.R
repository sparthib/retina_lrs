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


# Save the transcript info as an RDS file
output_rds_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/transcript_info"
if (!dir.exists(output_rds_dir)) {
  dir.create(output_rds_dir, recursive = TRUE)
}


# Generate coverage plots for each gc content bin
output_plot_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/transcript_biotype"
if (!dir.exists(output_plot_dir)) {
  dir.create(output_plot_dir, recursive = TRUE)
}



