#!/usr/bin/env Rscript
# 08_isoform_SR_LR_comparison.R
# Compare KNOWN (ENST) isoforms detected in short-read (salmon) vs long-read (bambu),
# and identify isoforms detected only in long reads.
#
# Two comparisons:
#   (1) All known isoforms      : SR-detected ENST (all biotypes) vs LR-detected ENST (all biotypes)
#   (2) Protein-coding isoforms : both sides restricted to protein_coding ENST
# In each, the LR-only set = ENST detected in LR but not in SR.
#
# Only known GENCODE isoforms (ENST...) are comparable: the short-read salmon quant
# used a plain GENCODE release-46 index (ENST only), while long-read bambu used an
# EXTENDED annotation that also contains novel isoforms (Bambu... IDs) which short
# reads cannot see by construction. Novel LR isoforms are reported separately.
#
# LR reference = the broader common-isoform-filtered bambu matrix (detection set, all
# biotypes present), NOT the protein-coding+expression-filtered DE matrix, so both the
# all-biotype and protein-coding comparisons come from one detection-based matrix.
#
# Transcript biotype comes from references .../primary_assembly/transcript_meta_info.rds
# (columns: isoform_id, GC_percent, transcript_length, transcript_biotype).
# Detection rule (both platforms): >=1 estimated count in >=1 sample.
# Run with: module load conda_R/4.5.x && Rscript 08_isoform_SR_LR_comparison.R

.libPaths("/users/sparthib/R/4.5.x")

data_dir <- "/dcs04/hicks/data/sparthib/retina_lrs"
ref_dir  <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
sr_file  <- file.path(data_dir, "10_short_reads/salmon/transcript_counts_matrix.csv")
lr_file  <- file.path(data_dir, "06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS")
meta_file<- file.path(ref_dir,  "transcript_meta_info.rds")
out_dir  <- file.path(data_dir, "10_short_reads/salmon")

strip <- function(x) sub("\\.[0-9]+$", "", x)

## ENST -> biotype map
meta <- readRDS(meta_file)
meta$isoform_id <- strip(as.character(meta$isoform_id))
biotype <- setNames(as.character(meta$transcript_biotype), meta$isoform_id)
cat("meta transcripts:", nrow(meta),
    " protein_coding:", sum(meta$transcript_biotype == "protein_coding", na.rm = TRUE), "\n")

## SR detected known isoforms (>=1 count in >=1 of 28 samples)
sr <- as.matrix(read.csv(sr_file, row.names = 1, check.names = FALSE))
rownames(sr) <- strip(rownames(sr))
sr <- sr[grepl("^ENST", rownames(sr)), , drop = FALSE]
sr_det <- rownames(sr)[rowSums(sr >= 1) >= 1]
cat("SR: known(ENST) annotated", nrow(sr), " detected", length(sr_det), "\n")

## LR detected known isoforms (>=1 count in >=1 of 7 samples)
lr <- as.matrix(readRDS(lr_file))
rownames(lr) <- strip(rownames(lr))
lr_enst   <- rownames(lr)[grepl("^ENST", rownames(lr))]
lr_novel  <- rownames(lr)[grepl("^Bambu", rownames(lr), ignore.case = TRUE)]
lr_det    <- rownames(lr)[rowSums(lr >= 1) >= 1]
lr_enst_det  <- intersect(lr_det, lr_enst)
lr_novel_det <- intersect(lr_det, lr_novel)
cat("LR: known(ENST) rows", length(lr_enst), " detected", length(lr_enst_det),
    " | novel(Bambu) detected", length(lr_novel_det), "\n")

compare <- function(sr_set, lr_set, tag) {
  both    <- intersect(sr_set, lr_set)
  lr_only <- setdiff(lr_set, sr_set)
  sr_only <- setdiff(sr_set, lr_set)
  cat("\n== ", tag, " ==\n", sep = "")
  cat("SR detected:", length(sr_set), " LR detected:", length(lr_set), "\n")
  cat("in BOTH:", length(both), " LR-only:", length(lr_only), " SR-only:", length(sr_only), "\n")
  df <- data.frame(transcript_id = lr_only, biotype = biotype[lr_only], row.names = NULL)
  write.csv(df, file.path(out_dir, paste0("isoform_LRonly_", tag, ".csv")), row.names = FALSE)
  c(SR = length(sr_set), LR = length(lr_set), both = length(both),
    LR_only = length(lr_only), SR_only = length(sr_only))
}

r1 <- compare(sr_det, lr_enst_det, "all_known")
pc <- names(biotype)[biotype == "protein_coding" & !is.na(biotype)]
r2 <- compare(intersect(sr_det, pc), intersect(lr_enst_det, pc), "protein_coding")

summ <- rbind(data.frame(comparison = "all_known", t(r1)),
              data.frame(comparison = "protein_coding", t(r2)))
summ$LR_novel_detected <- length(lr_novel_det)
write.csv(summ, file.path(out_dir, "isoform_SR_LR_overlap_summary.csv"), row.names = FALSE)
cat("\nSUMMARY\n"); print(summ); cat("\nDONE\n")
