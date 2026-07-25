#!/usr/bin/env Rscript
# 08_isoform_SR_LR_comparison.R
# Compare KNOWN (ENST) isoforms detected in short-read (salmon) vs long-read (bambu),
# identify isoforms detected ONLY in long reads, and annotate them with their genes.
#
# Two comparisons:
#   (1) All known isoforms      : SR-detected ENST (all biotypes) vs LR-detected ENST (all biotypes)
#   (2) Protein-coding isoforms : both sides restricted to protein_coding ENST
# In each, LR-only = ENST detected in LR but not in SR. LR-only lists are annotated
# with gene_id (ENSG) and gene_name, and summarized per gene.
#
# Only known GENCODE isoforms (ENST...) are comparable: short-read salmon quant used a
# plain GENCODE release-46 index (ENST only); long-read bambu used an EXTENDED annotation
# that also contains novel isoforms (Bambu...) short reads cannot see. Novel LR isoforms
# are reported separately.
#
# LR reference = broader common-isoform-filtered bambu matrix (detection set, all biotypes).
# Biotype from transcript_meta_info.rds; gene mapping from bambu counts_transcript.txt
# (TXNAME->GENEID) and gene symbols from release-46 transcript FASTA headers.
# Detection rule (both platforms): >=1 estimated count in >=1 sample.
# Run with: module load conda_R/4.5.x && Rscript 08_isoform_SR_LR_comparison.R

.libPaths("/users/sparthib/R/4.5.x")

data_dir <- "/dcs04/hicks/data/sparthib/retina_lrs"
ref_dir  <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
sr_file  <- file.path(data_dir, "10_short_reads/salmon/transcript_counts_matrix.csv")
lr_file  <- file.path(data_dir, "06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS")
meta_file<- file.path(ref_dir,  "transcript_meta_info.rds")
bambu_tx <- file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/counts_transcript.txt")
fa_file  <- file.path(ref_dir,  "release_46_all_transcripts.fa")
out_dir  <- file.path(data_dir, "10_short_reads/salmon")

strip <- function(x) sub("\\.[0-9]+$", "", x)

## ENST -> biotype
meta <- readRDS(meta_file)
meta$isoform_id <- strip(as.character(meta$isoform_id))
biotype <- setNames(as.character(meta$transcript_biotype), meta$isoform_id)

## ENST -> ENSG (from bambu transcript table: TXNAME, GENEID); covers all LR transcripts
b <- read.table(bambu_tx, header = TRUE, sep = "\t", comment.char = "",
                colClasses = c("character","character", rep("NULL", 11)))
tx2g <- setNames(strip(b$GENEID), strip(b$TXNAME))

## ENSG -> gene_name  and  ENST -> gene_name  (release-46 FASTA header field 6)
## header: >ENST|ENSG|OTTHUMG|OTTHUMT|txname|GENE_NAME|length|biotype|
hdr <- system(paste0("grep '^>' ", fa_file, " | sed 's/^>//'"), intern = TRUE)
p   <- strsplit(hdr, "|", fixed = TRUE)
h_enst <- strip(vapply(p, `[`, "", 1))
h_ensg <- strip(vapply(p, function(z) if (length(z) >= 2) z[2] else NA_character_, ""))
h_name <- vapply(p, function(z) if (length(z) >= 6) z[6] else NA_character_, "")
ensg2name <- setNames(h_name, h_ensg)          # ENSG -> symbol
cat("maps: biotype", length(biotype), " tx2gene", length(tx2g), " ensg2name", length(unique(h_ensg)), "\n")

## SR detected known isoforms
sr <- as.matrix(read.csv(sr_file, row.names = 1, check.names = FALSE))
rownames(sr) <- strip(rownames(sr))
sr <- sr[grepl("^ENST", rownames(sr)), , drop = FALSE]
sr_det <- rownames(sr)[rowSums(sr >= 1) >= 1]

## LR detected known isoforms
lr <- as.matrix(readRDS(lr_file)); rownames(lr) <- strip(rownames(lr))
lr_enst  <- rownames(lr)[grepl("^ENST", rownames(lr))]
lr_novel <- rownames(lr)[grepl("^Bambu", rownames(lr), ignore.case = TRUE)]
lr_det   <- rownames(lr)[rowSums(lr >= 1) >= 1]
lr_enst_det  <- intersect(lr_det, lr_enst)
lr_novel_det <- intersect(lr_det, lr_novel)
cat("SR known detected", length(sr_det), " LR known detected", length(lr_enst_det),
    " LR novel detected", length(lr_novel_det), "\n")

annotate <- function(ids) {
  g <- tx2g[ids]
  data.frame(transcript_id = ids,
             biotype   = biotype[ids],
             gene_id   = unname(g),
             gene_name = unname(ensg2name[g]),
             row.names = NULL)
}

compare <- function(sr_set, lr_set, tag) {
  both <- intersect(sr_set, lr_set); lr_only <- setdiff(lr_set, sr_set); sr_only <- setdiff(sr_set, lr_set)
  cat("\n== ", tag, " ==\n", sep = "")
  cat("SR:", length(sr_set), " LR:", length(lr_set), " both:", length(both),
      " LR-only:", length(lr_only), " SR-only:", length(sr_only), "\n")
  df <- annotate(lr_only)
  write.csv(df, file.path(out_dir, paste0("isoform_LRonly_", tag, ".csv")), row.names = FALSE)
  # gene-level summary of the LR-only isoforms
  gt <- as.data.frame(table(gene_id = df$gene_id), stringsAsFactors = FALSE)
  gt <- gt[order(-gt$Freq), ]
  gt$gene_name <- ensg2name[gt$gene_id]
  colnames(gt)[colnames(gt) == "Freq"] <- "n_LRonly_isoforms"
  write.csv(gt, file.path(out_dir, paste0("isoform_LRonly_", tag, "_by_gene.csv")), row.names = FALSE)
  cat("  unique genes with LR-only isoforms:", nrow(gt), "\n")
  cat("  top genes:", paste(head(paste0(gt$gene_name, "(", gt$n_LRonly_isoforms, ")"), 8), collapse = ", "), "\n")
  c(SR = length(sr_set), LR = length(lr_set), both = length(both),
    LR_only = length(lr_only), SR_only = length(sr_only), LRonly_genes = nrow(gt))
}

r1 <- compare(sr_det, lr_enst_det, "all_known")
pc <- names(biotype)[biotype == "protein_coding" & !is.na(biotype)]
r2 <- compare(intersect(sr_det, pc), intersect(lr_enst_det, pc), "protein_coding")

summ <- rbind(data.frame(comparison = "all_known", t(r1)),
              data.frame(comparison = "protein_coding", t(r2)))
summ$LR_novel_detected <- length(lr_novel_det)
write.csv(summ, file.path(out_dir, "isoform_SR_LR_overlap_summary.csv"), row.names = FALSE)
cat("\nSUMMARY\n"); print(summ); cat("\nDONE\n")
