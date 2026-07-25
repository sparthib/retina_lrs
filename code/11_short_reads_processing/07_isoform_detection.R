#!/usr/bin/env Rscript
# 07_isoform_detection.R
# Short-read isoform detection and per-gene isoform diversity.
#
# Input : transcript-level salmon counts matrix (tximport, txOut=TRUE),
#         251,955 GENCODE release-46 transcripts x 28 samples, and the
#         transcript-to-gene map used to build the salmon index.
# Output: dataset-, sample-, and timepoint-level counts of detected isoforms,
#         and the distribution of detected isoforms per gene.
#
# "Detected" = an annotated transcript with expression evidence. This is
# short-read isoform *quantification* against the annotation, not isoform
# *discovery*; multi-mapping reads are distributed across a gene's annotated
# isoforms. The confident set (>=5 counts in >=3 samples) is the more
# conservative view. Run with: module load conda_R/4.5.x
#
# Run: module load conda_R/4.5.x && Rscript 07_isoform_detection.R

.libPaths("/users/sparthib/R/4.5.x")
suppressMessages(library(ggplot2))

salmon_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon"
tx_file    <- file.path(salmon_dir, "transcript_counts_matrix.csv")
t2g_file   <- file.path(salmon_dir, "tx2gene_release46.tsv")
out_dir    <- salmon_dir   # csv/png outputs (copied into the repo separately)

tx  <- as.matrix(read.csv(tx_file, row.names = 1, check.names = FALSE))
t2g <- read.table(t2g_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
strip <- function(x) sub("\\.[0-9]+$", "", x)
rownames(tx) <- strip(rownames(tx))
t2g$V1 <- strip(t2g$V1); t2g$V2 <- strip(t2g$V2)
gene <- setNames(t2g$V2, t2g$V1)[rownames(tx)]

## ---- 1. dataset / sample / timepoint detection ----
samp <- colnames(tx); day <- sub("_.*", "", samp)
cat("matrix:", nrow(tx), "transcripts x", ncol(tx), "samples\n")
cat("detected (>0 in >=1 sample):", sum(rowSums(tx > 0)  >= 1), "\n")
cat(">=1 count in >=1 sample:    ", sum(rowSums(tx >= 1) >= 1), "\n")
cat(">=1 count in >=3 samples:   ", sum(rowSums(tx >= 1) >= 3), "\n")
cat(">=5 count in >=3 samples:   ", sum(rowSums(tx >= 5) >= 3), "\n")

per_samp <- colSums(tx >= 1)
write.csv(data.frame(sample = samp, day = day, isoforms_detected_ge1 = per_samp,
                     row.names = NULL),
          file.path(out_dir, "transcript_detection_per_sample.csv"), row.names = FALSE)

days <- unique(day)
day_union <- sapply(days, function(d) sum(rowSums(tx[, day == d, drop = FALSE] >= 1) >= 1))
write.csv(data.frame(day = days, isoforms_detected_union_ge1 = as.integer(day_union),
                     row.names = NULL),
          file.path(out_dir, "transcript_detection_per_day.csv"), row.names = FALSE)

## ---- 2. isoforms detected per gene ----
det  <- rowSums(tx >= 1) >= 1              # detected: >=1 count in >=1 sample
detC <- rowSums(tx >= 5) >= 3              # confident: >=5 counts in >=3 samples
g_det <- table(gene[det]); g_ann <- table(gene[!is.na(gene)])
v <- as.integer(g_det); vc <- as.integer(table(gene[detC]))

cat("\n== isoforms per gene (detected set) ==\n")
cat("n genes:", length(v), " median:", median(v), " mean:", round(mean(v), 2),
    " max:", max(v), "\n")
cat("genes with exactly 1 detected isoform:", sum(v == 1),
    "frac", round(mean(v == 1), 3), "\n")
cat("confident set: n genes", length(vc), " median", median(vc), "\n")

per_gene <- data.frame(gene = names(g_det),
                       n_detected_iso  = as.integer(g_det),
                       n_annotated_iso = as.integer(g_ann[names(g_det)]),
                       row.names = NULL)
per_gene$frac_detected <- round(per_gene$n_detected_iso / per_gene$n_annotated_iso, 3)
write.csv(per_gene, file.path(out_dir, "isoforms_per_gene_detected.csv"), row.names = FALSE)

tabd <- as.data.frame(table(n_iso = v)); colnames(tabd) <- c("n_detected_isoforms", "n_genes")
write.csv(tabd, file.path(out_dir, "isoforms_per_gene_distribution.csv"), row.names = FALSE)

dfp <- data.frame(n = v); dfp$ncap <- pmin(dfp$n, 15)
sub_txt <- paste0(length(v), " genes with >=1 detected isoform; median=",
                  median(v), ", mean=", round(mean(v), 2))
p1 <- ggplot(dfp, aes(x = ncap)) +
  geom_bar(fill = "#3182BD", color = "white") +
  scale_x_continuous(breaks = 1:15, labels = c(as.character(1:14), "15+")) +
  labs(x = "Detected isoforms per gene", y = "Number of genes",
       title = "Distribution of detected isoforms per gene", subtitle = sub_txt) +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "isoforms_per_gene_distribution.png"), p1,
       width = 8, height = 5, dpi = 200)
cat("\nDONE\n")
