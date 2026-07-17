#!/usr/bin/env Rscript
# Run with: module load conda_R/4.5.x
# Aggregate 28 salmon quant.sf (transcript-level) -> gene-level raw counts.
# ENSG version suffixes stripped to match the long-read bambu gene matrix.
# Column names use underscores (D100_5, D280_A2, ...) and are ordered by timepoint.
.libPaths("/users/sparthib/R/4.5.x")
suppressMessages(library(tximport))
OUT <- "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon"
quant_dir <- file.path(OUT, "quants")
tx2gene <- read.table(file.path(OUT, "tx2gene_release46.tsv"),
                      header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(tx2gene) <- c("TXNAME","GENEID")
tx2gene$GENEID <- gsub("\\.[0-9]+$", "", tx2gene$GENEID)   # strip ENSG version

samples <- list.dirs(quant_dir, recursive=FALSE, full.names=FALSE)
files <- file.path(quant_dir, samples, "quant.sf"); names(files) <- samples
stopifnot(all(file.exists(files)))

txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=FALSE)
cts <- round(txi$counts)
colnames(cts) <- gsub("-", "_", colnames(cts))            # D100-5 -> D100_5
day <- as.integer(sub("^D([0-9]+).*", "\\1", colnames(cts)))
cts <- cts[, order(day, colnames(cts))]

cat("matrix:", nrow(cts), "genes x", ncol(cts), "samples\n")
write.csv(cts, file.path(OUT, "gene_counts_matrix.csv"))
saveRDS(txi, file.path(OUT, "txi_salmon.rds"))


# --- Transcript-level counts (txOut=TRUE, no gene aggregation) ---
txi_tx <- tximport(files, type="salmon", txOut=TRUE)
tx <- round(txi_tx$counts)
colnames(tx) <- gsub("-", "_", colnames(tx))
tx <- tx[, order(as.integer(sub("^D([0-9]+).*","\\1", colnames(tx))), colnames(tx))]
rownames(tx) <- gsub("\\.[0-9]+$", "", rownames(tx))   # bare ENST
write.csv(tx, file.path(OUT, "transcript_counts_matrix.csv"))
saveRDS(txi_tx, file.path(OUT, "txi_salmon_transcript.rds"))
