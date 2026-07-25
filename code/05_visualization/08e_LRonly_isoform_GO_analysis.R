#!/usr/bin/env Rscript
# 08e_LRonly_isoform_GO_analysis.R
# GO over-representation analysis (Biological Process) of genes bearing isoforms
# detected ONLY in long reads (not detected in the short-read salmon quant).
#
# Two gene sets (from 08_isoform_SR_LR_comparison.R):
#   - all_known      : genes with any LR-only known (ENST) isoform (13,801 genes)
#   - protein_coding : genes with LR-only protein-coding isoforms  (2,776 genes)
#
# Follows the lab GO convention (code/05_visualization/06_go_analysis.R + helper.R):
#   enrichGO(OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP", pAdjustMethod="fdr",
#            minGSSize=100, pvalueCutoff=0.01, qvalueCutoff=0.01, readable=TRUE)
#   dotplot(ego, showCategory=15).
#
# Background (universe): genes DETECTED in long reads (the pool LR-only is drawn from),
# so enrichment reflects "LR-only vs all LR-detected", not vs the whole genome.
# Run with: module load conda_R/4.4.x && Rscript 08e_LRonly_isoform_GO_analysis.R

.libPaths(c("/users/sparthib/R/4.4.x", .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})

data_dir <- "/dcs04/hicks/data/sparthib/retina_lrs"
iso_dir  <- file.path(data_dir, "10_short_reads/salmon")
lr_file  <- file.path(data_dir, "06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS")
bambu_tx <- file.path(data_dir, "06_quantification/bambu/all_samples_extended_annotation_track_reads/counts_transcript.txt")
out_dir  <- file.path(iso_dir, "go")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
strip <- function(x) sub("\\.[0-9]+$", "", x)

## universe = genes with any DETECTED known isoform in LR
b <- read.table(bambu_tx, header = TRUE, sep = "\t", comment.char = "",
                colClasses = c("character","character", rep("NULL", 11)))
tx2g <- setNames(strip(b$GENEID), strip(b$TXNAME))
lr <- as.matrix(readRDS(lr_file)); rownames(lr) <- strip(rownames(lr))
lr_det_enst <- rownames(lr)[grepl("^ENST", rownames(lr)) & rowSums(lr >= 1) >= 1]
universe <- unique(tx2g[lr_det_enst]); universe <- universe[grepl("^ENSG", universe)]
cat("universe (LR-detected genes):", length(universe), "\n")

run_go <- function(by_gene_csv, tag) {
  gt <- read.csv(by_gene_csv, stringsAsFactors = FALSE)
  genes <- unique(gt$gene_id[grepl("^ENSG", gt$gene_id)])
  cat("\n==", tag, ": ", length(genes), "genes ==\n")
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "fdr",
                  minGSSize = 100, pvalueCutoff = 0.01, qvalueCutoff = 0.01,
                  readable = TRUE)
  ed <- as.data.frame(ego)
  cat("enriched BP terms:", nrow(ed), "\n")
  if (nrow(ed) > 0) {
    write.csv(ed, file.path(out_dir, paste0("GO_BP_LRonly_", tag, ".csv")), row.names = FALSE)
    p <- dotplot(ego, showCategory = 15) +
      ggtitle(paste0("GO:BP — genes with LR-only ", gsub("_", " ", tag), " isoforms"))
    ggsave(file.path(out_dir, paste0("GO_BP_LRonly_", tag, ".png")), p,
           width = 9, height = 7, dpi = 200)
    cat("top terms:", paste(head(ed$Description, 6), collapse = " | "), "\n")
  }
  invisible(ed)
}

run_go(file.path(iso_dir, "isoform_LRonly_all_known_by_gene.csv"),      "all_known")
run_go(file.path(iso_dir, "isoform_LRonly_protein_coding_by_gene.csv"), "protein_coding")
cat("\nDONE\n")
