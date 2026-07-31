#!/usr/bin/env Rscript
# 08f_platform_DE_GO_analysis.R
# (1) Annotate the platform-DE stage tables with gene symbols via biomaRt.
# (2) GO:BP over-representation of the platform-DE gene sets, per stage and direction.
#
# Gene sets (same significance rule as the DE call): FDR < 0.05 AND |log2FC| >= 1.
#   up_in_LR : logFC >= +1   (higher in long reads)
#   up_in_SR : logFC <= -1   (higher in short reads)
# Background universe for ORA = genes tested in that stage (rows of the DE table).
#
# Follows the lab GO convention (code/05_visualization/06_go_analysis.R + helper.R
# ora_plot): enrichGO(OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP",
# pAdjustMethod="fdr", minGSSize=100, pvalueCutoff=0.01, qvalueCutoff=0.01,
# readable=TRUE); dotplot(showCategory=15).
#
# Run with: module load conda_R/4.4.x && Rscript 08f_platform_DE_GO_analysis.R

.libPaths(c("/users/sparthib/R/4.4.x", .libPaths()))
suppressMessages({library(biomaRt); library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})

de_dir <- "/users/sparthib/retina_lrs/processed_data/short_reads/02_platform_DE"
out_dir <- file.path(de_dir, "go")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stages <- c("W_S1", "W_S2", "W_S3")
strip <- function(x) sub("\\.[0-9]+$", "", x)

## ---- gene symbol map via biomaRt (all stage genes at once) ----
files <- setNames(file.path(de_dir, paste0("platform_DE_stage_", stages, ".csv")), stages)
tabs  <- lapply(files, read.csv, stringsAsFactors = FALSE)
all_ensg <- unique(strip(unlist(lapply(tabs, function(d) d$gene))))
cat("unique genes across stages:", length(all_ensg), "\n")

sym <- NULL
for (host in c("https://www.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")) {
  sym <- tryCatch({
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = host)
    getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
          filters = "ensembl_gene_id", values = all_ensg, mart = mart)
  }, error = function(e) { message("biomaRt via ", host, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(sym) && nrow(sym) > 0) { cat("biomaRt OK via", host, "-> mapped", nrow(sym), "genes\n"); break }
}
if (is.null(sym) || nrow(sym) == 0) stop("biomaRt returned no mappings from any mirror")
g2s <- setNames(sym$external_gene_name, sym$ensembl_gene_id)

## ---- annotate + rewrite each DE table with gene_name (2nd column) ----
for (s in stages) {
  d <- tabs[[s]]
  d$gene_name <- unname(g2s[strip(d$gene)])
  d <- d[, c("gene", "gene_name", setdiff(colnames(d), c("gene", "gene_name")))]
  write.csv(d, files[s], row.names = FALSE)
  cat(s, ": annotated", sum(!is.na(d$gene_name) & d$gene_name != ""), "/", nrow(d), "genes with a symbol\n")
  tabs[[s]] <- d
}

## ---- GO ORA per stage x direction ----
run_go <- function(genes, universe, tag) {
  cat("\n==", tag, ":", length(genes), "genes ==\n")
  if (length(genes) < 10) { cat("  too few genes, skipping\n"); return(invisible(NULL)) }
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "fdr",
                  minGSSize = 100, pvalueCutoff = 0.01, qvalueCutoff = 0.01, readable = TRUE)
  ed <- as.data.frame(ego)
  cat("  enriched BP terms:", nrow(ed), "\n")
  if (nrow(ed) > 0) {
    write.csv(ed, file.path(out_dir, paste0("GO_BP_", tag, ".csv")), row.names = FALSE)
    p <- dotplot(ego, showCategory = 15) + ggtitle(paste0("GO:BP  ", tag))
    ggsave(file.path(out_dir, paste0("GO_BP_", tag, ".png")), p, width = 9, height = 7, dpi = 200)
    cat("  top:", paste(head(ed$Description, 5), collapse = " | "), "\n")
  }
  invisible(ed)
}

for (s in stages) {
  d <- tabs[[s]]
  universe <- strip(d$gene)
  up_lr <- strip(d$gene[d$FDR < 0.05 & d$logFC >=  1])
  up_sr <- strip(d$gene[d$FDR < 0.05 & d$logFC <= -1])
  run_go(up_lr, universe, paste0(s, "_up_in_LR"))
  run_go(up_sr, universe, paste0(s, "_up_in_SR"))
}
cat("\nDONE\n")
