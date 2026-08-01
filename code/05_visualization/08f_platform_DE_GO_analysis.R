#!/usr/bin/env Rscript
# 08f_platform_DE_GO_analysis.R
# GO:BP over-representation of the platform-DE gene sets, per stage and direction.
# NOTE: as of the LR-salmon redo, the DE tables (02_platform_DE/platform_DE_stage_*.csv)
# already carry gene_name (from the release-46 FASTA header map in 08h), so this
# script no longer calls biomaRt — it reads the existing tables directly.
#
# Gene sets (same rule as the DE call): FDR < 0.05 AND |log2FC| >= 1.
#   up_in_LR : logFC >= +1   (higher in long reads, salmon)
#   up_in_SR : logFC <= -1   (higher in short reads, salmon)
# Background universe = genes tested in that stage. Lab GO convention:
# enrichGO(ENSEMBL/BP/fdr, minGSSize=100, p/q<0.01, readable), dotplot n=15.
# Run: module load conda_R/4.4.x
.libPaths(c("/users/sparthib/R/4.4.x", .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
de_dir <- "/users/sparthib/retina_lrs/processed_data/short_reads/02_platform_DE"
out_dir <- file.path(de_dir, "go"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
stages <- c("W_S1","W_S2","W_S3"); strip <- function(x) sub("\\.[0-9]+$","",x)
files <- setNames(file.path(de_dir, paste0("platform_DE_stage_", stages, ".csv")), stages)
tabs  <- lapply(files, read.csv, stringsAsFactors=FALSE)

run_go <- function(genes, universe, tag) {
  cat("\n==", tag, ":", length(genes), "genes ==\n")
  if (length(genes) < 10) { cat("  too few, skip\n"); return(invisible()) }
  ego <- enrichGO(gene=genes, universe=universe, OrgDb=org.Hs.eg.db, keyType="ENSEMBL",
                  ont="BP", pAdjustMethod="fdr", minGSSize=100, pvalueCutoff=0.01,
                  qvalueCutoff=0.01, readable=TRUE)
  ed <- as.data.frame(ego); cat("  BP terms:", nrow(ed), "\n")
  if (nrow(ed) > 0) {
    write.csv(ed, file.path(out_dir, paste0("GO_BP_", tag, ".csv")), row.names=FALSE)
    p <- dotplot(ego, showCategory=15) + ggtitle(paste0("GO:BP  ", tag))
    ggsave(file.path(out_dir, paste0("GO_BP_", tag, ".png")), p, width=9, height=7, dpi=200)
    cat("  top:", paste(head(ed$Description,5), collapse=" | "), "\n")
  }
  invisible(ed)
}
for (s in stages) {
  d <- tabs[[s]]; universe <- strip(d$gene)
  up_lr <- strip(d$gene[d$FDR<0.05 & d$logFC>= 1])
  up_sr <- strip(d$gene[d$FDR<0.05 & d$logFC<=-1])
  run_go(up_lr, universe, paste0(s,"_up_in_LR"))
  run_go(up_sr, universe, paste0(s,"_up_in_SR"))
}
cat("\nDONE\n")
