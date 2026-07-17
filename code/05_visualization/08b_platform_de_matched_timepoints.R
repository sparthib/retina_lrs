#!/usr/bin/env Rscript
# 08b_platform_de_matched_timepoints.R
# Run with: module load conda_R/4.5.x
#
# PLATFORM (technology) differential expression: long-read (direct cDNA / bambu)
# vs short-read (Illumina / salmon) at matched developmental windows.
#
# IMPORTANT INTERPRETATION: platform is fully confounded with technology, so this
# is NOT a biological DE. It cannot be batch-corrected (removing the platform
# effect removes the entire comparison). Results describe TECHNICAL
# concordance/discordance between the two platforms, not biology. The very high
# fraction of significant genes (55-68%) is the expected signature of a
# technology comparison (length bias, molecule-vs-fragment, 3' bias).
#
# Design (user rulings):
#   LR D45  vs SR D0 + D10 + D25   (window W_D45)
#   LR D100 vs SR D65 + D100       (window W_D100; D100 is the one exact match)
#   LR D200 vs SR D180 + D280      (window W_D200)
# Method: within-platform TMM normalization; post-normalization CPM>=1 filter
# (depth-independent, applied in BOTH platforms); edgeR quasi-likelihood F-test
# on the platform factor. Gene universe = LR protein-coding + expression-filtered
# set intersected with SR genes; ENSG version suffixes stripped on both sides.
#

.libPaths("/users/sparthib/R/4.5.x")
suppressMessages({library(edgeR)})

sr_path <- "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/salmon/gene_counts_matrix.csv"
lr_path <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype/filtered_gene_counts.RDS"
outdir  <- "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/08b_platform_de"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

## ---- load raw counts, strip ENSG versions on both sides ----
sr <- as.matrix(read.csv(sr_path, row.names=1, check.names=FALSE))
rownames(sr) <- sub("\\.[0-9]+$","", rownames(sr))
lr <- as.matrix(readRDS(lr_path)); storage.mode(lr) <- "numeric"
rownames(lr) <- sub("\\.[0-9]+$","", rownames(lr))
cat("SR:", nrow(sr),"x",ncol(sr),"  LR:", nrow(lr),"x",ncol(lr),"\n")

lr_tp <- c(H9_CRX_ROs_D45="D45", EP1_WT_ROs_D45="D45",
           EP1_WT_hRO_2="D100", H9_BRN3B_hRO_2="D100", H9_CRX_hRO_2="D100",
           EP1_BRN3B_RO="D200", H9_BRN3B_RO="D200")
stopifnot(all(colnames(lr) %in% names(lr_tp)))

sr_day <- as.integer(sub("^D([0-9]+).*","\\1", colnames(sr)))
sr_win <- ifelse(sr_day %in% c(0,10,25), "W_D45",
          ifelse(sr_day %in% c(65,100), "W_D100",
          ifelse(sr_day %in% c(180,280), "W_D200", NA)))
names(sr_win) <- colnames(sr)

windows <- list(
  W_D45  = list(lr=names(lr_tp)[lr_tp=="D45"],  sr=names(sr_win)[sr_win=="W_D45"]),
  W_D100 = list(lr=names(lr_tp)[lr_tp=="D100"], sr=names(sr_win)[sr_win=="W_D100"]),
  W_D200 = list(lr=names(lr_tp)[lr_tp=="D200"], sr=names(sr_win)[sr_win=="W_D200"])
)

summ <- data.frame()
for (w in names(windows)) {
  lrs <- windows[[w]]$lr; srs <- windows[[w]]$sr
  common <- intersect(rownames(lr), rownames(sr))
  lr_w <- lr[common, lrs, drop=FALSE]; sr_w <- sr[common, srs, drop=FALSE]

  cpm_lr <- cpm(lr_w); cpm_sr <- cpm(sr_w)
  keep <- (rowSums(cpm_lr >= 1) >= length(lrs)) & (rowSums(cpm_sr >= 1) >= length(lrs))
  lr_w <- lr_w[keep,,drop=FALSE]; sr_w <- sr_w[keep,,drop=FALSE]

  counts <- cbind(lr_w, sr_w)
  platform <- factor(c(rep("LR",ncol(lr_w)), rep("SR",ncol(sr_w))), levels=c("SR","LR"))
  dge <- DGEList(counts=counts, group=platform)
  d_lr <- calcNormFactors(DGEList(lr_w)); d_sr <- calcNormFactors(DGEList(sr_w))
  dge$samples$norm.factors <- c(d_lr$samples$norm.factors, d_sr$samples$norm.factors)

  design <- model.matrix(~platform)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef="platformLR")
  tt <- topTags(qlf, n=Inf)$table
  tt$gene <- rownames(tt)
  tt <- tt[,c("gene","logFC","logCPM","F","PValue","FDR")]
  write.csv(tt, file.path(outdir, paste0("08b_platform_DE_", w, ".csv")), row.names=FALSE)

  n_sig <- sum(tt$FDR < 0.05, na.rm=TRUE)
  summ <- rbind(summ, data.frame(window=w, n_LR=length(lrs), n_SR=length(srs),
                genes_tested=nrow(tt), sig_FDR05=n_sig,
                up_in_LR=sum(tt$FDR<0.05 & tt$logFC>0, na.rm=TRUE),
                up_in_SR=sum(tt$FDR<0.05 & tt$logFC<0, na.rm=TRUE)))
  cat(sprintf("%s: LR=%d SR=%d tested=%d sig=%d\n", w, length(lrs), length(srs), nrow(tt), n_sig))
}
write.csv(summ, file.path(outdir, "08b_platform_DE_summary.csv"), row.names=FALSE)
cat("DONE\n"); print(summ)
