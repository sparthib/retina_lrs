#!/usr/bin/env Rscript
# 08h_platform_DE_LRsalmon.R   (run: module load conda_R/4.5.x)
# Platform DE, QUANTIFIER-CONTROLLED: short-read salmon vs LONG-READ SALMON
# (alignment-mode, script 11) — replaces the SR-salmon vs LR-bambu contrast in
# 02_platform_DE so both platforms use the same quantifier. Method mirrors the
# committed stage DE: within-platform TMM, post-normalization CPM>=1 filter (both
# platforms), edgeR QLF on the platform factor. Significance FDR<0.05 AND
# |log2FC|>=1. LR-salmon subset to protein-coding to match bambu scope. Stage
# tables carry gene_name + per-gene median length/GC.
.libPaths("/users/sparthib/R/4.5.x")
suppressMessages({library(edgeR)})
strip <- function(x) sub("\\.[0-9]+$","",x)
dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
sal <- file.path(dd,"10_short_reads/salmon")
ref <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
outdir <- file.path(sal,"platform_DE_LRsalmon"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

sr <- as.matrix(read.csv(file.path(sal,"gene_counts_matrix.csv"), row.names=1, check.names=FALSE))
rownames(sr) <- strip(rownames(sr))
lr <- as.matrix(read.csv(file.path(sal,"LRsalmon_gene_counts_matrix.csv"), row.names=1, check.names=FALSE))
rownames(lr) <- strip(rownames(lr)); storage.mode(lr) <- "numeric"

gl <- system(paste0("awk -F'\t' '$3==\"gene\"' ", file.path(ref,"release_46_primary_assembly.gtf"),
      " | sed -E 's/.*gene_id \"([^\"]+)\".*gene_type \"([^\"]+)\".*/\\1\t\\2/'"), intern=TRUE)
gm <- do.call(rbind, strsplit(gl,"\t")); gbt <- setNames(gm[,2], strip(gm[,1]))
pc <- names(gbt)[gbt=="protein_coding"]
lr <- lr[rownames(lr) %in% pc,,drop=FALSE]
cat("SR:",nrow(sr),"x",ncol(sr),"  LR-salmon(PC):",nrow(lr),"x",ncol(lr),"\n")

hdr <- system(paste0("grep '^>' ", file.path(ref,"release_46_all_transcripts.fa"),
       " | awk -F'|' '{print $2\"\t\"$6}' | sort -u"), intern=TRUE)
hm <- do.call(rbind, strsplit(hdr,"\t")); ensg2sym <- setNames(hm[,2], strip(hm[,1]))
meta <- readRDS(file.path(ref,"transcript_meta_info.rds"))
t2g <- read.table(file.path(sal,"tx2gene_release46.tsv"),header=FALSE,sep="\t",stringsAsFactors=FALSE)
iso2gene <- setNames(strip(t2g$V2), strip(t2g$V1))
meta$gene <- iso2gene[strip(meta$isoform_id)]
g <- meta$gene; ok <- !is.na(g)
lenmap <- tapply(meta$transcript_length[ok], g[ok], median)
gcmap  <- round(tapply(meta$GC_percent[ok],      g[ok], median), 2)
nmap   <- tapply(meta$isoform_id[ok],           g[ok], length)
mxmap  <- tapply(meta$transcript_length[ok],    g[ok], max)

# LR stage from colname (hyphen or underscore tolerant)
lr_stage <- c("EP1_WT_ROs_D45"="S1","H9_CRX_ROs_D45"="S1",
              "EP1_WT_hRO_2"="S2","H9_BRN3B_hRO_2"="S2","H9_CRX_hRO_2"="S2",
              "EP1_BRN3B_RO"="S3","H9_BRN3B_RO"="S3")
lrkey <- gsub("-","_",colnames(lr))
lr_win <- list(S1=colnames(lr)[lr_stage[lrkey]=="S1"],
               S2=colnames(lr)[lr_stage[lrkey]=="S2"],
               S3=colnames(lr)[lr_stage[lrkey]=="S3"])
sr_day <- as.integer(sub("^D([0-9]+).*","\\1", colnames(sr)))
sr_win <- list(S1=colnames(sr)[sr_day %in% c(25,65)],
               S2=colnames(sr)[sr_day %in% c(65,100)],
               S3=colnames(sr)[sr_day %in% c(180,280)])
maplab <- c(S1="Stage1_LRd45_vs_SR_D25_D65", S2="Stage2_LRd100_vs_SR_D65_D100", S3="Stage3_LRd200_vs_SR_D180_D280")

summ <- data.frame()
for (w in c("S1","S2","S3")) {
  lrs <- lr_win[[w]]; srs <- sr_win[[w]]
  common <- intersect(rownames(lr), rownames(sr))
  lr_w <- lr[common, lrs, drop=FALSE]; sr_w <- sr[common, srs, drop=FALSE]
  cpm_lr <- cpm(lr_w); cpm_sr <- cpm(sr_w)
  keep <- (rowSums(cpm_lr>=1) >= length(lrs)) & (rowSums(cpm_sr>=1) >= length(lrs))
  lr_w <- lr_w[keep,,drop=FALSE]; sr_w <- sr_w[keep,,drop=FALSE]
  counts <- cbind(lr_w, sr_w)
  platform <- factor(c(rep("LR",ncol(lr_w)), rep("SR",ncol(sr_w))), levels=c("SR","LR"))
  dge <- DGEList(counts=counts, group=platform)
  d_lr <- calcNormFactors(DGEList(lr_w)); d_sr <- calcNormFactors(DGEList(sr_w))
  dge$samples$norm.factors <- c(d_lr$samples$norm.factors, d_sr$samples$norm.factors)
  design <- model.matrix(~platform)
  dge <- estimateDisp(dge, design); fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef="platformLR"); tt <- topTags(qlf, n=Inf)$table
  tt$gene <- rownames(tt)
  tt$gene_name <- ensg2sym[tt$gene]
  tt$median_tx_length <- lenmap[tt$gene]; tt$median_gc_pct <- gcmap[tt$gene]
  tt$n_transcripts <- nmap[tt$gene]; tt$max_tx_length <- mxmap[tt$gene]
  tt <- tt[,c("gene","gene_name","logFC","logCPM","F","PValue","FDR",
              "median_tx_length","median_gc_pct","n_transcripts","max_tx_length")]
  write.csv(tt, file.path(outdir, paste0("platform_DE_stage_W_",w,".csv")), row.names=FALSE)
  sig <- tt$FDR<0.05
  summ <- rbind(summ, data.frame(window=paste0("W_",w), mapping=maplab[w],
    n_LR=length(lrs), n_SR=length(srs), genes_tested=nrow(tt),
    sig_FDR05=sum(sig,na.rm=TRUE),
    sig_FDR05_absLFC1=sum(sig & abs(tt$logFC)>=1, na.rm=TRUE),
    up_in_LR_LFC1=sum(sig & tt$logFC>=1, na.rm=TRUE),
    up_in_SR_LFC1=sum(sig & tt$logFC<=-1, na.rm=TRUE)))
  cat(sprintf("W_%s: LR=%d SR=%d tested=%d sigFDR=%d sig+LFC=%d upLR=%d upSR=%d\n",
      w, length(lrs), length(srs), nrow(tt), sum(sig,na.rm=TRUE),
      sum(sig&abs(tt$logFC)>=1,na.rm=TRUE), sum(sig&tt$logFC>=1,na.rm=TRUE), sum(sig&tt$logFC<=-1,na.rm=TRUE)))
}
write.csv(summ, file.path(outdir,"platform_DE_stage_summary.csv"), row.names=FALSE)
cat("DONE\n"); print(summ, row.names=FALSE)
