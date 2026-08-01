#!/usr/bin/env Rscript
# 08i_correlation_and_lengthgc_LRsalmon.R  (module load conda_R/4.5.x)
# QUANTIFIER-CONTROLLED redo of 01_correlation and 03_length_gc using LR-SALMON.
#  (A) Correlation: SR-salmon vs LR-salmon gene CPM over all common protein-coding
#      genes (Spearman), heatmap with the committed white->pink->red scheme.
#  (B) Length/GC: distributions of median gene GC% and tx length for genes up- vs
#      down-in-LR at each platform-DE stage, from the LR-salmon DE tables (08h).
.libPaths("/users/sparthib/R/4.5.x")
suppressMessages({library(edgeR); library(ggplot2); library(reshape2); library(dplyr); library(patchwork)})
strip <- function(x) sub("\\.[0-9]+$","",x)
dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
sal <- file.path(dd,"10_short_reads/salmon")
ref <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
stg <- file.path(sal,"platform_DE_LRsalmon")
out <- file.path(sal,"corr_lengthgc_LRsalmon"); dir.create(out, showWarnings=FALSE)

## protein-coding gene ids
gl <- system(paste0("awk -F'\t' '$3==\"gene\"' ", file.path(ref,"release_46_primary_assembly.gtf"),
      " | sed -E 's/.*gene_id \"([^\"]+)\".*gene_type \"([^\"]+)\".*/\\1\t\\2/'"), intern=TRUE)
gm <- do.call(rbind, strsplit(gl,"\t")); pc <- strip(gm[gm[,2]=="protein_coding",1])

## ---- (A) correlation ----
sr <- as.matrix(read.csv(file.path(sal,"gene_counts_matrix.csv"), row.names=1, check.names=FALSE))
rownames(sr) <- strip(rownames(sr))
lr <- as.matrix(read.csv(file.path(sal,"LRsalmon_gene_counts_matrix.csv"), row.names=1, check.names=FALSE))
rownames(lr) <- strip(rownames(lr)); storage.mode(lr) <- "numeric"
sr_cpm <- cpm(DGEList(sr)); lr_cpm <- cpm(DGEList(lr))
common <- Reduce(intersect, list(rownames(sr_cpm), rownames(lr_cpm), pc))
common <- setdiff(common, "ENSG00000210082")   # drop MT-rRNA outlier (as committed)
cat("common protein-coding genes for correlation:", length(common), "\n")
srC <- sr_cpm[common,,drop=FALSE]; lrC <- lr_cpm[common,,drop=FALSE]
cormat <- cor(srC, lrC, method="spearman")     # rows SR, cols LR
write.csv(cormat, file.path(out,"SRsalmon_vs_LRsalmon_spearman_corr_allPC.csv"))
cat("Spearman range:", round(range(cormat),3), "\n")

cd <- melt(cormat); colnames(cd) <- c("SR","LR","value")
sr_ord <- c("D0_1","D0_2","D0_3","D10_3","D10_4","D10_5","D10_6","D25_1","D25_2","D25_3","D25_4",
            "D65_2","D65_3","D65_5","D65_6","D100_3","D100_4","D100_5","D100_6",
            "D180_1","D180_2","D180_3","D180_4","D180_5","D180_6","D280_A2","D280_B1","D280_C1")
cd$SR <- factor(cd$SR, levels=intersect(sr_ord, unique(as.character(cd$SR))))
p <- ggplot(cd, aes(SR, LR, fill=value)) + geom_tile() +
  geom_text(aes(label=sprintf("%.2f", value)), color="black", size=2.6) +
  scale_fill_gradient2(low="white", mid="pink", high="red", midpoint=0.5,
                       limits=c(0,1), name="Spearman") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Short-read (salmon) samples", y="Long-read (salmon) samples",
       title="SR-salmon vs LR-salmon gene-expression correlation (all common protein-coding genes)")
ggsave(file.path(out,"SRsalmon_vs_LRsalmon_correlation_heatmap.png"), p, width=11, height=5.5, dpi=200)
cat("SAVED correlation heatmap\n")

## ---- (B) length / GC of platform-DE genes (LR-salmon DE) ----
stages <- c(W_S1="Stage 1 (LR D45 vs SR D25+D65)",
            W_S2="Stage 2 (LR D100 vs SR D65+D100)",
            W_S3="Stage 3 (LR D200 vs SR D180+D280)")
LFC<-1; FDRc<-0.05; all_df<-list(); summ<-data.frame()
for (w in names(stages)) {
  de <- read.csv(file.path(stg, paste0("platform_DE_stage_", w, ".csv")), stringsAsFactors=FALSE)
  de$direction <- with(de, ifelse(FDR<FDRc & logFC>= LFC, "Up in long read",
                            ifelse(FDR<FDRc & logFC<=-LFC, "Down in long read","Not sig")))
  de$stage <- stages[[w]]
  keep <- de[de$direction!="Not sig" & !is.na(de$median_gc_pct), ]
  all_df[[w]] <- keep
  for (d in c("Up in long read","Down in long read")) {
    s <- keep[keep$direction==d, ]
    summ <- rbind(summ, data.frame(stage=stages[[w]], direction=d, n=nrow(s),
                  median_gc=round(median(s$median_gc_pct),2),
                  median_len=round(median(s$median_tx_length),1)))
  }
}
df <- bind_rows(all_df)
df$direction <- factor(df$direction, levels=c("Up in long read","Down in long read"))
df$stage <- factor(df$stage, levels=unname(stages))
write.csv(summ, file.path(out,"length_gc_distribution_summary.csv"), row.names=FALSE)
# also write the per-gene table used
write.csv(df[,c("gene","gene_name","stage","direction","median_gc_pct","median_tx_length","logFC","FDR")],
          file.path(out,"gene_length_gc.csv"), row.names=FALSE)
cat("length/GC summary:\n"); print(summ)
cols <- c("Up in long read"="#D7301F","Down in long read"="#2166AC")
p_gc <- ggplot(df, aes(median_gc_pct, fill=direction, colour=direction)) +
  geom_density(alpha=0.35, linewidth=0.5) + facet_wrap(~stage, nrow=1) +
  scale_fill_manual(values=cols)+scale_colour_manual(values=cols) +
  labs(x="Median gene GC content (%)", y="Density", fill=NULL, colour=NULL,
       title="A. GC content of platform-DE genes (SR-salmon vs LR-salmon, |logFC|>=1, FDR<0.05)") +
  theme_bw(base_size=10)+theme(legend.position="top", panel.grid.minor=element_blank())
p_len <- ggplot(df, aes(median_tx_length, fill=direction, colour=direction)) +
  geom_density(alpha=0.35, linewidth=0.5) + facet_wrap(~stage, nrow=1) + scale_x_log10() +
  scale_fill_manual(values=cols)+scale_colour_manual(values=cols) +
  labs(x="Median transcript length (bp, log10)", y="Density", fill=NULL, colour=NULL,
       title="B. Transcript length of platform-DE genes (SR-salmon vs LR-salmon, |logFC|>=1, FDR<0.05)") +
  theme_bw(base_size=10)+theme(legend.position="top", panel.grid.minor=element_blank())
ggsave(file.path(out,"length_gc_distributions.png"), p_gc/p_len, width=11, height=7, dpi=200)
cat("SAVED length/gc plot\nDONE\n")
