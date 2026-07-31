#!/usr/bin/env Rscript
# 09_isoforms_per_gene_SR_LR.R
# Unique number of DETECTED isoforms per gene, for short reads (salmon) and long
# reads (bambu), and a side-by-side comparison.
#
# Detection rule (both platforms): >= 1 estimated count in >= 1 sample.
# Only genes with >= 1 detected isoform are counted. Transcript->gene mapping:
#   SR : release-46 tx2gene (ENST -> ENSG); LR : bambu counts_transcript.txt (TXNAME -> GENEID).
# IDs version-stripped. LR includes novel (Bambu) isoforms in its per-gene counts.
# Run with: module load conda_R/4.5.x && Rscript 09_isoforms_per_gene_SR_LR.R

.libPaths("/users/sparthib/R/4.5.x")
dd <- "/dcs04/hicks/data/sparthib/retina_lrs"
sr_file  <- file.path(dd, "10_short_reads/salmon/transcript_counts_matrix.csv")
tx2g_sr  <- file.path(dd, "10_short_reads/salmon/tx2gene_release46.tsv")
lr_file  <- file.path(dd, "06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS")
bambu_tx <- file.path(dd, "06_quantification/bambu/all_samples_extended_annotation_track_reads/counts_transcript.txt")
out_dir  <- file.path(dd, "10_short_reads/salmon")
strip <- function(x) sub("\\.[0-9]+$", "", x)

summarize <- function(n_per_gene, label) {
  q <- quantile(n_per_gene, c(.25,.5,.75,.9,.95))
  data.frame(platform=label, n_genes=length(n_per_gene),
             min=min(n_per_gene), median=median(n_per_gene),
             mean=round(mean(n_per_gene),3), q75=unname(q[3]), q90=unname(q[4]),
             q95=unname(q[5]), max=max(n_per_gene),
             pct_1iso=round(100*mean(n_per_gene==1),1),
             n_ge2=sum(n_per_gene>=2), n_ge5=sum(n_per_gene>=5), n_ge10=sum(n_per_gene>=10))
}

## ---- short read ----
sr <- as.matrix(read.csv(sr_file, row.names=1, check.names=FALSE))
rownames(sr) <- strip(rownames(sr))
sr_det <- rownames(sr)[rowSums(sr >= 1) >= 1]
t2g <- read.table(tx2g_sr, header=FALSE, sep="\t", stringsAsFactors=FALSE)
sr_map <- setNames(strip(t2g$V2), strip(t2g$V1))
sr_genes <- sr_map[sr_det]; sr_genes <- sr_genes[!is.na(sr_genes)]
sr_per_gene <- as.integer(table(sr_genes))
cat("SR detected isoforms:", length(sr_det), " genes:", length(sr_per_gene), "\n")

## ---- long read ----
lr <- as.matrix(readRDS(lr_file)); rownames(lr) <- strip(rownames(lr))
lr_det <- rownames(lr)[rowSums(lr >= 1) >= 1]
b <- read.table(bambu_tx, header=TRUE, sep="\t", comment.char="",
                colClasses=c("character","character", rep("NULL",11)))
lr_map <- setNames(strip(b$GENEID), strip(b$TXNAME))
lr_genes <- lr_map[lr_det]; lr_genes <- lr_genes[!is.na(lr_genes)]
lr_per_gene <- as.integer(table(lr_genes))
cat("LR detected isoforms:", length(lr_det), " genes:", length(lr_per_gene), "\n")

## ---- summary table ----
summ <- rbind(summarize(sr_per_gene,"short_read"), summarize(lr_per_gene,"long_read"))
write.csv(summ, file.path(out_dir,"isoforms_per_gene_SR_LR_summary.csv"), row.names=FALSE)
cat("\nSUMMARY\n"); print(summ)

## ---- distribution table (count of genes at each isoform number, capped bins) ----
bin <- function(v){ b<-pmin(v,20); t<-table(factor(b, levels=1:20)); as.integer(t) }
dist <- data.frame(n_isoforms=c(as.character(1:19),"20+"),
                   short_read=bin(sr_per_gene), long_read=bin(lr_per_gene))
write.csv(dist, file.path(out_dir,"isoforms_per_gene_SR_LR_distribution.csv"), row.names=FALSE)

## ---- comparison plot (overlaid histograms, capped at 20) ----
library(ggplot2)
df <- rbind(data.frame(platform="short read", n=pmin(sr_per_gene,20)),
            data.frame(platform="long read",  n=pmin(lr_per_gene,20)))
p <- ggplot(df, aes(x=n, fill=platform)) +
  geom_histogram(position="identity", alpha=0.55, binwidth=1, colour="grey30") +
  scale_fill_manual(values=c("short read"="#D7301F","long read"="#2166AC")) +
  scale_x_continuous(breaks=c(1,5,10,15,20), labels=c("1","5","10","15","20+")) +
  labs(x="Detected isoforms per gene (capped at 20+)", y="Number of genes",
       title="Detected isoforms per gene: short read (salmon) vs long read (bambu)",
       fill=NULL) +
  theme_bw(base_size=13) + theme(legend.position=c(0.85,0.85))
ggsave(file.path(out_dir,"isoforms_per_gene_SR_LR.png"), p, width=8, height=5.5, dpi=200)
cat("\nDONE\n")
