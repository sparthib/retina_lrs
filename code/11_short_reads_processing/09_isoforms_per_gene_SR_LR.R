#!/usr/bin/env Rscript
# 09_isoforms_per_gene_SR_LR.R
# Unique number of DETECTED isoforms per gene, short reads (salmon) vs long reads
# (bambu), for (A) all genes and (B) protein-coding genes only.
#
# Detection rule (both platforms): >= 1 estimated count in >= 1 sample.
# Protein-coding gene = a gene with >= 1 protein_coding transcript in the GENCODE
# release-46 annotation (transcript_meta_info.rds). For those genes ALL detected
# isoforms are counted (a coding gene may also express non-coding isoforms), and on
# the LR side that includes novel (Bambu) isoforms belonging to the gene.
# Transcript->gene: SR = release-46 tx2gene (ENST->ENSG); LR = bambu
# counts_transcript.txt (TXNAME->GENEID). IDs version-stripped.
# Run with: module load conda_R/4.5.x && Rscript 09_isoforms_per_gene_SR_LR.R

.libPaths("/users/sparthib/R/4.5.x")
dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
ref <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
sr_file  <- file.path(dd, "10_short_reads/salmon/transcript_counts_matrix.csv")
tx2g_sr  <- file.path(dd, "10_short_reads/salmon/tx2gene_release46.tsv")
lr_file  <- file.path(dd, "06_quantification/counts_matrices/bambu/ROs/filtered/isoform_counts.RDS")
bambu_tx <- file.path(dd, "06_quantification/bambu/all_samples_extended_annotation_track_reads/counts_transcript.txt")
meta_f   <- file.path(ref, "transcript_meta_info.rds")
out_dir  <- file.path(dd, "10_short_reads/salmon")
strip <- function(x) sub("\\.[0-9]+$", "", x)

## ---- transcript->gene maps ----
t2g <- read.table(tx2g_sr, header=FALSE, sep="\t", stringsAsFactors=FALSE)
sr_map <- setNames(strip(t2g$V2), strip(t2g$V1))
b <- read.table(bambu_tx, header=TRUE, sep="\t", comment.char="",
                colClasses=c("character","character", rep("NULL",11)))
lr_map <- setNames(strip(b$GENEID), strip(b$TXNAME))

## ---- protein-coding gene set (gene has >=1 pc transcript) ----
meta <- readRDS(meta_f); meta$isoform_id <- strip(as.character(meta$isoform_id))
pc_enst <- meta$isoform_id[meta$transcript_biotype == "protein_coding" & !is.na(meta$transcript_biotype)]
pc_genes <- unique(sr_map[pc_enst]); pc_genes <- pc_genes[!is.na(pc_genes)]
cat("protein-coding genes (>=1 pc transcript):", length(pc_genes), "\n")

## ---- detected isoforms -> per-gene counts ----
sr <- as.matrix(read.csv(sr_file, row.names=1, check.names=FALSE)); rownames(sr) <- strip(rownames(sr))
sr_det <- rownames(sr)[rowSums(sr >= 1) >= 1]
sr_g <- sr_map[sr_det]; sr_g <- sr_g[!is.na(sr_g)]

lr <- as.matrix(readRDS(lr_file)); rownames(lr) <- strip(rownames(lr))
lr_det <- rownames(lr)[rowSums(lr >= 1) >= 1]
lr_g <- lr_map[lr_det]; lr_g <- lr_g[!is.na(lr_g)]
lr_novel_det <- sum(grepl("^Bambu", lr_det, ignore.case=TRUE))

per_gene <- function(genes_vec, keep=NULL){
  if(!is.null(keep)) genes_vec <- genes_vec[genes_vec %in% keep]
  as.integer(table(genes_vec))
}
summarize <- function(v,label){
  q<-quantile(v,c(.75,.9,.95))
  data.frame(platform=label,n_genes=length(v),min=min(v),median=median(v),
             mean=round(mean(v),3),q75=unname(q[1]),q90=unname(q[2]),q95=unname(q[3]),
             max=max(v),pct_1iso=round(100*mean(v==1),1),
             n_ge2=sum(v>=2),n_ge5=sum(v>=5),n_ge10=sum(v>=10))
}
bin <- function(v){ as.integer(table(factor(pmin(v,20),levels=1:20))) }

## ---- (A) all genes ----
sr_all <- per_gene(sr_g); lr_all <- per_gene(lr_g)
summ_all <- rbind(summarize(sr_all,"short_read"), summarize(lr_all,"long_read"))
write.csv(summ_all, file.path(out_dir,"isoforms_per_gene_SR_LR_summary.csv"), row.names=FALSE)
dist_all <- data.frame(n_isoforms=c(as.character(1:19),"20+"), short_read=bin(sr_all), long_read=bin(lr_all))
write.csv(dist_all, file.path(out_dir,"isoforms_per_gene_SR_LR_distribution.csv"), row.names=FALSE)

## ---- (B) protein-coding genes only ----
sr_pc <- per_gene(sr_g, pc_genes); lr_pc <- per_gene(lr_g, pc_genes)
summ_pc <- rbind(summarize(sr_pc,"short_read"), summarize(lr_pc,"long_read"))
write.csv(summ_pc, file.path(out_dir,"isoforms_per_gene_SR_LR_summary_protein_coding.csv"), row.names=FALSE)
dist_pc <- data.frame(n_isoforms=c(as.character(1:19),"20+"), short_read=bin(sr_pc), long_read=bin(lr_pc))
write.csv(dist_pc, file.path(out_dir,"isoforms_per_gene_SR_LR_distribution_protein_coding.csv"), row.names=FALSE)

cat("\n== ALL GENES ==\n"); print(summ_all)
cat("\n== PROTEIN-CODING GENES ==\n"); print(summ_pc)
cat("\nLR novel (Bambu) isoforms among detected:", lr_novel_det, "\n")

## ---- plots ----
library(ggplot2)
mkplot <- function(sr_v,lr_v,title,fn){
  df <- rbind(data.frame(platform="short read",n=pmin(sr_v,20)),
              data.frame(platform="long read", n=pmin(lr_v,20)))
  p <- ggplot(df,aes(x=n,fill=platform)) +
    geom_histogram(position="identity",alpha=0.55,binwidth=1,colour="grey30") +
    scale_fill_manual(values=c("short read"="#D7301F","long read"="#2166AC")) +
    scale_x_continuous(breaks=c(1,5,10,15,20),labels=c("1","5","10","15","20+")) +
    labs(x="Detected isoforms per gene (capped at 20+)",y="Number of genes",title=title,fill=NULL) +
    theme_bw(base_size=13) + theme(legend.position=c(0.85,0.85))
  ggsave(file.path(out_dir,fn),p,width=8,height=5.5,dpi=200)
}
mkplot(sr_all,lr_all,"Detected isoforms per gene: SR (salmon) vs LR (bambu) — all genes","isoforms_per_gene_SR_LR.png")
mkplot(sr_pc, lr_pc, "Detected isoforms per gene: SR vs LR — protein-coding genes","isoforms_per_gene_SR_LR_protein_coding.png")
cat("\nDONE\n")
