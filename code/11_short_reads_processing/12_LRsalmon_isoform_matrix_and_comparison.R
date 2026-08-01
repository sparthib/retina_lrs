#!/usr/bin/env Rscript
# 12_LRsalmon_isoform_matrix_and_comparison.R
# Quantifier-controlled isoform comparison. The 7 retinal-organoid LONG reads were
# quantified with salmon alignment-mode (script 11) against the SAME release-46
# transcriptome as the short reads. This script:
#   (1) builds the LR-salmon PC + expression-filtered isoform matrix using the SAME
#       01c recipe as bambu/SR (protein_coding isoforms, filterByExpr min.count=2, TMM);
#   (2) compares isoforms-per-gene across three FILTERED matrices, holding the
#       filter constant so the only variable is quantifier x platform:
#         SR-salmon (28s) | LR-salmon (7s) | LR-bambu (7s).
# Detection rule: >=1 count in >=1 sample within the filtered matrix.
# Run: module load conda_R/4.5.x && Rscript 12_LRsalmon_isoform_matrix_and_comparison.R

.libPaths("/users/sparthib/R/4.5.x")
suppressMessages(library(edgeR))
strip <- function(x) sub("\\.[0-9]+$", "", x)
dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
ref <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
sal <- file.path(dd, "10_short_reads/salmon")
out_dir <- file.path(sal, "filtered_by_counts_and_biotype")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- protein-coding genes + isoform->gene map ----
gl <- system(paste0("awk -F'\t' '$3==\"gene\"' ", file.path(ref,"release_46_primary_assembly.gtf"),
      " | sed -E 's/.*gene_id \"([^\"]+)\".*gene_type \"([^\"]+)\".*/\\1\t\\2/'"), intern=TRUE)
gm <- do.call(rbind, strsplit(gl,"\t")); gbt <- setNames(gm[,2], strip(gm[,1]))
pc_genes <- names(gbt)[gbt=="protein_coding"]
t2g <- read.table(file.path(sal,"tx2gene_release46.tsv"),header=FALSE,sep="\t",stringsAsFactors=FALSE)
iso2gene <- setNames(strip(t2g$V2), strip(t2g$V1))
pc_iso <- names(iso2gene)[iso2gene %in% pc_genes]

## ---- (1) build LR-salmon filtered isoform matrix (same recipe as 01c/SR) ----
lrs <- as.matrix(read.csv(file.path(sal,"LRsalmon_transcript_counts_matrix.csv"),row.names=1,check.names=FALSE))
rownames(lrs) <- strip(rownames(lrs))
lrs <- round(lrs)                                    # salmon estimated counts -> integers for edgeR
# stage group matches the bambu RO design (column order = D45,D45,D100,D100,D100,D200,D200)
lr_group <- factor(c("Stage_1","Stage_1","Stage_2","Stage_2","Stage_2","Stage_3","Stage_3"))
stopifnot(ncol(lrs)==7)
m <- lrs[rownames(lrs) %in% pc_iso, , drop=FALSE]
dge <- DGEList(counts=m); keep <- filterByExpr(dge, group=lr_group, min.count=2)
dge <- dge[keep,]; dge <- calcNormFactors(dge)
lrs_f <- m[keep,,drop=FALSE]
saveRDS(lrs_f, file.path(out_dir,"LRsalmon_filtered_isoform_counts.RDS"))
saveRDS(cpm(dge,normalized.lib.sizes=TRUE), file.path(out_dir,"LRsalmon_filtered_isoform_cpm.RDS"))
cat("LR-salmon PC isoform subset:", nrow(m), "-> filtered:", nrow(lrs_f), "\n")

## ---- load the other two FILTERED matrices ----
srs_f <- readRDS(file.path(out_dir,"filtered_isoform_counts.RDS"))                       # SR-salmon 117762
lrb_f <- readRDS(file.path(dd,"06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype/filtered_isoform_counts.RDS")) # LR-bambu 55708

## ---- isoforms per gene (detected in filtered matrix) ----
b <- read.table(file.path(dd,"06_quantification/bambu/all_samples_extended_annotation_track_reads/counts_transcript.txt"),
                header=TRUE,sep="\t",comment.char="",colClasses=c("character","character",rep("NULL",11)))
bmap <- setNames(strip(b$GENEID), strip(b$TXNAME))
perg <- function(mat, mp){ det<-rownames(mat)[rowSums(mat>=1)>=1]; g<-mp[strip(det)]; g<-g[!is.na(g)]; as.integer(table(g)) }
srs_pg <- perg(srs_f, iso2gene); lrs_pg <- perg(lrs_f, iso2gene); lrb_pg <- perg(lrb_f, bmap)

summ <- function(v,l) data.frame(dataset=l, n_genes=length(v), n_iso=sum(v),
  median=median(v), mean=round(mean(v),3), q90=quantile(v,.9), max=max(v),
  pct_1=round(100*mean(v==1),1), n_ge5=sum(v>=5), n_ge10=sum(v>=10))
out <- rbind(summ(srs_pg,"SR-salmon (28s)"), summ(lrs_pg,"LR-salmon (7s)"), summ(lrb_pg,"LR-bambu (7s)"))
print(out, row.names=FALSE)
write.csv(out, file.path(sal,"isoforms_per_gene_salmon_vs_bambu_FILTERED_summary.csv"), row.names=FALSE)
bin <- function(v) as.integer(table(factor(pmin(v,20),levels=1:20)))
dist <- data.frame(n_iso=c(as.character(1:19),"20+"), SR_salmon=bin(srs_pg), LR_salmon=bin(lrs_pg), LR_bambu=bin(lrb_pg))
write.csv(dist, file.path(sal,"isoforms_per_gene_salmon_vs_bambu_FILTERED_distribution.csv"), row.names=FALSE)
cat("\nDISTRIBUTION\n"); print(dist, row.names=FALSE)
cat("\nDONE\n")
