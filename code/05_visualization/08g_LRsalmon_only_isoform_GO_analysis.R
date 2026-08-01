#!/usr/bin/env Rscript
# 08g_LRsalmon_only_isoform_GO_analysis.R
# QUANTIFIER-CONTROLLED, FILTER-BASED "LR-only" genes: isoforms present in the
# LR-salmon FILTERED matrix (protein-coding + filterByExpr min.count=2) but NOT in
# the SR-salmon FILTERED matrix. Both platforms quantified with salmon vs the same
# release-46 transcriptome and put through the identical 01c filter, so this is a
# real platform difference over expression-filtered protein-coding isoforms, not a
# quantifier or detection-threshold artifact. GO:BP ORA, lab convention.
# Universe = LR-salmon filtered genes. Run under conda_R/4.4.x.
.libPaths("/users/sparthib/R/4.4.x")
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)})
strip <- function(x) sub("\\.[0-9]+$","",x)
dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
sal <- file.path(dd,"10_short_reads/salmon")
fbcb<- file.path(sal,"filtered_by_counts_and_biotype")
out <- file.path(sal,"go_LRsalmon_only"); dir.create(out, showWarnings=FALSE)

## maps
t2g <- read.table(file.path(sal,"tx2gene_release46.tsv"),header=FALSE,sep="\t",stringsAsFactors=FALSE)
iso2gene <- setNames(strip(t2g$V2), strip(t2g$V1))
# gene symbol map from release-46 FASTA header (field6), same source as 08e
hdr <- system(paste0("grep '^>' ",
       "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_all_transcripts.fa",
       " | awk -F'|' '{print $2\"\t\"$6}' | sort -u"), intern=TRUE)
hm <- do.call(rbind, strsplit(hdr,"\t")); ensg2sym <- setNames(hm[,2], strip(hm[,1]))

## filtered matrices -> LR-only isoforms
srf <- readRDS(file.path(fbcb,"filtered_isoform_counts.RDS"))
lrf <- readRDS(file.path(fbcb,"LRsalmon_filtered_isoform_counts.RDS"))
sr_iso <- strip(rownames(srf)); lr_iso <- strip(rownames(lrf))
lronly <- setdiff(lr_iso, sr_iso)
genes  <- unique(iso2gene[lronly]); genes <- genes[!is.na(genes)]
universe <- unique(iso2gene[lr_iso]); universe <- universe[!is.na(universe)]
cat("LR-only isoforms:", length(lronly), " genes:", length(genes), " universe:", length(universe), "\n")

## gene list with symbol + n LR-only isoforms
pgc <- table(iso2gene[lronly])
gl <- data.frame(gene_id=names(pgc), gene_symbol=ensg2sym[strip(names(pgc))],
                 n_LRonly_isoforms=as.integer(pgc))
gl <- gl[order(-gl$n_LRonly_isoforms),]
write.csv(gl, file.path(out,"LRsalmonOnly_genes.csv"), row.names=FALSE)
write.csv(data.frame(isoform_id=lronly, gene_id=iso2gene[lronly],
          gene_symbol=ensg2sym[strip(iso2gene[lronly])]),
          file.path(out,"LRsalmonOnly_isoforms.csv"), row.names=FALSE)

## GO:BP ORA (lab convention)
eg <- enrichGO(gene=genes, universe=universe, OrgDb=org.Hs.eg.db, keyType="ENSEMBL",
               ont="BP", pAdjustMethod="fdr", minGSSize=100, pvalueCutoff=0.01,
               qvalueCutoff=0.01, readable=TRUE)
df <- as.data.frame(eg)
cat("GO:BP terms:", nrow(df), "\n")
if(nrow(df)>0){
  write.csv(df, file.path(out,"GO_BP_LRsalmonOnly.csv"), row.names=FALSE)
  p <- dotplot(eg, showCategory=15) + ggtitle("GO:BP — LR-salmon-only isoform genes (filtered)")
  ggsave(file.path(out,"GO_BP_LRsalmonOnly.png"), p, width=9, height=7, dpi=200)
  cat("top:", paste(head(df$Description,5),collapse=" | "), "\n")
}
cat("DONE\n")
