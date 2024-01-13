#Code adapted from 

#formatting packages
BiocManager::install("Rsubread")
BiocManager::install("tximport")
install.packages("wesanderson")
BiocManager::install("EDASeq")
BiocManager::install("NOISeq")
BiocManager::install("Rvisdiff")
library(Rsubread)
library(limma)
library(ggplot2)
library(readr)
library(wesanderson)
library(cowplot)
library(edgeR)
library(tximport)
library(GenomicFeatures)
library(tidyr)
library('biomaRt')
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)
library(EDASeq)
library(NOISeq)
### OUTPUT FOLDER ####
# OUT=Sys.getenv(OUTPUT_FOLDER)
# 

### Local desktop Isoquant  ### 
# EP1_BRN3B_transcript_counts <- read_tsv("processed_data/01_IsoQuant_output/EP1_BRN3B_RO/OUT.transcript_tpm.tsv")
# H9_BRN3B_transcript_counts <- read_tsv("processed_data/01_IsoQuant_output/H9_BRN3B_RO/OUT.transcript_tpm.tsv")

# 
# # 9-BRN3B vs EP1-BRN3B ROs
EP1_BRN3B_transcript_counts <- read_tsv("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/EP1-BRN3B-RO/OUT/OUT.transcript_tpm.tsv",
                               col_names=TRUE)
H9_BRN3B_transcript_counts <- read_tsv("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/H9-BRN3B-RO/OUT/OUT.transcript_tpm.tsv",
                            col_names=TRUE)

EP1_BRN3B_transcript_counts$EP1_BRN3B_tpm <- EP1_BRN3B_transcript_counts$TPM
H9_BRN3B_transcript_counts$H9_BRN3B_tpm <- H9_BRN3B_transcript_counts$TPM

EP1_BRN3B_transcript_counts<- EP1_BRN3B_transcript_counts |> select(-TPM)
H9_BRN3B_transcript_counts<- H9_BRN3B_transcript_counts |> select(-TPM)



counts <- merge(EP1_BRN3B_transcript_counts, H9_BRN3B_transcript_counts, by = "#feature_id", all = FALSE)

head(counts)
rownames(counts) <- counts$`#feature_id`
# deal with gene names



gtf <- "/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf"
txdb <- makeTxDbFromGFF(gtf)
# # saveRDS(txdb, "/dcs04/hicks/data/sparthib/txdb.RDS")
# txdb <- readRDS("/dcs04/hicks/data/sparthib/txdb.RDS")


txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID") 
tab <- table(txdf$GENEID) 
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

## DE analysis using limmac


x <- DGEList(counts = counts, group=c("Emb", "HI"))
x$genes <- txdf[match(rownames(x), txdf$TXNAME),]

x$samples
keep <- filterByExpr(x, group=c("Emb", "HI"))
x <- x[keep, , keep.lib.sizes=FALSE]
x <- normLibSizes(x)

##GC content offset 

gene_id <- strsplit2(x$genes$GENEID, "\\.")[,1]
gene_id = as.vector(gene_id)


gene_length_GCcontent <- getGeneLengthAndGCContent(gene_id, org="hsa", mode=c("biomart", "org.db"))
gene_length_GCcontent <- as.data.frame(gene_length_GCcontent)


# > nrow(gene_length_GCcontent)
# [1] 6938

##                  dataset              description    version
## 80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14


counts_and_features <- cbind(x$genes, gene_length_GCcontent)
counts_and_features <- cbind(counts_and_features, counts[,2:3])
counts_and_features$GENEID <- gene_id

gene_counts <- counts_and_features |> group_by(GENEID) |>
  summarise(EP1_BRN3B_gene_tpm = sum(EP1_BRN3B_tpm),
            H9_BRN3B_gene_tpm=sum(H9_BRN3B_tpm),
            length=mean(length),
            gc=mean(gc))

rownames(gene_counts) <- gene_counts$GENEID

bcv <- 0.4
x <- DGEList(counts = gene_counts[,2:3], group=1:2)
x$genes <- gene_counts[,1]

et <- exactTest(x, dispersion=bcv^2)

tt <- topTags(et, n = Inf, adjust.method = "BH", sort.by = "none")



hist(x = gene_counts$length, 
    y = gene_counts$EP1_BRN3B_gene_tpm, breaks = 'FD',
     col = 'skyblue', xlab = 'length', ylab = 'EP1_BRN3B_gene_tpm')





mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

G_list <- getBM( attributes= c("ensembl_gene_id","external_gene_name"),
                values=gene_id,mart= mart)





library(RColorBrewer)
col <- brewer.pal(5, "Set1")

x <- calcNormFactors(x)
cpm <- cpm(x, log=TRUE)
plotMDS(cpm, labels = x$samples$group)
x$samples

human.gene <- grep("^ENST", rownames(x))
sequin.gene <- grep("^R", rownames(x))
cpm.human <- cpm(x[human.gene, ], log=TRUE)
cpm.sequin <- cpm(x[sequin.gene, ], log=TRUE)
# par(mfrow=c(1,2))
pdf(file.path(OUT, "mds_human.pdf"), height = 4, width = 8)
par(mar=c(5.1, 5.1, 3.1, 6.1), xpd=TRUE)
plotMDS(cpm.human, pch = rep(c(1, 2), c(6, 9)),
        col = rep(col, rep(3, 5)), main="ONT human",
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
legend("topright", inset=c(-0.16,0), c("000", "025", "050", "075", "100"), 
       pch = c(1, 2, 2, 2, 1),
       col = col[c(1, 5, 4, 3, 2)],
       title = "Group", title.col = "black",
       text.col = col[c(1, 5, 4, 3, 2)], cex = 1.25)
dev.off()
pdf(file.path(OUT, "mds_human_noLeg.pdf"), height = 4, width = 3.5)
plotMDS(cpm.human, pch = rep(c(1, 2), c(6, 9)),
        col = rep(col, rep(3, 5)), main="ONT human",
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
dev.off()
pdf(file.path(OUT, "mds_sequin.pdf"), height = 5, width = 8)
plotMDS(cpm.sequin, pch = rep(c(1, 2), c(6, 9)),
        col = rep(col, rep(3, 5)), main = "ONT sequin",
        cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5)
legend("topleft", c("000", "025", "050", "075", "100"), 
       text.col = col[c(1, 5, 4, 3, 2)], cex = 1.25, bty = "n")
dev.off()


### Fit a separate model for human genes


x.human <- x[human.gene,]
x.human <- calcNormFactors(x.human)
tt.human <- lapply(comparison, function(x){
  x.tmp <- x.human[, x.human$samples$group %in% x]
  design <- model.matrix(~x.tmp$samples$group)
  v.human <- voom(x.tmp, design = design, plot=TRUE)
  fit.human <- lmFit(v.human)
  # fit.human <- contrasts.fit(fit.human, contrasts=contr)
  efit.human <- eBayes(fit.human)
  # dt.human <- decideTests(efit.human)
  topTable(efit.human, coef=ncol(design), n=Inf)
})
for(i in 1:4){
  write.table(tt.human[[i]], file = paste(OUT, "/topTableHuman", names(tt.human)[i], ".tsv", sep=""),
              sep = "\t")
}

Get DE genes


DE.human.limma <- lapply(tt.human, function(x){
  rownames(x)[x$adj.P.Val < 0.05]
})
DE.sequin.limma <- lapply(tt.sequin, function(x){
  rownames(x)[x$adj.P.Val < 0.05]
})
tx.human.limma <- lapply(tt.human, function(x){
  return(x$TXNAME)
})
tx.sequin.limma <- lapply(tt.sequin, function(x){
  return(x$TXNAME)
})


