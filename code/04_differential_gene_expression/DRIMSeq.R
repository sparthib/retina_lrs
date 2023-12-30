library(ggplot2)
library(DRIMSeq)
# library(DEXSeq)
library(limma)
library(edgeR)
# library(satuRn)
library(SummarizedExperiment)
library(stringr)
library(tximport)
library(readr)

samples <- list.files("/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level/ensembl_fa/")

quant <- file.path("/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level/ensembl_fa/", samples, "quant.sf")

txi <- tximport(quant, type="salmon", txOut=TRUE, countsFromAbundance = "no")

tx2gene <-  read_tsv("/dcs04/hicks/data/sparthib/transcript_lengths_sorted.tsv")
tx2gene <- tx2gene |> dplyr::select(tx_id, gene_id)

samples <- data.frame(
  sample_id = samples,
  group = c("RGC", "RO","RO", "RGC"),
  stringsAsFactors = FALSE)

samples$group <- as.factor(samples$group)
counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)


genes <- rownames(counts)
genes <- limma::strsplit2(genes, "|", fixed=TRUE)
rownames(counts) <- genes[,1]
counts <- counts[match(tx2gene$tx_id, rownames(counts)),]
counts <- cbind(tx2gene$gene_id, tx2gene$tx_id, counts)
colnames(counts)[3:ncol(counts)] <- make.names(colnames(counts)[3:ncol(counts)])
colnames(counts) <- c("gene_id", "feature_id", samples$sample_id)


# DRIMSeq analysis for human genes

ddata <- dmDSdata(counts = counts, samples = samples)

table(DRIMSeq::samples(ddata)$group)
  
  
d <- dmFilter(ddata, min_samps_gene_expr = 2, min_samps_feature_expr = 2,
              min_gene_expr = 10, min_feature_expr = 10, run_gene_twice = T)
  
plotData(d)
  
design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))

set.seed(123)
d <- dmPrecision(d, design_full, BPPARAM = BiocParallel::bpparam())
d <- dmFit(d, design = design, verbose=1, BPPARAM = BiocParallel::bpparam())
dmTest(d, coef = colnames(design_full)[ncol(design_full)])



# head(DRIMSeq::results(d))
# head(DRIMSeq::results(d, level = "feature"))

drres.human <- lapply(d.human, DRIMSeq::results)
drres.txp.human <- lapply(d.human, 
                          function(x){DRIMSeq::results(x, level = "feature")})

DTU.gene.human.DRIMSeq <- lapply(drres.human, function(x){
  na.omit(x$gene_id[x$adj_pvalue < 0.05])
})
DTU.tx.human.DRIMSeq <- lapply(drres.txp.human, function(x){
  na.omit(x$feature_id[x$adj_pvalue < 0.05])
})
gene.human.DRIMSeq <- lapply(drres.human, function(x){
  return(x$gene_id)
})
tx.human.DRIMSeq <- lapply(drres.txp.human, function(x){
  return(x$feature_id)
})
```

```{r DRIMSeq2, echo=FALSE, eval=TRUE}
samples <- data.frame(
  sample_id = samples,
  group = rep(c("000", "100", "075", "050", "025"), rep(3, 5)),
  stringsAsFactors = FALSE
)
samples$group <- as.factor(samples$group)
counts <- as.data.frame(txi$counts, stringAsFactors = FALSE)
# deal with gene names
genes <- rownames(counts)
genes <- strsplit2(genes, "|", fixed=TRUE)
rownames(counts) <- genes[,1]
counts <- counts[match(txdf$TXNAME, rownames(counts)),]
counts <- cbind(txdf$GENEID, txdf$TXNAME, counts)
colnames(counts)[3:ncol(counts)] <- make.names(colnames(counts)[3:ncol(counts)])
colnames(counts) <- c("gene_id", "feature_id", samples$sample_id)

comparison <- list(
  c100vs0=c("100", "000"),
  c75vs25=c("075", "025"),
  c50vs25=c("050", "025"),
  c75vs50=c("075", "050")
)

human.gene <- grep("^ENSG", counts$gene_id)
sequin.gene <- grep("^R", counts$gene_id)

d.human <- readRDS("d.human.RDS")
drres.human <- lapply(d.human, DRIMSeq::results)
drres.txp.human <- lapply(d.human, 
                          function(x){DRIMSeq::results(x, level = "feature")})

DTU.gene.human.DRIMSeq <- lapply(drres.human, function(x){
  na.omit(x$gene_id[x$adj_pvalue < 0.05])
})
DTU.tx.human.DRIMSeq <- lapply(drres.txp.human, function(x){
  na.omit(x$feature_id[x$adj_pvalue < 0.05])
})
gene.human.DRIMSeq <- lapply(drres.human, function(x){
  return(x$gene_id)
})
tx.human.DRIMSeq <- lapply(drres.txp.human, function(x){
  return(x$feature_id)
})



###Limma

x <- DGEList(counts = as.matrix(counts(d)[,-c(1:2)]),
                     samples = samples,
                     genes = counts(d)[,c(1,2)])
rownames(x) <- counts(d)$feature_id
x <- calcNormFactors(x)

design <- model.matrix(~group, x$samples)

v <- voom(x, design, plot=T)
fit.human <- lmFit(v)
# efit.human <- eBayes(fit.human)
diffSplice(fit.human, geneid = x$genes$gene_id, exonid=x$genes$feature_id)

