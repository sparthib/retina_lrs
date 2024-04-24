library(DEXSeq)
library(Rsamtools)
library(IsoformSwitchAnalyzeR)
library(GenomicFeatures)
library(GenomicRanges)
library(here)
library(readr)
library(DRIMSeq)
#Load the data

#read txt file in 

cts <- read.table(here("processed_data", "bambu", 
                  "CPM_transcript.txt"), header = TRUE, sep = "\t")

dim(cts)

# > colnames(cts)
# [1] "TXNAME"              "GENEID"             
# [3] "DG.WT.hRGC"          "EP1.BRN3B.RO_sorted"
# [5] "EP1.WT_ROs_D45"      "H9.BRN3B.RO"        
# [7] "H9.CRX_ROs_D45"      "H9.FT_1"            
# [9] "H9.FT_2"             "H9.hRGC_1"          
# [11] "H9.hRGC_2"           "hRGC"               
# [13] "YZ.15T_hRGC"         "YZ.3T_hRGC"   

colnames(cts) <- c( "feature_id" , "gene_id" ,"DG.WT.hRGC"  ,"EP1.BRN3B.RO", "EP1.WT_ROs_D45",
                     "H9.BRN3B.RO","H9.CRX_ROs_D45" ,"H9.FT_1","H9.FT_2" ,  "H9.hRGC_1" ,         
                     "H9.hRGC_2", "hRGC" , "YZ.15T_hRGC", "YZ.3T_hRGC")


selected_columns <- c(3, 4,5,6,7,10,11,12)
cts <- cts[rowSums(cts[selected_columns]) > 0,]
nrow(cts)

samps <- data.frame( sample_id = colnames(cts[selected_columns]) ,
                     condition = c("RGC", "RO_D209", "RO_D45", "RO_D209",
                                   "RO_D45","hRGC", "hRGC", "hRGC"))

counts <- cts[, c(1,2, selected_columns)]


d <- dmDSdata(counts=counts, samples=samps)

methods(class=class(d))

n <- 12
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_gene_expr=10)

table(table(counts(d)$gene_id))

design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)


library(DEXSeq)
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])


library(stageR)

strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])


stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

head(dex.padj)
