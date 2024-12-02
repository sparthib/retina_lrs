library(readxl)


mmc4 <- readxl::read_xlsx("/users/sparthib/retina_lrs/raw_data/mmc4.xlsx",
                  skip=2)

sce <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/sce_files/01_post_qc_sce.rds")


library(scater)
lib_sce <- librarySizeFactors(sce)
summary(lib_sce)


library(scran)
clust <- quickCluster(sce) 
sce <- scran::computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)


##### model Gene Var #######

dec <- modelGeneVar(sce)

fit <- metadata(dec)
pdf("/users/sparthib/retina_lrs/single_cell_plots/mean_log_expression.pdf")
plot(fit$mean, fit$var, xlab="Mean of log gene expression across samples",
     ylab="Variance of log gene expression across samples")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dev.off()


#get all values of fit$mean > 6 
fit$mean[fit$mean > 6]

# ENSG00000251562 ENSG00000198712 ENSG00000198804 ENSG00000198886 ENSG00000198899 
# 8.851495        6.228218        6.725954        6.414235        6.988745 
# ENSG00000198938 ENSG00000210082 
# 6.713921        6.536870 

#get the variance for all these genes 
fit$var[fit$mean > 6]
# ENSG00000251562 ENSG00000198712 ENSG00000198804 ENSG00000198886 ENSG00000198899 
# 1.241735        3.423114        3.416316        3.028909        3.315936 
# ENSG00000198938 ENSG00000210082 
# 3.482173        3.034801 


dec[order(dec$bio, decreasing=TRUE),] 


# Taking the top 1000 genes here:
hvg <- getTopHVGs(dec, n=1000)
str(hvg)

# Performing PCA only on the chosen HVGs.
sce <- runPCA(sce, subset_row=hvg)
reducedDimNames(sce)
rowSubset(sce, "hvg_1000") <- hvg
colnames(rowData(sce))


##### dimensionality reduction #######
set.seed(100)
sce <- fixedPCA(sce, subset.row=hvg) 
reducedDimNames(sce)

dim(reducedDim(sce, "PCA"))


percent.var <- attr(reducedDim(sce), "percentVar")
# plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")

library(patchwork)
pdf("/users/sparthib/retina_lrs/single_cell_plots/PCA_variance_explained.pdf",
    width=10, height=5)
p1 <- plotReducedDim(sce, dimred="PCA", colour_by="day")
p2 <- plotReducedDim(sce, dimred="PCA", colour_by="discard")
p3 <- plotReducedDim(sce, dimred="PCA", colour_by="sample_base")
p <- p1 | p2 | p3
print(p)
dev.off()


sce <- runTSNE(sce, dimred="PCA")
pdf("/users/sparthib/retina_lrs/single_cell_plots/tSNE.pdf",
    width=10, height=5)
p1 <- plotReducedDim(sce, dimred="TSNE", colour_by="day")
p2 <- plotReducedDim(sce, dimred="TSNE", colour_by="discard")
p3 <- plotReducedDim(sce, dimred="TSNE", colour_by="sample_base")
p <- p1 | p2 | p3
print(p)
dev.off()


##### UMAP #######

set.seed(1100101001)
sce <- runUMAP(sce , dimred="PCA")
pdf("/users/sparthib/retina_lrs/single_cell_plots/UMAP.pdf",
    width=10, height=5)
p1 <- plotReducedDim(sce, dimred="UMAP", colour_by="day")
p2 <- plotReducedDim(sce, dimred="UMAP", colour_by="discard")
p3 <- plotReducedDim(sce, dimred="UMAP", colour_by="sample_base")
p <- p1 | p2 | p3
print(p)
dev.off()


##### Clustering #######
library(scran)
clust.50<- clusterCells(sce, use.dimred="PCA", 
                            BLUSPARAM=SNNGraphParam(k=50, type="rank", cluster.fun="walktrap"))


# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 134 132 147 205 242 156 171 357 106  55 367  77 152  34  45  45  46 811  86  69 
# 21  22  23  24  25  26  27  28  29  30  31 
# 49  67  83  72  24  29  33  74 229  34  54 

colLabels(sce) <- clust.50
pdf("/users/sparthib/retina_lrs/single_cell_plots/PCA_TSNE_20nn_clusters.pdf",
    width=10, height=5)
plotReducedDim(sce, "TSNE", colour_by="label")

dev.off()


#### mmc4

colnames(mmc4)

mmc4_cell_types  <- mmc4 |> dplyr::select("Cones", "RGCs", "Neuorgenic Cells", "MG/RPCs","Amacrine Cells",
                      "Horizontal Cells", "Photoreceptors (Rod/Cones)",
                      "Photoreceptor Precursors/Bipolar Cells/Photoreceptors"
                      )


#get ensembl ids for the genes in the mmc4 data
rgc_genes <- mmc4$RGCs
rgc_genes <- as.character(rgc_genes)
rgc_genes<- rgc_gene_ids[!is.na(rgc_genes)]

# get the gene names for the ensembl ids

library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

rgc_gene_ids <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                        filters="external_gene_name",
                        values=rgc_genes,
                        mart=mart)


#### Marker detection 
marker.info <- scoreMarkers(sce, colLabels(sce))

colnames(marker.info[["1"]])

chosen <- marker.info[["1"]]
ordered <- chosen[order(chosen$rank.logFC.cohen),]

top.ranked <- ordered[ordered$rank.logFC.cohen <= 5,]
rownames(top.ranked)

de_genes <- rowData(sce)[rownames(top.ranked),"gene_symbol"] 
# [1] "JUN"     "TF"      "PTN"     "CLU"     "RDH10"   "VIM"     "WIF1"   
# [8] "TMED10"  "UQCR11"  "MT-ND4"  "TMSB4X"  "BTG2"    "SPP1"    "GPM6A"  
# [15] "TALAM1"  "CRABP1"  "COX5A"   "RCVRN"   "MT-CYB"  "FRZB"    "SFRP2"  
# [22] "FOS"     "ZFP36L1" "RPL4"    "TTC3"    "MT-ND5"  "ITM2B"   "RPL23"  
# [29] "RPL38"   "APOE"    "APP"     "MT-CO1"  "MT-ATP6" "MT-CO3"  "GPM6B"  
# [36] "ENO1"    "HES1"    "FABP7"   "ACTB"    "FIS1"    "PTGDS"   "CD81"   
# [43] "IER2"    "MT-CO2"  "MT-ND3"  "BEX1"   


#convert to entrez IDs
de_entrez_gene_ids <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters="external_gene_name",
                        values=de_genes,
                        mart=mart)




library(limma)
go.out <- goana(unique(de_entrez_gene_ids$entrezgene_id), species="Hs")

# Only keeping biological process terms that are not overly general.
go.out <- go.out[order(go.out$P.DE),]
go.useful <- go.out[go.out$Ont=="BP" & go.out$N <= 200,]
head(go.useful[,c(1,3,4)], 18)



### cluster based on rgc_gene_ids 
rgc_gene_ids <- rgc_gene_ids$ensembl_gene_id
rgc_gene_ids <- as.character(rgc_gene_ids)

rgc_gene_ids <- rgc_gene_ids[!is.na(rgc_gene_ids)]

rgc_gene_ids <- intersect(rownames(sce), rgc_gene_ids)
set.seed(100) 
sce <-  fixedPCA(sce, subset.row = rgc_gene_ids)


percent.var <- attr(reducedDim(sce), "percentVar")
sce <- runTSNE(sce, dimred="PCA")

pdf("/users/sparthib/retina_lrs/single_cell_plots/PCA_TSNE_variance_explained_rgcs.pdf",
    width=10, height=5)
plotReducedDim(sce, dimred="TSNE", colour_by="day")
dev.off()




