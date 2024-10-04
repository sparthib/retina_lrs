library(scuttle)
library(data.table)
library(Matrix)
library(rtracklayer)
library(scater)
library(sessioninfo)
library(bluster)
library(scran)


sce_list <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/transcriptome/01_quality_controlled/sce_list.rds")



##reduced dimensions 
sce_list <- lapply(sce_list, function(x) {
    x <- scater::runPCA(x)
    x <- scater::runTSNE(x, perplexity = 0.1)
    # Perplexity should be lower than K!
    u <- uwot::umap(t(logcounts(x)), n_neighbors = 2)
    reducedDim(x, "UMAP_uwot") <- u
    x # Now stored in the object.
})


set.seed(100)
##size factors and normalization
sce_list <- lapply(sce_list, function(x) {
  cluster_x <- quickCluster(x) 
  x <- computeSumFactors(x, cluster=cluster_x, min.mean=0.1)
  x <- logNormCounts(x)
  x
})



#### sce_list

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
