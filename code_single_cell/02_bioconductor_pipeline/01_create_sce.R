library(scuttle)
library(data.table)
library(Matrix)
library(rtracklayer)
library(scater)
library(sessioninfo)
library(bluster)
library(scran)



sce_10x_D100_EP1_A1 <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10x_D100-EP1_A1/sce.rds")

counts(sce_10x_D100_EP1_A1)[1:5, 1:5]

print(paste0( "Number of cells found: " , dim(counts(sce_10x_D100_EP1_A1))[2]))

sample_names <- c( "10x_D100-EP1_A1",
                   "10x_D100-EP1_A2",
                   "10x_D100-EP1_B1",
                   "10x_D100-EP1_B2",
                   "10x_D200-EP1-1_A1",
                   "10x_D200-EP1-1_A2",
                   "10x_D200-EP1-1_B1",
                   "10x_D200-EP1-1_B2",
                   "10x_D200-EP1-2_A1",
                   "10x_D200-EP1-2_A2",
                   "10x_D200-EP1-2_B1",
                   "10x_D200-EP1-2_B2")

    
sce <- scuttle::logNormCounts(sce)
sce$sample_name <- sample

#add QC cols to sce
sce <- scuttle::addPerCellQC(sce)
sce <- scuttle::addPerFeatureQC(sce)

rowData(sce)    
   
    
    # Feature selection.
    dec <- modelGeneVar(sce)
    hvg <- getTopHVGs(dec, prop=0.1)
    
    ## dimension reduction 
    sce <- scater::runPCA(sce, ncomponents=25, subset_row=hvg)
    dim(reducedDim(sce, "PCA"))
    
    # Clustering.
    colLabels(sce) <- clusterCells(sce, use.dimred='PCA',
                                   BLUSPARAM=NNGraphParam(cluster.fun="louvain"))    
    
    sce <- scater::runTSNE(sce, perplexity = 0.1)
    sce <- runUMAP(sce, dimred = 'PCA')
    #plot umap
    p <- plotUMAP(sce, colour_by="label")
    pdf(file.path(plots_dir, paste0(sample, "_umap.pdf")))
    print(p)
    dev.off()
    
    #save sce object
    saveRDS(sce, file.path(output_dir, paste0(sample,'_sce.rds')))
}

lapply(sample_names, create_sce)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
