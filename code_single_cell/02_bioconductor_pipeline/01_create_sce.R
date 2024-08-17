library(scuttle)
library(data.table)
library(Matrix)
library(rtracklayer)


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

output_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/"
plots_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/01_umaps"
# dir.create(output_dir, showWarnings = FALSE)
create_sce <- function(sample){
  # input counts path 
  counts_path <- paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/",
                sample, "/transcript_count.csv.gz")
  if (file.exists(counts_path)){
    mat <- fread(counts_path)
    
    rownames <- paste(mat[[1]], mat[[2]], sep = "_")
    transcript_ids <- mat[[1]]
    gene_ids <- mat[[2]]
    
    # Convert the remaining columns to a matrix
    mat <- as.matrix(mat[, -c(1, 2), with = FALSE])
    rownames(mat) <- rownames
    sparse.mat <- Matrix(mat, sparse = TRUE)
    rm(mat)
    
    sce <- SingleCellExperiment(
      assays = list(counts = sparse.mat)
    )
    
    sce <- scuttle::logNormCounts(sce)
    
    #add sample name to colData
    sce$sample_name <- sample
    
    #add QC cols to sce
    sce <- scuttle::addPerCellQC(sce)
    
    ## on rowData
    rowData(sce)$id <- rownames
    rowData(sce)$transcript_id <- transcript_ids
    rowData(sce)$gene_id <- gene_ids
    sce <- scuttle::addPerFeatureQC(sce)
    rowData(sce)
    
    
    ## dimension reduction 
    sce <- scater::runPCA(sce)
    dim(reducedDim(sce, "PCA"))
    
    sce <- scater::runTSNE(sce, perplexity = 0.1)
    
    sce <- runUMAP(sce, dimred = 'PCA')
    p <- plotUMAP(sce, colour_by="label")
    
    pdf(file.path(plots_dir, paste0(sample, "_umap.pdf")))
    #save sce object
    saveRDS(sce, file.path(plots_dir, paste0(sample,'_sce.rds')))
    dev.off()
    
  } else {
    message(paste0(sample, " counts file does not exist"))
  }
  
}

lapply(sample_names, create_sce)


## add genomics ranges to sce 


