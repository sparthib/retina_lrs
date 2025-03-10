library(scran)
library(scater)
sce <- readRDS(file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/04_hvg_sce.rds")


rowData(sce) |> colnames()
# hvg_top1000ratio

# get all top HVGs from rowData

hvg.CV2 <- rownames(sce)[which(rowData(sce)$hvg_top1000ratio %in% TRUE)]

set.seed(100)
sce <- fixedPCA(sce, subset.row=hvg.CV2) 
reducedDimNames(sce)

# dim(reducedDim(sce, "PCA"))
# 3740   50
#Top 50 PCs and 3740 high quality cells 

plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/04_dimensionality_reduction/"
dir.create(plot_dir, showWarnings = FALSE)

pdf(paste0(plot_dir, "01_pca_variance_explained.pdf"))
percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
dev.off()

sum(percent.var[1:25])

### the plot shows first two PCs ###
pdf(paste0(plot_dir, "02a_first_2PCs_day.pdf"))
plotReducedDim(sce, dimred="PCA", colour_by="day")
dev.off()

pdf(paste0(plot_dir, "02b_first_2PCs_sample.pdf"))
plotReducedDim(sce, dimred="PCA", colour_by="sample")
dev.off()

### plots top 4 PCs ###
pdf(paste0(plot_dir, "03a_first_4PCs_day.pdf"))
plotReducedDim(sce, dimred="PCA", ncomponents=4,
               colour_by="day")
dev.off()

pdf(paste0(plot_dir, "03b_first_4PCs_sample.pdf"))
plotReducedDim(sce, dimred="PCA", ncomponents=4,
               colour_by="sample")
dev.off()



set.seed(100)

# runTSNE() stores the t-SNE coordinates in the reducedDims
# for re-use across multiple plotReducedDim() calls.
sce <- runTSNE(sce, dimred="PCA")
pdf(paste0(plot_dir, "04a_tsne_day.pdf"))
plotReducedDim(sce, dimred="TSNE", colour_by="day")
dev.off()

pdf(paste0(plot_dir, "04b_tsne_sample.pdf"))
plotReducedDim(sce, dimred="TSNE", colour_by="sample")
dev.off()

#### run UMAP
sce  <- runUMAP(sce, dimred="PCA")
pdf(paste0(plot_dir, "05a_umap_day.pdf"))
plotReducedDim(sce, dimred="UMAP", colour_by="day")
dev.off()

pdf(paste0(plot_dir, "05b_umap_sample.pdf"))
plotReducedDim(sce, dimred="UMAP", colour_by="sample")
dev.off()

saveRDS(sce, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/05_dr_sce.rds")




