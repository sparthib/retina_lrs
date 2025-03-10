# It is common to select highly variable genes (HVGs) 
# for downstream processing, for example for constructing a PCA.
# Why not just use the variance? What’s the problem?
# We want to account for the mean-variance relationship.
library(scran)

sce <- readRDS(file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/03_normalized_sce.rds")


dec <- modelGeneVar(sce)
fit <- metadata(dec)

# defaults can occasionally yield an overfitted trend when the 
#few high-abundance genes are also highly variable. 
#In such cases, users can reduce the contribution of those
#high-abundance genes by turning off density weights, 
dec.noweight <- modelGeneVar(sce, density.weights=FALSE)
fit.noweight <- metadata(dec.noweight)

plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/03_feature_selection/"
dir.create(plot_dir, showWarnings = FALSE)

pdf(paste0(plot_dir, "01_mean_variance_relationship.pdf"))
plot(fit$mean, fit$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
curve(fit.noweight$trend(x), col="red", add=TRUE, lwd=2)
dev.off()


### Plot CV2 vs mean

dec.CV2 <- modelGeneCV2(sce)
fit.CV2 <- metadata(dec.CV2)
pdf(paste0(plot_dir, "02_mean_cv2_relationship.pdf"))
plot(fit.CV2$mean, fit.CV2$cv2, log="xy")
curve(fit.CV2$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dev.off()

## here we select HVGs based on the largest ratios  
dec.CV2[order(dec.CV2$ratio, decreasing=TRUE),]

hvg.CV2 <- getTopHVGs(dec.CV2, 
                      var.field="ratio", n=1000)
str(hvg.CV2)

# add HVG to rowData 
rowData(sce)$hvg_top1000ratio[rownames(sce) %in% hvg.CV2] <- TRUE

saveRDS(sce, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/04_hvg_sce.rds")
