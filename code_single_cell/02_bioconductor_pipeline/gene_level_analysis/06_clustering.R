library(scran)
library(scater)
library(bluster)
library(dendextend) # to make a prettier dendrogram
BiocManager::install("dynamicTreeCut")
library(dynamicTreeCut) # to cut the dendrogram)
sce <- readRDS(file = 
                 "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/05_dr_sce.rds")


# higher k yields a more inter-connected graph and broader clusters.
# Users can exploit this by experimenting with different values of k
# to obtain a satisfactory resolution.


nn.clusters <- clusterCells(sce, use.dimred="PCA", 
                             BLUSPARAM=SNNGraphParam(k= 20, type="rank", cluster.fun="walktrap"))
table(nn.clusters)


colLabels(sce) <- nn.clusters

plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/05_clustering/"
dir.create(plot_dir, showWarnings = FALSE)
pdf(paste0(plot_dir, "01_20NN_clusters_day.pdf"))
p <- plotReducedDim(sce, "TSNE", colour_by="label",
                    text_by = "label")
print(p)
dev.off()

#### hierarchical clustering 

hclust <- clusterCells(sce, use.dimred="PCA",
                            BLUSPARAM=HclustParam(method="ward.D2"), full=TRUE)
tree <- hclust$objects$hclust

tree$labels <- seq_along(tree$labels)
dend <- as.dendrogram(tree, hang=0.1)

labels_colors(dend) <- c(
  "D100"="blue",
  "D200"="red"
)[sce$day][order.dendrogram(dend)]


pdf(paste0(plot_dir, "02_hclust_day.pdf"))
plot(dend)
dev.off()

hclust.dyn <- clusterCells(sce, use.dimred="PCA",
                           BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE,
                                  cut.params=list(minClusterSize=10, deepSplit=1)))
table(hclust.dyn)
# 1    2    3    4    5    6    7    8    9 
# 1569  521  481  383  301  173  168   99   45 



color_palette <- c("red", "blue", "green", "purple", "orange", "pink", "cyan", "brown", "black")

# Assign colors based on cluster labels
labels_colors(dend) <- color_palette[as.integer(hclust.dyn)[order.dendrogram(dend)]]


pdf(paste0(plot_dir, "03_hclust_dynamic.pdf"))
plot(dend)
dev.off()

sce$label_hclust_dyn <- factor(hclust.dyn)

pdf(paste0(plot_dir, "04_hclust_dynamic_tsne.pdf"))
p <- plotReducedDim(sce, "TSNE", colour_by="label_hclust_dyn",
                    text_by = "label_hclust_dyn")
print(p)
dev.off()

saveRDS(sce, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/06_clustered_sce.rds")
