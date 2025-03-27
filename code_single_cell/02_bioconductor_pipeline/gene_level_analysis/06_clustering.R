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

sce_day100 <- sce[, sce$day == "D100"]
dim(sce_day100)
sce_day200 <- sce[, sce$day == "D200"]
dim(sce_day200)
# > dim(sce_day100)
# [1] 14054  2332
# > dim(sce_day200)
# [1] 14054  1408

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




##### separated by day
hclust.dyn <- clusterCells(sce_day100, use.dimred="PCA",
                           BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE,
                                                 cut.params=list(minClusterSize=10, deepSplit=1)))
table(hclust.dyn)

# 1   2   3   4   5   6   7   8   9 
# 788 424 224 209 199 194 161  94  39 
sce_day100$label_hclust_dyn <- factor(hclust.dyn)

pdf(paste0(plot_dir, "04b_hclust_dynamic_tsne_day100.pdf"))
p <- plotReducedDim(sce_day100, "TSNE", colour_by="label_hclust_dyn",
                    text_by = "label_hclust_dyn")
print(p)
dev.off()


hclust.dyn <- clusterCells(sce_day200, use.dimred="PCA",
                           BLUSPARAM=HclustParam(method="ward.D2", cut.dynamic=TRUE,
                                                 cut.params=list(minClusterSize=10, deepSplit=1)))
table(hclust.dyn)

# 1   2   3   4   5   6   7 
# 429 267 258 190 165  61  38 

sce_day200$label_hclust_dyn <- factor(hclust.dyn)

pdf(paste0(plot_dir, "04c_hclust_dynamic_tsne_day200.pdf"))
p <- plotReducedDim(sce_day200, "TSNE", colour_by="label_hclust_dyn",
                    text_by = "label_hclust_dyn")
print(p)
dev.off()


saveRDS(sce_day100, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/06_clustered_sce_day100.rds")
saveRDS(sce_day200, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/06_clustered_sce_day200.rds")
