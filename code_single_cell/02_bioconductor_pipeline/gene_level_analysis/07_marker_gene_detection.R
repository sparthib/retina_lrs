library(scran)
library(scater)


# Load the clustered sce object
sce <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/06_clustered_sce_day100.rds")


rowData(sce)$ENSEMBL_ID <- rownames(sce) 
rownames(sce) <- paste0(rowData(sce)$symbol, "_", rowData(sce)$ENSEMBL_ID)


# Marker gene detection

marker.info <- scoreMarkers(sce, 
                            sce$label_hclust_dyn)
 # statistics for cluster 1.
# [1] "self.average"          "other.average"         "self.detected"        
# [4] "other.detected"        "mean.logFC.cohen"      "min.logFC.cohen"      
# [7] "median.logFC.cohen"    "max.logFC.cohen"       "rank.logFC.cohen"     
# [10] "mean.AUC"              "min.AUC"               "median.AUC"           
# [13] "max.AUC"               "rank.AUC"              "mean.logFC.detected"  
# [16] "min.logFC.detected"    "median.logFC.detected" "max.logFC.detected"   
# [19] "rank.logFC.detected"  


# Section 6.3
# The AUC or Cohen’s d is usually the best choice 
# for general purpose marker detection,
# as they are effective regardless
# of the magnitude of the expression values. 
# metric <- "mean.AUC"

metric <- "min.AUC"
plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/06_marker_gene_detection/sce_d100/"
dir.create(plot_dir, showWarnings = FALSE)

cluster_num <- 9

chosen <- marker.info[[as.character(cluster_num)]]

if (is.null(chosen) || nrow(chosen) == 0) {
  stop("No marker info found for cluster: ", cluster_num)
}

ordered <- chosen[order(chosen[[metric]], decreasing=TRUE),]

# if (nrow(ordered) < 10) {
#   warning("Fewer than 10 markers found!")
# }

top_markers <- rownames(ordered)[seq_len(min(10, nrow(ordered)))]

# Check if markers exist in `sce`
missing_features <- setdiff(top_markers, rownames(sce))
if (length(missing_features) > 0) {
  warning("Some markers not found in `sce`: ", paste(missing_features, collapse=", "))
}

file_name <- paste0(plot_dir, "cluster_", cluster_num, "_top_10_markers_", metric, ".pdf")

pdf(file_name)

# Ensure `plotExpression()` actually receives valid input
if (length(top_markers) > 0) {
  plotExpression(sce, features = top_markers, x = "label_hclust_dyn", colour_by = "label_hclust_dyn")
} else {
  warning("No valid markers to plot.")
}

dev.off()
  


  

