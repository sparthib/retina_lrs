library("ggplot2")
library("DGEobj.utils")


## PCA with transcript level tpm or cpm

#bambu

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
counts <- read.table(file.path(isoquant_dir, "OUT.transcript_grouped_counts.tsv"), header = TRUE)
tpm <- read.table(file.path(isoquant_dir, "OUT.transcript_grouped_tpm.tsv"), header = TRUE)

samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                  "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
                  "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")
colnames(counts) <- c("gene_id", sample_names)
colnames(tpm) <- c("gene_id", sample_names)

RO_tpm <- tpm[,1:8]

#remove rows that have the same value across all columns
removeZeroVarRows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}


RO_tpm  <- removeZeroVarRows(RO_tpm)

#samples
# FT_vs_RGC_bambu_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
# ROs_bambu_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
#                               "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")
# FT_vs_RGC_isoquant_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")

ROs_isoquant_samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                                             "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2")


ROs_isoquant_tpm <- convertCounts(as.matrix(ROs_isoquant_tpm),
                                        unit = "CPM",normalize = "TMM"
                                        )

ROs_isoquant_groups <- factor(c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"))

# make pc plot function
output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/pca/"
pc_plot <- function(tpm, samples, groups, name, samples_name) {
  # Perform PCA
  pc <- prcomp(t(log2(tpm + 1)), scale = TRUE)
  pcr <- data.frame(pc$x[, 1:2]) # PC scores of samples on PC1 and PC2 
  
  # Specify the PDF output file and dimensions
  pdf(paste0(output_plots_dir, name, "/", samples_name, "isoform_level_pca.pdf"), width = 12, height = 6)
  
  # Create the PCA plot
  p <- ggplot(pcr, aes(PC1, PC2, color = groups)) +  
    geom_point(size = 2) +  # Reduced point size
    theme_bw() +
    ggtitle("PCA on Isoform Expression") +  
    xlab(paste("PC1 (", round(pc$sdev[1]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +
    ylab(paste("PC2 (", round(pc$sdev[2]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +  
    geom_label(aes(label = rownames(pcr)), size = 2, fill = "white", alpha = 0.7) 
  
  # Variance explained for Scree plot
  var_explained = pc$sdev[1:10]^2 / sum(pc$sdev^2)
  var_explained_df <- data.frame(
    Principal_Component = factor(1:length(var_explained), levels = 1:length(var_explained)),
    Variance_Explained = var_explained
  )
  
  # Create the scree plot
  q <- ggplot(var_explained_df, aes(x = Principal_Component, y = Variance_Explained)) + 
    geom_bar(stat = "identity", fill = "skyblue") +  # Using bars
    xlab("Principal Component") +  
    ylab("Variance Explained") +  
    ggtitle("Scree Plot") +  
    ylim(0, 1) +
    theme_minimal()  # Optional: change the theme for better aesthetics
  
  # Print the combined plot to the PDF device
  plt <- p + q
  print(plt)  
  
  # Close the PDF device
  dev.off()  
}


# pc_plot(FT_vs_RGC_bambu_tpm, FT_vs_RGC_bambu_samples, FT_vs_RGC_bambu_groups, "bambu")
# pc_plot(ROs_bambu_tpm, ROs_bambu_samples, ROs_bambu_groups, "bambu")
# pc_plot(FT_vs_RGC_isoquant_tpm, FT_vs_RGC_isoquant_samples, FT_vs_RGC_isoquant_groups, "isoquant")
pc_plot(ROs_isoquant_tpm, ROs_isoquant_samples, ROs_isoquant_groups, "isoquant")



#make gene level pca plot
isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"
counts <- read.table(file.path(isoquant_dir, "OUT.gene_grouped_counts.tsv"), header = TRUE)
tpm <- read.table(file.path(isoquant_dir, "OUT.gene_grouped_tpm.tsv"), header = TRUE)

samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
             "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
             "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")
colnames(counts) <- c("gene_id", sample_names)
colnames(tpm) <- c("gene_id", sample_names)

RO_tpm <- tpm[,1:8]

#convert bambu gene counts to cpm

ROs_isoquant_gene_tpm <- convertCounts(as.matrix(RO_tpm),
                                              unit = "cpm", normalize = "TMM")


#remove rows that have the same value across all columns
ROs_isoquant_gene_tpm <- removeZeroVarRows(ROs_isoquant_gene_tpm)

#samples
# FT_vs_RGC_bambu_gene_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
# ROs_bambu_gene_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
#                                    "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")
# FT_vs_RGC_isoquant_gene_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_isoquant_gene_samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                               "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2")

#reset column names
# colnames(FT_vs_RGC_bambu_gene_tpm) <-  FT_vs_RGC_bambu_gene_samples
# colnames(ROs_bambu_gene_tpm) <- ROs_bambu_gene_samples
# colnames(FT_vs_RGC_isoquant_gene_tpm) <- FT_vs_RGC_isoquant_gene_samples
colnames(ROs_isoquant_gene_tpm) <- ROs_isoquant_gene_samples

#groups
# FT_vs_RGC_bambu_gene_groups <- c("FT", "FT", "RGC", "RGC")
# ROs_bambu_gene_groups <- c("RO_D200", "RO_D100", "RO_D45", 
#                                   "RO_D100", "RO_D200", "RO_D100", "RO_D45")
# FT_vs_RGC_isoquant_gene_groups <- c("FT", "FT", "RGC", "RGC")
ROs_isoquant_gene_groups <- c("RO_D200", "RO_D45", "RO_D100",
                                     "RO_D200", "RO_D100", "RO_D45", "RO_D100")
ROs_isoquant_groups <- factor(c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"))


# make pc plot function
# Load necessary library


pc_plot_gene <- function(tpm, samples, groups, name){ 
  pc <- prcomp(t(log2(tpm + 1)), scale = TRUE)
  pcr <- data.frame(pc$x[, 1:2])
  
  # Specify the PDF output file and dimensions
  pdf(paste0(output_plots_dir, name, "/gene_level_pca.pdf"), width = 12, height = 6)
  
  # Create the plot
  p <- ggplot(pcr, aes(PC1, PC2, color = groups)) + 
    geom_point(size = 4) + 
    theme_bw() +
    ggtitle("PCA on Gene Expression") + 
    geom_text(aes(label = rownames(pcr)), hjust = 0, vjust = 0, size = 2) +
    geom_jitter()
  
  # Print the plot to the PDF device
  print(p)
  
  # Close the PDF device
  dev.off()
}


# 
# pc_plot_gene(FT_vs_RGC_bambu_gene_tpm, FT_vs_RGC_bambu_gene_samples, FT_vs_RGC_bambu_gene_groups, "bambu")
# pc_plot_gene(ROs_bambu_gene_tpm, ROs_bambu_gene_samples, ROs_bambu_gene_groups, "bambu")
# pc_plot_gene(FT_vs_RGC_isoquant_gene_tpm, FT_vs_RGC_isoquant_gene_samples, FT_vs_RGC_isoquant_gene_groups, "isoquant")
pc_plot_gene(ROs_isoquant_gene_tpm, ROs_isoquant_gene_samples, ROs_isoquant_gene_groups, "isoquant")
