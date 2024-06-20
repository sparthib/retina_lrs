library("ggplot2")
library("DGEobj.utils")


## PCA with transcript level tpm or cpm

#bambu

FT_vs_RGC_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/CPM_transcript.txt",
                                  sep = "\t", row.names = 1, header = TRUE)
ROs_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/CPM_transcript.txt",
                                   sep = "\t", row.names = 1, header = TRUE)
FT_vs_RGC_isoquant_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/OUT.transcript_model_grouped_tpm.tsv",
                                    sep = "\t", row.names = 1)
ROs_isoquant_tpm <-read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/OUT.transcript_model_grouped_tpm.tsv",
                                    sep = "\t", row.names = 1)
FT_vs_RGC_bambu_tpm <- FT_vs_RGC_bambu_tpm[,2:5]    
ROs_bambu_tpm <- ROs_bambu_tpm[,2:8]

#remove rows that have the same value across all columns
removeZeroVarRows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}

FT_vs_RGC_bambu_tpm <- removeZeroVarRows(FT_vs_RGC_bambu_tpm)
ROs_bambu_tpm <- removeZeroVarRows(ROs_bambu_tpm)
FT_vs_RGC_isoquant_tpm <- removeZeroVarRows(FT_vs_RGC_isoquant_tpm)
ROs_isoquant_tpm <- removeZeroVarRows(ROs_isoquant_tpm)



#samples
FT_vs_RGC_bambu_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_bambu_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
                              "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")
FT_vs_RGC_isoquant_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_isoquant_samples <-c("EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
                                "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2")

#reset column names 
colnames(FT_vs_RGC_bambu_tpm) <-  FT_vs_RGC_bambu_samples
colnames(ROs_bambu_tpm) <- ROs_bambu_samples
colnames(FT_vs_RGC_isoquant_tpm) <- FT_vs_RGC_isoquant_samples
colnames(ROs_isoquant_tpm) <- ROs_isoquant_samples

FT_vs_RGC_bambu_tpm <- convertCounts(as.matrix(FT_vs_RGC_bambu_tpm),
                                     unit = "CPM",normalize = "TMM"
                                     )
ROs_bambu_tpm <- convertCounts(as.matrix(ROs_bambu_tpm),
                                      unit = "CPM",normalize = "TMM"
                                      )
FT_vs_RGC_isoquant_tpm <- convertCounts(as.matrix(FT_vs_RGC_isoquant_tpm),
                                       unit = "CPM",normalize = "TMM"
                                       )
ROs_isoquant_tpm <- convertCounts(as.matrix(ROs_isoquant_tpm),
                                        unit = "CPM",normalize = "TMM"
                                        )
#groups 
FT_vs_RGC_bambu_groups <- c("FT", "FT", "RGC", "RGC")
ROs_bambu_groups <- c("RO_D200", "RO_D100", "RO_D45", 
                             "RO_D100", "RO_D200", "RO_D100", "RO_D45")
FT_vs_RGC_isoquant_groups <- c("FT", "FT", "RGC", "RGC")
ROs_isoquant_groups <- c("RO_D200", "RO_D45", "RO_D100",
                                "RO_D200", "RO_D100", "RO_D45", "RO_D100")

# make pc plot function
output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/pca/"
pc_plot <- function(tpm, samples, groups, name){ 
  pc <- prcomp(t(log2(tpm+1)), scale = TRUE)
  pcr <- data.frame(pc$x[,1:2])
  
  pdf(paste0(output_plots_dir, name, "/isoform_level_pca.pdf"), width = 10, height = 5)
  p <- ggplot(pcr, aes(PC1, PC2, color=groups)) + 
    geom_point(size=4) + 
    theme_bw() +
    ggtitle("PCA on Isoform Expression of all genes") + 
    geom_text(aes(label=rownames(pcr)), hjust=0, vjust=0, size=2) +
    geom_jitter()
  print(p)
  dev.off()
  }

pc_plot(FT_vs_RGC_bambu_tpm, FT_vs_RGC_bambu_samples, FT_vs_RGC_bambu_groups, "bambu")
pc_plot(ROs_bambu_tpm, ROs_bambu_samples, ROs_bambu_groups, "bambu")
pc_plot(FT_vs_RGC_isoquant_tpm, FT_vs_RGC_isoquant_samples, FT_vs_RGC_isoquant_groups, "isoquant")
pc_plot(ROs_isoquant_tpm, ROs_isoquant_samples, ROs_isoquant_groups, "isoquant")



#make gene level pca plot
FT_vs_RGC_bambu_gene_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/counts_gene.txt",
                                      sep = "\t", row.names = 1, header = TRUE)
ROs_bambu_gene_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/counts_gene.txt",
                                       sep = "\t", row.names = 1, header = TRUE)
FT_vs_RGC_isoquant_gene_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/OUT.gene_grouped_tpm.tsv",
                                        sep = "\t", row.names = 1)
ROs_isoquant_gene_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/OUT.gene_grouped_tpm.tsv",
                                         sep = "\t", row.names = 1)

#convert bambu gene counts to cpm
FT_vs_RGC_bambu_gene_tpm <- convertCounts(as.matrix(FT_vs_RGC_bambu_gene_tpm),
                                          unit = "cpm", normalize = "TMM")
ROs_bambu_gene_tpm <- convertCounts(as.matrix(ROs_bambu_gene_tpm),
                                           unit = "cpm", normalize = "TMM")
FT_vs_RGC_isoquant_gene_tpm <- convertCounts(as.matrix(FT_vs_RGC_isoquant_gene_tpm),
                                             unit = "cpm", normalize = "TMM")
ROs_isoquant_gene_tpm <- convertCounts(as.matrix(ROs_isoquant_gene_tpm),
                                              unit = "cpm", normalize = "TMM")



#remove rows that have the same value across all columns
FT_vs_RGC_bambu_gene_tpm <- removeZeroVarRows(FT_vs_RGC_bambu_gene_tpm)
ROs_bambu_gene_tpm <- removeZeroVarRows(ROs_bambu_gene_tpm)
FT_vs_RGC_isoquant_gene_tpm <- removeZeroVarRows(FT_vs_RGC_isoquant_gene_tpm)
ROs_isoquant_gene_tpm <- removeZeroVarRows(ROs_isoquant_gene_tpm)

#samples
FT_vs_RGC_bambu_gene_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_bambu_gene_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
                                   "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")
FT_vs_RGC_isoquant_gene_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_isoquant_gene_samples <- c("EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
                                      "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2")

#reset column names
colnames(FT_vs_RGC_bambu_gene_tpm) <-  FT_vs_RGC_bambu_gene_samples
colnames(ROs_bambu_gene_tpm) <- ROs_bambu_gene_samples
colnames(FT_vs_RGC_isoquant_gene_tpm) <- FT_vs_RGC_isoquant_gene_samples
colnames(ROs_isoquant_gene_tpm) <- ROs_isoquant_gene_samples

#groups
FT_vs_RGC_bambu_gene_groups <- c("FT", "FT", "RGC", "RGC")
ROs_bambu_gene_groups <- c("RO_D200", "RO_D100", "RO_D45", 
                                  "RO_D100", "RO_D200", "RO_D100", "RO_D45")
FT_vs_RGC_isoquant_gene_groups <- c("FT", "FT", "RGC", "RGC")
ROs_isoquant_gene_groups <- c("RO_D200", "RO_D45", "RO_D100",
                                     "RO_D200", "RO_D100", "RO_D45", "RO_D100")

# make pc plot function
pc_plot_gene <- function(tpm, samples, groups, name){ 
  pc <- prcomp(t(log2(tpm+1)), scale = TRUE)
  pcr <- data.frame(pc$x[,1:2])
  
  pdf(paste0(output_plots_dir, name, "/gene_level_pca.pdf"), width = 10, height = 5)
  p <- ggplot(pcr, aes(PC1, PC2, color=groups)) + 
    geom_point(size=4) + 
    theme_bw() +
    ggtitle("PCA on Gene Expression") + 
    geom_text(aes(label=rownames(pcr)), hjust=0, vjust=0, size=2) +
    geom_jitter()
  print(p)
  dev.off()
}

pc_plot_gene(FT_vs_RGC_bambu_gene_tpm, FT_vs_RGC_bambu_gene_samples, FT_vs_RGC_bambu_gene_groups, "bambu")
pc_plot_gene(ROs_bambu_gene_tpm, ROs_bambu_gene_samples, ROs_bambu_gene_groups, "bambu")
pc_plot_gene(FT_vs_RGC_isoquant_gene_tpm, FT_vs_RGC_isoquant_gene_samples, FT_vs_RGC_isoquant_gene_groups, "isoquant")
pc_plot_gene(ROs_isoquant_gene_tpm, ROs_isoquant_gene_samples, ROs_isoquant_gene_groups, "isoquant")
