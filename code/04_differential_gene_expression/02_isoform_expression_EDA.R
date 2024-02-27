library("readr")
library("here")
library("dplyr")
library("matrixStats")
library("ggplot2")
library("scater")
TMM <- read_tsv(here("processed_data",
                    "dtu",
                    "IsoformSwitchAnalyzeR",
                    "salmon_alignment_mode_high_mapq",
                    "isoform_abundance.tsv"))

counts <- read_tsv(here("processed_data",
                        "dtu",
                        "IsoformSwitchAnalyzeR",
                        "salmon_alignment_mode_high_mapq",
                        "isoform_counts.tsv"))



DEXSeq_switch_list <- read_tsv(here("processed_data",
                                    "dtu",
                                    "IsoformSwitchAnalyzeR",
                                    "salmon_alignment_mode_high_mapq",
                                    "DEXSeqSwitchList.tsv"))

isoform_gene <- DEXSeq_switch_list |> select(isoform_id, gene_name)

TMM_gene_id <- merge(TMM, isoform_gene, by = "isoform_id", all.x = TRUE)
TMM_gene_id <- unique(TMM_gene_id)
TMM_gene_id$isoform_gene <- paste0(TMM_gene_id$isoform_id, "_",
                                  TMM_gene_id$gene_name)
rownames(TMM_gene_id) <- TMM_gene_id$isoform_gene

# 
# [1] "isoform_id"     "DG-WT-hRGC"     "EP1-BRN3B-RO"   "EP1-WT_ROs_D45"
# [5] "H9-BRN3B-RO"    "H9-CRX_ROs_D45" "H9-hRGC_1"      "H9-hRGC_2"     
# [9] "hRGC"

gr <- c("RGC","RO_D209" , "RO_D45",
        "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC")

# 
# nrow(DEXSeq_switch_list)
# [1] 23133

# 
# nrow(DEXSeq_switch_list |> filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)) 
# [1] 2096


nrow(TMM_gene_id)

# gr is group name of your samples based on columns in counts matrix.
pc<-prcomp(log2(TMM_gene_id[,2:9]+1))
summary(pc)

# Importance of components:
#   PC1     PC2     PC3     PC4     PC5     PC6     PC7  PC8
# Standard deviation     4.1309 1.40349 0.92211 0.54930 0.45167 0.41241 0.37443 0.10458
# Proportion of Variance 0.8239 0.09511 0.04105 0.01457 0.00985 0.00821 0.00677 0.00053
# Cumulative Proportion  0.8239 0.91902 0.96007 0.97464 0.98449 0.99270 0.99947 1.00000


pcr <- data.frame(pc$r[,1:2], Group=gr)
pcr$short_sample_name <- c("DG_WT", "EP1_BRN3B", "EP1_WT", "H9_BRN3B","H9_CRX",
                           "H9_hRGC_1","H9_hRGC_2", "hRGC")

pdf(here("plots", "de", "switch_analyzer",
         "alignment_mode_high_mapq", "PCA_all_isoforms_log_transformed.pdf"))
p <- ggplot(pcr, aes(PC1, PC2, color=Group))+geom_point(size=4)+theme_bw()+
  ggtitle("PCA on Isoform Expression of all genes") + 
  geom_jitter() 
  geom_text(aes(label = short_sample_name), size = 2 , vjust = -1.5 ) +
  
print(p)
dev.off()

#### transposed #####
pc_transpose<-prcomp(t(log2(TMM_gene_id[,2:9]+1)), scale = TRUE)
summary(pc_transpose)

pcr_transpose <- data.frame(pc_transpose$r[,1:2])

pdf(here("plots", "de", "switch_analyzer",
         "alignment_mode_high_mapq", "PCA_transpose_all_isoforms_log_transformed.pdf"))
p <- ggplot(pcr, aes(PC1, PC2, color=Group))+geom_point(size=4)+theme_bw()+
  ggtitle("PCA on Isoform Expression of all genes") + 
  geom_jitter() 
print(p)
dev.off()


colData <- data.frame(group = c("RGC","RO_D209" , "RO_D45",
                                "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC"),
                      sample_name=c( "DG-WT-hRGC" , "EP1-BRN3B-RO", "EP1-WT_ROs_D45",  "H9-BRN3B-RO", 
                                     "H9-CRX_ROs_D45", "H9-hRGC_1", "H9-hRGC_2" , "hRGC"))

se <- SummarizedExperiment(assays=list(counts=TMM_gene_id[,2:9],
                                       logcounts=log2(TMM_gene_id[,2:9]+1)),
                           colData=colData,
                           rowData = rownames(TMM_gene_id[,2:9]))


library(scater) 
library(scran)
library(PCAtools)
gene_var <- modelGeneVar(se)

p <- gene_var |> 
  # convert to tibble for ggplot
  as_tibble() |> 
  # make the plot
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")
p


hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs) 

se <- se |> logNormCounts() |> scater::runPCA(rank = 50) 
se <-  runPCA(se)

s <- scater::calculatePCA(log2(TMM_gene_id[,2:9]+1))

scater::plotReducedDim(se, dimred = "PCA", colour_by = "")

# 
# By default, the prcomp() function does the centering but not the scaling.
# See the ?prcomp help to see how to change this default behaviour.
# In this case, because we are using the transformed data, this is not too much of an issue, 
# but try re-running the PCA with the centered and scaled data to see how it changes it.
# Note that scaling is particularly recommended if your variables are on different scales.
# 

colRanges(as.matrix(log2(TMM[,2:9]+1)))
# [,1]      [,2]
# DG-WT-hRGC        0  9.384963
# EP1-BRN3B-RO      0  8.653219
# EP1-WT_ROs_D45    0  9.992974
# H9-BRN3B-RO       0  9.925430
# H9-CRX_ROs_D45    0  9.669771
# H9-hRGC_1         0  9.588048
# H9-hRGC_2         0 10.047250
# hRGC              0  9.704026

## scaling not necessary when range of the columns are similar. 


#### Hierarchical Clustering #####

# make heatmap of 200 most DE genes 

nrow(DEXSeq_switch_list)

significant_DTUs <- DEXSeq_switch_list |> dplyr::filter(isoform_switch_q_value < 0.05 &
                                      abs(dIF) > 0.1 ) 

DTU_isoforms <- unique(significant_DTUs$isoform_id)

TMM_significant_isoforms <- TMM_gene_id |> dplyr::filter(isoform_id %in% 
                                                   DTU_isoforms)

# > nrow(TMM)
# [1] 7727
# > nrow(TMM_significant_isoforms)
# [1] 1798


# get isoform_ids of 200 DTU isoforms with the smallest q vals. 
top_200_dTU_isoforms <- significant_DTUs |> arrange( isoform_switch_q_value) |> head(n = 250) |> 
  dplyr::select(isoform_id) |> unique() |> head(n = 200)



### TMM of top 200 isoforms. 
TMM_top_200 <- TMM_significant_isoforms |> filter(isoform_id %in% top_200_dTU_isoforms$isoform_id)

TMM_top_200 <- TMM_top_200 |> select(-c(isoform_id, gene_name))


#####clustering based on only top isoforms that show DTU. 
library(ComplexHeatmap)
library(circlize)


col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue","white", "red"))
col_fun(seq(-3, 3))


TMM_matrix <- as.matrix(TMM_top_200[1:8])
TMM_matrix <- t(scale(t(TMM_matrix), center = T, scale = T))
rownames(TMM_matrix) <-rownames(TMM_top_200)

ha = HeatmapAnnotation(type = gr, annotation_name_side = "left",
                       col = list(type = c("RGC" = "seagreen", "RO_D209" = "purple", "RO_D45" = "orange")
                       ))

pdf(here("plots", "de", "switch_analyzer",
         "alignment_mode_high_mapq", "top_200_isoforms_heatmap.pdf"))
ht_list = Heatmap(TMM_matrix, name = "scaled TMM expression of top 176 DTU isoforms", row_km = 5, 
                  col = col_fun,
                  top_annotation = ha, show_row_names = TRUE,
                  show_column_names = TRUE, row_title = "isoforms",
                  row_names_gp = gpar(fontsize = 2.2),
                  column_names_gp = gpar(fontsize = 5),
                  show_row_dend = T, show_column_dend = T) 
draw(ht_list)
dev.off()

#### try heatmap for all significant genes 

TMM_matrix <- as.matrix(TMM_significant_isoforms[2:9])
rownames(TMM_matrix) <-TMM_significant_isoforms$isoform_id
ha = HeatmapAnnotation(type = gr, annotation_name_side = "left",
                       col = list(type = c("RGC" = "seagreen", "RO_D209" = "purple", "RO_D45" = "orange")
                       ))
pdf(here("plots", "de", "switch_analyzer",
         "alignment_mode_high_mapq", "all_significant_DTU_isoforms_heatmap.pdf"))
ht_list = Heatmap(TMM_matrix, name = "TMM expression of all significant DTU isoforms", row_km = 5, 
                  col = col_fun,
                  top_annotation = ha, show_row_names = FALSE,
                  show_column_names = TRUE, row_title = "isoforms",
                  column_names_gp = gpar(fontsize = 5),
                  show_row_dend = T, show_column_dend = T) 
draw(ht_list)
dev.off()
