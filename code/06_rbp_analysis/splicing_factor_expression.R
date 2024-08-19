library(readr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
# install.packages("DGEobj.utils")
library("DGEobj.utils")
library(grid)


###### ADD ENSEMBL ID TO THE DATA AND SAVE ######
# head(splicing_factors)
# 
# #convert gene symbol to ENSEMBL ID
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ids <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                filters = "hgnc_symbol", values = splicing_factors$`Gene Symbol`,
#                mart = ensembl)
# 
# 
# # Merge the data
# splicing_factors <- merge(splicing_factors, gene_ids, by.x = "Gene Symbol", by.y = "hgnc_symbol")
# 
# # Save the data
# write_csv(splicing_factors, here("raw_data", "GeneCards-Pathway-Splicing.csv"))


#read data 
splicing_factors <- read_csv("/users/sparthib/retina_lrs/raw_data/GeneCards-Pathway-Splicing.csv")
# Load FT vs RGC gene counts matrix 

##### ADD GENE COUNTS #####
FT_vs_RGC_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/counts_gene.txt",
                                  sep = "\t", row.names = 1, header = TRUE)
ROs_bambu_tpm <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/counts_gene.txt",
                            sep = "\t", row.names = 1, header = TRUE)


removeZeroVarRows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}

FT_vs_RGC_bambu_tpm <- convertCounts(as.matrix(FT_vs_RGC_bambu_tpm),
                                     unit = "CPM",normalize = "TMM"
)
ROs_bambu_tpm <- convertCounts(as.matrix(ROs_bambu_tpm),
                               unit = "CPM",normalize = "TMM"
)

FT_vs_RGC_bambu_tpm <- removeZeroVarRows(FT_vs_RGC_bambu_tpm)
ROs_bambu_tpm <- removeZeroVarRows(ROs_bambu_tpm)


#samples
FT_vs_RGC_bambu_samples <- c("H9_FT_1", "H9_FT_2", "H9_RGC_1", "H9_RGC_2")
ROs_bambu_samples <- c("EP1_BRN3B_RO", "EP1_WT_hRO_2", "EP1_WT_ROs_D45",
                       "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")


#remove "_primary_over_30_chr_only_sorted" in column names 


colnames(FT_vs_RGC_bambu_tpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(FT_vs_RGC_bambu_tpm))
#remove version number from rownames
rownames(FT_vs_RGC_bambu_tpm) <- gsub("\\..*", "", rownames(FT_vs_RGC_bambu_tpm))
head(FT_vs_RGC_bambu_tpm)

colnames(ROs_bambu_tpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(ROs_bambu_tpm))
rownames(ROs_bambu_tpm) <- gsub("\\..*", "", rownames(ROs_bambu_tpm))
head(ROs_bambu_tpm)

FT_vs_RGC_bambu_groups <- factor(c("FT", "FT", "RGC", "RGC"))
ROs_bambu_groups <- factor(c("RO_D200", "RO_D100", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45"))

splicing_factors_df <- splicing_factors |> select(ensembl_gene_id, `Gene Symbol`)
colnames(splicing_factors_df) <- c("gene_id", "gene_name")
#remove duplicate rows in splicing_factors_df
splicing_factors_df <- splicing_factors_df[!duplicated(splicing_factors_df$gene_id),]

#only keep the gene_id and the counts found in splicing_factors$ensembl_gene_id 
FT_vs_RGC_bambu_tpm <- FT_vs_RGC_bambu_tpm[rownames(FT_vs_RGC_bambu_tpm) %in% splicing_factors$ensembl_gene_id,]
ROs_bambu_tpm <- ROs_bambu_tpm[rownames(ROs_bambu_tpm) %in% splicing_factors$ensembl_gene_id,]


quant_name <- "bambu"
compare <- "ROs"
tmm <- ROs_bambu_tpm
groups <- ROs_bambu_groups
output_plots_dir <- "/users/sparthib/retina_lrs/plots/splicing_factor_analysis"


dir.create("/users/sparthib/retina_lrs/plots/splicing_factor_analysis", showWarnings = FALSE)
plot_dir <- "/users/sparthib/retina_lrs/plots/splicing_factor_analysis"

#counts number of duplicate rownames
sum(duplicated(rownames(FT_vs_RGC_bambu_tpm)))
sum(duplicated(rownames(ROs_bambu_tpm)))

sum(duplicated(TMM_significant_genes$gene))

plot_heatmap <- function(dge_dir, quant_name, compare, tmm, groups, output_plots_dir) {
  
  tmm <- as.data.frame(tmm)
  tmm$gene_id <- rownames(tmm)
  TMM_significant_genes <- tmm |>
    inner_join(splicing_factors_df, by = "gene_id")
  
  # Create a new column with gene name and isoform name
  
  TMM_significant_genes <- TMM_significant_genes |> 
    mutate(gene = paste(gene_name, gene_id, sep = "_")) |>
    column_to_rownames(var = "gene") |>
    select(-gene_id, -gene_name)
  
  
  # Convert to matrix and scale
  TMM_matrix <- as.matrix(TMM_significant_genes)
  TMM_matrix <- t(scale(t(TMM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Heatmap annotation and 

  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown"),
                                       annotation_name_gp = gpar(fontsize = 0.75)
                            ))
                            
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen"),
                                       annotation_name_gp = gpar(fontsize = 0.75)
                            ))
  }
  
  # Create PDF output
  pdf(file.path(plot_dir, paste0(compare, "_splicing_factors_heatmap.pdf")))
  
  # Generate and draw heatmap
  ht_list <- Heatmap(TMM_matrix, name = paste0("scaled TMM expression of splicing factors"), row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 1),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
}

plot_heatmap(plot_dir, "bambu", "ROs", ROs_bambu_tpm, ROs_bambu_groups, plot_dir)
plot_heatmap(plot_dir, "bambu", "FT_vs_RGC", FT_vs_RGC_bambu_tpm, FT_vs_RGC_bambu_groups, plot_dir)




