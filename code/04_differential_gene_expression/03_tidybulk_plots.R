# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(here)

# tidyverse-friendly packages
library(plotly)
library(ggrepel)
library(GGally)
library(tidybulk)
library(tidyHeatmap)
library(tidySummarizedExperiment)
colData <- data.frame(group = c("RGC","RO_D209" , "RO_D45",
                                "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC"),
                     sample_name=c( "DG-WT-hRGC" , "EP1-BRN3B-RO", "EP1-WT_ROs_D45",  "H9-BRN3B-RO", 
                                    "H9-CRX_ROs_D45", "H9-hRGC_1", "H9-hRGC_2" , "hRGC"))

geneCountMatrix <- read_tsv(here("processed_data/dtu/IsoformSwitchAnalyzeR/salmon_alignment_mode_high_mapq/extracted_gene_counts.tsv"))
genes  <- geneCountMatrix$gene_id
geneCountMatrix <- geneCountMatrix[,3:10]
geneCountMatrix <- as.matrix(geneCountMatrix)
rownames(geneCountMatrix) <- genes 


se <- SummarizedExperiment(assays=list(counts=geneCountMatrix),
                           colData=colData)

## filtering
se |> filter(`.feature` == "BSG")

# .feature .sample        counts group   sample_name   
# <chr>    <chr>           <dbl> <chr>   <chr>         
#   1 BSG      DG-WT-hRGC      128.  RGC     DG-WT-hRGC    
# 2 BSG      EP1-BRN3B-RO    570.  RO_D209 EP1-BRN3B-RO  
# 3 BSG      EP1-WT_ROs_D45  164.  RO_D45  EP1-WT_ROs_D45
# 4 BSG      H9-BRN3B-RO     559.  RO_D209 H9-BRN3B-RO   
# 5 BSG      H9-CRX_ROs_D45  345.  RO_D45  H9-CRX_ROs_D45
# 6 BSG      H9-hRGC_1        28.0 RGC     H9-hRGC_1     
# 7 BSG      H9-hRGC_2       193.  RGC     H9-hRGC_2     
# 8 BSG      hRGC            265.  RGC     hRGC  

##combine group_by and summarise to calculate the total number
##of rows (transcripts) for each sample.

se |> group_by(`.sample`) |> summarise(n=n())

# .sample            n
# <chr>          <int>
#   1 DG-WT-hRGC      1205
# 2 EP1-BRN3B-RO    1205
# 3 EP1-WT_ROs_D45  1205
# 4 H9-BRN3B-RO     1205
# 5 H9-CRX_ROs_D45  1205
# 6 H9-hRGC_1       1205
# 7 H9-hRGC_2       1205
# 8 hRGC            1205



##### Pre-processing ####

# We filter out lowly expressed genes using tidybulk `keep_abundant`
# or `identify_abundant`. 
# These functions can use the edgeR filterByExpr function described 
# in (Law et al. 2016) to automatically identify the genes with adequate abundance 
# for differential expression testing.

# Filtering counts
counts_scaled <- se |> tidybulk::identify_abundant(factor_of_interest = group) |>  tidybulk::scale_abundance()
# tidybulk says: the sample with largest library size hRGC was chosen as reference for scaling


## plotting scaled abundance density 

pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/alignment_mode_high_mapq/gene_counts_density.pdf")
p1 <- counts_scaled |>
  
  # Reshaping
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  
  # Plotting
  ggplot(aes(x = abundance + 1, color = sample_name)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10()
print(p1)
dev.off()



### Diff expression 
counts_scal_PCA <- counts_scaled |>
  reduce_dimensions(method = "PCA")

de_all <-
  
  counts_scal_PCA |>
  
  # edgeR QLT
  tidybulk::test_differential_abundance(
    ~ group ,
    method = "edgeR_quasi_likelihood",
    prefix = "edgerQLT_"
  ) |> 
  
  # edgeR LRT
  tidybulk::test_differential_abundance(
    ~ group ,
    method = "edgeR_likelihood_ratio",
    prefix = "edgerLR_"
  ) 

# 
# de_all |> select(c(`.feature`,  `counts_scaled`, 
#                  `TMM`, `edgerQLT_PValue`, `edgerQLT_FDR`)) |>
#   unique()
# 
# 
# 
# de_all |> select(c(`.feature`, `.sample`, `counts_scaled`, 
#                    `TMM`, `edgerLR_PValue`, `edgerLR_FDR`)) |>
#   filter(`.feature` == "A2M")

output_data_dir <- here("processed_data",
                        "dtu",
                    "IsoformSwitchAnalyzeR",
                    "salmon_alignment_mode_high_mapq") 



write_tsv(as_tibble(de_all), paste0(output_data_dir, "/tidy_bulk_dge.tsv"))
###### see what data looks like ####
de_all
# A SummarizedExperiment-tibble abstraction: 447,280 × 27
# Features=55910 | Samples=8 | Assays=counts, counts_scaled
# .feature  .sample    counts counts_scaled group sample_name   TMM multiplier
# <chr>     <chr>       <dbl>         <dbl> <chr> <chr>       <dbl>      <dbl>
#   1 A1BG      DG.WT.hRGC 316.          449.   RGC   DG-WT-hRGC   1.02       1.42
# 2 A1BG-AS1  DG.WT.hRGC 118.          167.   RGC   DG-WT-hRGC   1.02       1.42
# 3 A1CF      DG.WT.hRGC   7.00          9.94 RGC   DG-WT-hRGC   1.02       1.42
# 4 A2M       DG.WT.hRGC  94.9         135.   RGC   DG-WT-hRGC   1.02       1.42
# 5 A2M-AS1   DG.WT.hRGC  17.8          25.3  RGC   DG-WT-hRGC   1.02       1.42
# 6 A2ML1     DG.WT.hRGC  37.9          53.8  RGC   DG-WT-hRGC   1.02       1.42
# 7 A2ML1-AS1 DG.WT.hRGC   0             0    RGC   DG-WT-hRGC   1.02       1.42
# 8 A2ML1-AS2 DG.WT.hRGC   1.00          1.42 RGC   DG-WT-hRGC   1.02       1.42
# 9 A2MP1     DG.WT.hRGC  27.6          39.2  RGC   DG-WT-hRGC   1.02       1.42
# 10 A3GALT2   DG.WT.hRGC   2.15          3.05 RGC   DG-WT-hRGC   1.02       1.42
# ℹ 40 more rows
# ℹ 19 more variables: PC1 <dbl>, PC2 <dbl>, .abundant <lgl>,
#   edgerQLT_logFC <dbl>, edgerQLT_logCPM <dbl>, edgerQLT_F <dbl>,
#   edgerQLT_PValue <dbl>, edgerQLT_FDR <dbl>, edgerLR_logFC <dbl>,
#   edgerLR_logCPM <dbl>, edgerLR_LR <dbl>, edgerLR_PValue <dbl>,
#   edgerLR_FDR <dbl>, voom_logFC <dbl>, voom_AveExpr <dbl>, voom_t <dbl>,
#   voom_P.Value <dbl>, voom_adj.P.Val <dbl>, voom_B <dbl>
# ℹ Use `print(n = ...)` to see more rows

####

pdf("/users/sparthib/retina_lrs/plots/eda/tidy_bulk/de_all_plot.pdf") 
de_all_plot <- de_all |> 
  tidybulk::pivot_transcript() |> 
  select(edgerQLT_PValue, edgerLR_PValue, voom_P.Value) |> 
  ggpairs()
print(de_all_plot)
dev.off()

###########################################
counts_de <- counts_scal_PCA |> 
  tidybulk::test_differential_abundance(~ group,
                                        method = "edgeR_quasi_likelihood",
                                        prefix = "edgerQLT_")
#> =====================================
#> tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance 
#> or adjust_abundance have been calculated. Therefore, it is essential to add covariates 
#> such as batch effects (if applicable) in the formula.
#> =====================================
#> tidybulk says: The design column names are "(Intercept), dexuntrt, cellN061011, cellN080611, cellN61311"
#> tidybulk says: to access the raw results (fitted GLM) do `attr(..., "internals")$edgeR_quasi_likelihood`

# volcano plot

topgenes <-
  counts_de |>
  tidybulk::pivot_transcript() |>
  arrange(edgerQLT_FDR) |> select(`.feature`)


counts_de_for_volcano_annotation <- counts_de |> 
  tidybulk::pivot_transcript() |> 
  
  # Subset data
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) |> 
  mutate(symbol = case_when( `.feature` %in% topgenes_symbols$.feature ~ `.feature`, 
                             TRUE ~ "" )
)

####### volcano plot with top DE labels
pdf(paste0(plot_dir, "/annotated_gene_expression_volcano.pdf"))
annotated_volcano <-  counts_de_for_volcano_annotation |>
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, alpha = significant)) +
  geom_text(hjust=0, vjust=0) +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28"))
print(annotated_volcano)
dev.off()


######### Clustering ###########
counts_scaled_top_500 <- counts_scaled |> # extract 500 most variable genes
keep_variable( .abundance = counts_scaled, top = 500) |> as_tibble()

counts_scaled_top_500 |> heatmap( .column = `.sample`,
           .row = `.feature`,
           .value = counts_scaled,
           transform = log1p,
           row_names_gp = grid::gpar(fontsize = 0.8),
           column_names_gp = grid::gpar(fontsize = 5),
           column_title_gp = grid::gpar(fontsize = 7),
           row_title_gp = grid::gpar(fontsize = 7)) |> 
  add_tile(group,
           palette = c("#CC6677","#332288")) |>
  save_pdf(paste0( plot_dir, "/annotated_gene_expression_cluster.pdf")) 
  
  ## custom palette  add_tile(
# hp, 
# palette = c("red", "white", "blue")
# )

#> Getting the 500 most variable genes
#> Warning: The `annotation` argument of `heatmap()` is deprecated as of tidyHeatmap 1.1.0.
#> Please use the new annotation framework instead: heatmap(...) %>% add_tile(...) %>% add_point(...) %>% add_bar() %>% add_line() %>% ...
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
####### Hierarchical Clustering #########

# make heatmap of 200 most DE genes 

heatmap_df <- as_tibble(counts_de)
nrow(heatmap_df)
heatmap_df <- heatmap_df |> select(c(`.feature`,group,`.sample`,
                                     edgerQLT_logFC, 
                       edgerQLT_logCPM,
                       edgerQLT_FDR)) |> unique()

heatmap_df <- heatmap_df |> mutate(
  cpm = exp(edgerQLT_logCPM)
)
colnames(heatmap_df) <- c("gene_id", "group", "sample",
                          "logFC",  "FDR", "cpm")


significant_DGEs <- heatmap_df |> dplyr::filter(FDR < 0.05) |> 
  filter(abs(logFC) > 0.5 )


                                                  
# get genes_ids of 200 DTU isoforms with the smallest q vals. 
top_200_dge_genes <- significant_DGEs |> arrange(FDR) |> 
  dplyr::select(gene_id) |> unique() |> head(n = 200)



### TMM of top 200 isoforms. 
cpm_top_200 <- significant_DGEs |> filter(gene_id %in% top_200_dge_genes$gene_id)

log_cpm_top_200 <- log_cpm_top_200 |> select(gene_id, logCPM, sample) |> 
  pivot_wider(names_from = sample, values_from = logCPM)


#####clustering based on only top isoforms that show DTU. 
library(ComplexHeatmap)
library(circlize)


col_fun = colorRamp2(c(-2.5, 0, 2.5), c("blue","white", "red"))
col_fun(seq(-3, 3))

gene_id <- log_cpm_top_200$gene_id
log_cpm_top_200 <- as.matrix(log_cpm_top_200[2:9])
log_cpm_top_200 <- t(scale(t(log_cpm_top_200), center = T, scale = T))
rownames(log_cpm_top_200) <- gene_id 

gr <- c("RGC","RO_D209" , "RO_D45",
                "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC")
ha = HeatmapAnnotation(type = gr, annotation_name_side = "left",
                       col = list(type = c("RGC" = "seagreen", "RO_D209" = "purple", "RO_D45" = "orange")
                       ))

pdf(here("plots", "de", "switch_analyzer",
         "alignment_mode_high_mapq", "top_200_genes_heatmap.pdf"))
ht_list = Heatmap(log_cpm_top_200, name = "scaled logCPM expression of top 200 DTU isoforms", row_km = 5, 
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
