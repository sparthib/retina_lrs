# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# install.packages("markdown")
# BiocManager::install('org.Hs.eg.db')
# 
# 
# BiocManager::install("pcaExplorer")
# BiocManager::install("airway")

library("markdown")
library("pcaExplorer")
library("dplyr")
library("airway")
library("readr")
library("here")
library("ggplot2")
library("tidyr")
library("edgeR")
library("ggfortify")




#### Create Sample Metadata #####

group <- c("RGC","RO_D209" , "RO_D45",
        "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC")
sample_name <- c("DG-WT-hRGC" ,  "EP1-BRN3B-RO" , "EP1-WT_ROs_D45", 
                 "H9-BRN3B-RO", "H9-CRX_ROs_D45", "H9-hRGC_1" , "H9-hRGC_2" ,  
                 "hRGC")

counts <- read_tsv(here("processed_data",
                          "dtu",
                          "IsoformSwitchAnalyzeR",
                          "salmon_alignment_mode_high_mapq",
                          "extracted_gene_counts.tsv"))

annotation <- counts[,1:2]
write_tsv(annotation, here("processed_data",
                       "dtu",
                       "IsoformSwitchAnalyzeR",
                       "salmon_alignment_mode_high_mapq",
                       "gene_annotation.tsv"))

counts <- counts |> dplyr::select(-c(gene_name))
counts <- counts |>
  mutate_if(is.numeric, as.integer)


logcounts <- log10(as.matrix(counts[,2:9] + 1))
rownames(logcounts) <- counts$gene_id

#### convert to cpm #####
mat <- as.matrix(counts[2:9])
rownames(mat) <- counts$gene_id
d <- DGEList(mat)
d <- calcNormFactors(d)
d$samples$group <- group

d$logcpm <- cpm(d, prior.count=2, log=TRUE)

d_new <- d[-drop,] 
dim(d_new) # number of genes left



pr <- prcomp(t(d$logcpm))
summary(pr)
#variance explained by PCA
pc_eigenvalues <- pr$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) |>
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) |>
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

# PC    variance      pct pct_cum
# <fct>    <dbl>    <dbl>   <dbl>
#   1 1     4.27e+ 2 5.01e+ 1    50.1
# 2 2     2.25e+ 2 2.64e+ 1    76.4
# 3 3     9.27e+ 1 1.09e+ 1    87.3
# 4 4     5.14e+ 1 6.02e+ 0    93.3
# 5 5     2.45e+ 1 2.87e+ 0    96.2
# 6 6     1.86e+ 1 2.18e+ 0    98.3
# 7 7     1.41e+ 1 1.65e+ 0   100  
# 8 8     8.78e-28 1.03e-28   100  

pc_eigenvalues |> 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# 8 components are enough to virtually explain all of the variance in our dataset. 
# This makes sense in this case since we only have 8 biological samples.


#### Visualising samples on PC space ####

pc_scores <- pr$x


pc_scores <- pc_scores |> 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores
# sample           PC1    PC2     PC3    PC4    PC5     PC6     PC7      PC8
# <chr>          <dbl>  <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl>    <dbl>
#   1 DG-WT-hRGC     -20.2  -3.83  18.1     8.62 -0.955 -0.302  -1.06   2.62e-14
# 2 EP1-BRN3B-RO    18.5  16.9    1.84    1.79  7.95   5.05   -0.224  3.00e-14
# 3 EP1-WT_ROs_D45  19.0 -18.6   -2.40    1.48  4.93  -7.01   -0.0637 2.33e-14
# 4 H9-BRN3B-RO     16.9  23.8    1.40   -1.65 -6.52  -4.61   -0.0122 2.66e-14
# 5 H9-CRX_ROs_D45  22.1 -18.1   -2.53    1.26 -6.29   5.82    0.0799 2.40e-14
# 6 H9-hRGC_1      -24.4   4.85 -17.3     8.42 -0.407 -0.0128 -0.165  2.71e-14
# 7 H9-hRGC_2      -16.1  -2.85  -0.722 -11.8   0.608  0.600  -6.22   2.94e-14
# 8 hRGC           -15.9  -2.18   1.62   -8.14  0.680  0.465   7.67   2.83e-14


pc_scores$group <-  c("RGC","RO_D209" , "RO_D45",
                              "RO_D209" , "RO_D45" ,"RGC","RGC", "RGC") 
  # create the plot
pc_scores |> ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(group)))

# #Exploring correlation between genes and PCs
pc_loadings <- pr$rotation
pc_loadings <- pc_loadings |>
  as_tibble(rownames = "gene")
pc_loadings
# 
# # gene               PC1      PC2      PC3      PC4     PC5      PC6      PC7
# # <chr>            <dbl>    <dbl>    <dbl>    <dbl>   <dbl>    <dbl>    <dbl>
# #   1 ENSG00000002… -0.0330  -0.00661 -0.0433   0.00910  0.0115  0.0494  -0.0168 
# # 2 ENSG00000003…  0.0579   0.00420  0.00119 -0.00172 -0.0543 -0.0774   0.0754 
# # 3 ENSG00000003… -0.00276 -0.0120   0.00776 -0.00952  0.0203  0.00505 -0.00446
# # 4 ENSG00000004…  0.0441   0.0151  -0.0594   0.0350  -0.0525  0.137   -0.00699
# # 5 ENSG00000004…  0.0104   0.0229  -0.0223   0.0180  -0.0360  0.0217  -0.0180 
# # 6 ENSG00000004… -0.0598  -0.0290   0.00367 -0.0222  -0.0189 -0.0107  -0.0217 
# # 7 ENSG00000005… -0.0645   0.0262   0.0317   0.0225  -0.0359  0.110    0.0165 
# # 8 ENSG00000005…  0.0250  -0.0260  -0.0299  -0.0393  -0.0148 -0.0655   0.00769
# # 9 ENSG00000005… -0.0258  -0.00374  0.00804  0.00331 -0.0558  0.00110 -0.0200 
# # 10 ENSG00000005…  0.0274   0.00600 -0.0424  -0.0211  -0.0350 -0.0366   0.00129

# “What are the top 10 genes with highest loading on PC1 and PC2?”

top_genes <- pc_loadings |>
  # select only the PCs we are interested in
  select(gene, PC1, PC2) |>
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) |>
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) |>
  # pull the gene column as a vector
  pull(gene) |>
  # ensure only unique genes are retained
  unique()

top_genes <- annotation |> filter( gene_id %in% top_genes ) |> select(gene_name)
# 1 CNGB1    
# 2 FXYD3    
# 3 RP1      
# 4 SLC17A7  
# 5 CHN2     
# 6 AJUBA    
# 7 COL2A1   
# 8 TAGLN    
# 9 CCNB2    
# 10 CRYAA    
# 11 SUCLG2   
# 12 GPR160   
# 13 CABP4    
# 14 CNTN2    
# 15 TPM2     
# 16 L1CAM    
# 17 CRYGS 

# [1] "ENST00000380518" "ENST00000307944" "ENST00000416167" "ENST00000355897" "ENST00000288207"
# [6] "ENST00000329305" "ENST00000262713" "ENST00000220676" "ENST00000307227" "ENST00000242576"
# [11] "ENST00000311183" "ENST00000633708" "ENST00000436393" "ENST00000354190" "ENST00000528282"
# [16] "ENST00000423266" "ENST00000706164" "ENST00000608049" "ENST00000251102"


output_plots_dir <- "/Users/sparthib/Documents/retina_lrs/plots/de/switch_analyzer/alignment_mode_high_mapq/"

p <- autoplot(pr, data = pc_scores, color = "group") +
   geom_text(aes(label = sample, vjust=1.2), size = 1.5) +
  ggtitle("PCA on logCPM of genes")
pdf(paste0(output_plots_dir,"gene_logCPM_pca.pdf"), 
    width=8, height=4)
print(p)
dev.off()





