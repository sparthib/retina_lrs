library(scuttle)
library(scater)
library(ggplot2)
library(patchwork)
library(tidySingleCellExperiment)

sce <- readRDS(file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/01_preqc_sce.rds")

plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/01_qc/"

## remove lowly expressed genes

min_cells <- 10  # Set threshold for minimum cells expressing a gene
sce <- sce[rowSums(assay(sce) > 0) >= min_cells, ]

dim(assay(sce))
# > dim(assay(sce))
# [1] 14054  4228


# Identifying the mitochondrial transcripts in our SingleCellExperiment.
mito_genes <- grep("^MT-", 
                   rowData(sce)$symbol,
                   value = TRUE)
# [1] "MT-ND6"  "MT-CO2"  "MT-CYB"  "MT-ND2"  "MT-ND5"  "MT-CO1"  "MT-ND3" 
# [8] "MT-ND4"  "MT-ND1"  "MT-ATP6" "MT-CO3"  "MT-ND4L" "MT-ATP8"
mito <- which(rowData(sce)$symbol %in% mito_genes)
df <- perCellQCMetrics(sce, subsets=list(Mt=mito))


summary(df$sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 183    5290    8132    9918   12013   81489 


#get 90th percentile of the sum
quantile(df$sum, 0.99)
# 25892.05


summary(df$detected)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 78    1744    2320    2301    2905    6773 

summary(df$subsets_Mt_percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   5.512   8.410  10.579  12.093  97.336 

sce <- addPerCellQCMetrics(sce, 
                           subsets=list(Mito=mito))
#### QC with fized threshold ######

qc.lib <- df$sum < 1e3 | df$sum > 40011
qc.nexprs <- df$detected < 500
qc.mito <- df$subsets_Mt_percent > 20
discard <- qc.lib | qc.nexprs  | qc.mito

sce$discard_fixed <- discard

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(discard))

## adaptive thresholding ##

reasons <- perCellQCFilters(df, 
                            sub.fields=c("subsets_Mt_percent"))
colSums(as.matrix(reasons))
sce$discard_adaptive <- reasons$discard

attr(reasons$low_n_features, "thresholds")
# lower   higher 
# 768.7397      Inf 

attr(reasons$low_lib_size, "thresholds")
# > attr(reasons$low_lib_size, "thresholds")
# lower  higher 
# 1334.98     Inf 

attr(reasons$high_subsets_Mt_percent, "thresholds")

# > summary(reasons$discard)
# Mode   FALSE    TRUE 
# logical    3664     564

#### Diagnostic plot with fixed thresholding ####
pdf(file = paste0(plot_dir, "qc_plots_fixed_threshold.pdf"), width = 10, height = 10)
gridExtra::grid.arrange(
  plotColData(sce, x="day", y="sum", colour_by="discard_fixed") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="day", y="detected", colour_by="discard_fixed") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="day", y="subsets_Mito_percent", 
              colour_by="discard_fixed") + ggtitle("Mito percent"),
  ncol=1
)
dev.off()


### Diagnostic plot with adaptive thresholding ####
pdf(file = paste0(plot_dir, "qc_plots_adaptive_threshold.pdf"), width = 10, height = 10)
gridExtra::grid.arrange(
  plotColData(sce, x="day", y="sum", colour_by="discard_adaptive") + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="day", y="detected", colour_by="discard_adaptive") + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="day", y="subsets_Mito_percent", 
              colour_by="discard_adaptive") + ggtitle("Mito percent"),
  ncol=1
)
dev.off()

## save the sce object

saveRDS(sce, 
        file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/02_qc_sce.rds")



