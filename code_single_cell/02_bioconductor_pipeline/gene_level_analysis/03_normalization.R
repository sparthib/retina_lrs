library(scater)
library(ggplot2)


# For each cell: Normalized count = (raw count / total counts) x scale factor
# Log-transformation after adding a pseudo-count of 1

sce <- readRDS(file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/02_qc_sce.rds")
plot_dir <- "/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/02_normalization/"
dir.create(plot_dir, showWarnings = FALSE)

## remove discarded cells
sce <- sce[, !sce$discard_fixed]


# create library size factors for all the cells
lib.sce <- librarySizeFactors(sce)
summary(lib.sce)

# > summary(lib.sce)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09874 0.59089 0.83608 1.00000 1.20335 3.90715

pdf("/users/sparthib/retina_lrs/single_cell_plots/01_sce_eda/02_normalization/01_library_size_factors.pdf")
hist(log10(lib.sce),
     xlab="Log10[Size factor]",
     col='grey80')
dev.off()

# Normalization
sce <- logNormCounts(sce,
                     size.factors = lib.sce)

saveRDS(sce, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/03_normalized_sce.rds")


