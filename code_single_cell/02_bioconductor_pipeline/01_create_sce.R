library(scuttle)
library(data.table)
library(Matrix)
library(rtracklayer)
library(scater)
library(sessioninfo)
library(bluster)
library(scran)


sample_names <- c( "10x_D100-EP1_A1",
                   "10x_D100-EP1_A2",
                   "10x_D100-EP1_B1",
                   "10x_D100-EP1_B2",
                   "10x_D200-EP1-1_A1",
                   "10x_D200-EP1-1_A2",
                   "10x_D200-EP1-1_B1",
                   "10x_D200-EP1-1_B2",
                   "10x_D200-EP1-2_A1",
                   "10x_D200-EP1-2_A2",
                   "10x_D200-EP1-2_B1",
                   "10x_D200-EP1-2_B2")

sample_names <- "10x_D100-EP1_A1_downsampled"

# gene_file <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output",
#                        paste0(sample, "/gene_count.csv"))
# gene_mat <- scuttle::readSparseCounts(gene_file,sep = ",")
# # assign gene counts to assays slot
# assay(sample_sce, "gene_counts",withDimnames=FALSE) <- gene_mat 

load_sce <- function(sample) {
  # Replace invalid characters (e.g., '-') with '_'
  sanitized_sample <- gsub("-", "_", sample)
  sce_file <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/", paste0(sample, "/sce.rds"))
  sample_sce <- readRDS(sce_file)
  colData(sample_sce)$sample <- sanitized_sample
  # Dynamically create a variable named after the sanitized sample + "_sce"
  assign(paste0("sce_", sanitized_sample), sample_sce, envir = .GlobalEnv)
}

sanitized_sample_names <- gsub("-", "_", sample_names)

sce_list <- sample_names |> lapply(load_sce)

#add sample name to sce
sce_list <- Map(function(x, sample_name) {
  colData(x)$sample_replicate <- sample_name
  x
}, sce_list, sanitized_sample_names)

non_replicate_names <- c("10x_D100_EP1","10x_D100_EP1","10x_D100_EP1","10x_D100_EP1",
                         "10x_D200_EP1_1","10x_D200_EP1_1","10x_D200_EP1_1","10x_D200_EP1_1",
                         "10x_D200_EP1_2","10x_D200_EP1_2","10x_D200_EP1_2","10x_D200_EP1_2")

sce_list <- Map(function(x, sample_name) {
  colData(x)$sample_replicate <- sample_name
  x
}, sce_list, sanitized_sample_names)

sce_list <- Map(function(x, sample_name) {
  colData(x)$sample <- sample_name
  x
}, sce_list, non_replicate_names)


input_dimensions <- lapply(sce_list, function(x) {
  print(dim(x))
})

#### #####
# > input_dimensions
# $`10x_D100_EP1_A1`
# [1] 19098  2390
# 
# $`10x_D100_EP1_A2`
# [1] 2426 2287
# 
# $`10x_D100_EP1_B1`
# [1] 25329  2419
# 
# $`10x_D100_EP1_B2`
# [1] 8344 2343
# 
# $`10x_D200_EP1_1_A1`
# [1] 19975   987
# 
# $`10x_D200_EP1_1_A2`
# [1] 5486  902
# 
# $`10x_D200_EP1_1_B1`
# [1] 25686   992
# 
# $`10x_D200_EP1_1_B2`
# [1] 4840  903
# 
# $`10x_D200_EP1_2_A1`
# [1] 15095   802
# 
# $`10x_D200_EP1_2_A2`
# [1] 1693  614
# 
# $`10x_D200_EP1_2_B1`
# [1] 15328   809
# 
# $`10x_D200_EP1_2_B2`
# [1] 1959  619

names(input_dimensions) <- sanitized_sample_names
#### #####

#convert gene id to gene symbol using bioMart
library(biomaRt)

gene_id_to_symbol <- function(gene_id) {
  gene_id <- na.omit(gene_id)  # Remove NAs
  gene_id <- as.character(gene_id)  # Convert to character vector
  gene_id <- unlist(strsplit(gene_id, "\\."))  # Flatten list from strsplit
  
  # Initialize Ensembl Mart
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Fetch gene symbols using getBM
  gene_symbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = gene_id,
                       mart = mart)
  
  return(gene_symbol)
}



### Add QC metrics 
sce_list <- lapply(sce_list, function(x) {
  
  # Check for cells with all zeros
  x <- x[, colSums(counts(x)) > 0]
  x <- scuttle::logNormCounts(x)
  
  # Get gene IDs from rowData
  gene_ids <- rowData(x)$gene_id
  
  # Use gene_id_to_symbol to convert gene IDs to symbols
  gene_symbols <- gene_id_to_symbol(gene_ids)
  
  # Add clean gene IDs without the dot extensions to rowData
  rowData(x)$gene_id_clean <- unlist(lapply(rowData(x)$gene_id, function(id) strsplit(id, "\\.")[[1]][1]))
  
  # Merge gene_symbols with rowData, ensuring to keep all genes
  merged_data <- merge(as.data.frame(rowData(x)), gene_symbols, 
                       by.x = "gene_id_clean", by.y = "ensembl_gene_id", 
                       all.x = TRUE, all.y = FALSE)
  
  #keep only unique rows after merging 
  merged_data <- unique(merged_data)
  # Update rowData(x) with merged data
  rowData(x) <- DataFrame(merged_data)
  
  # Append "_MT" to rownames of x for mitochondrial genes
  mt_genes <- grep("MT", rowData(x)$external_gene_name)
  rownames(x)[mt_genes] <- paste0(rownames(x)[mt_genes], "_MT")
  
  # Identify mitochondrial genes and add per-cell QC metrics
  is.mito <- grepl("^MT-", rownames(x))
  x <- scuttle::addPerCellQC(x, subsets = list(Mito = is.mito))
  
  # Add per-feature QC metrics
  x <- scuttle::addPerFeatureQC(x)
  
  x
})

### get sum of counts 

sum_counts <- lapply(sce_list, function(x) {
  sum_counts <- sum(colData(x)$sum)
  sum_counts
})

names(sum_counts) <- sanitized_sample_names
# $`10x_D100_EP1_A1`
# [1] 5265807
# 
# $`10x_D100_EP1_A2`
# [1] 157188
# 
# $`10x_D100_EP1_B1`
# [1] 9000247
# 
# $`10x_D100_EP1_B2`
# [1] 1047044
# 
# $`10x_D200_EP1_1_A1`
# [1] 4835027
# 
# $`10x_D200_EP1_1_A2`
# [1] 493819
# 
# $`10x_D200_EP1_1_B1`
# [1] 6405364
# 
# $`10x_D200_EP1_1_B2`
# [1] 394517
# 
# $`10x_D200_EP1_2_A1`
# [1] 2962709
# 
# $`10x_D200_EP1_2_A2`
# [1] 85119
# 
# $`10x_D200_EP1_2_B1`
# [1] 3063850
# 
# $`10x_D200_EP1_2_B2`
# [1] 109183


### Identifying low-quality cells 

library(tibble)
# Combine all colData into one data frame
all_colData <- rbindlist(lapply(sce_list, function(x) as_tibble(colData(x))))


range(all_colData$subsets_Mito_percent)

range(all_colData$sum)
range(all_colData$detected)
# > range(all_colData$sum)
# [1]     1 36150
# > range(all_colData$detected)
# [1]    1 6667


sce_list <- lapply(sce_list, function(x) {
  x$discard <- perCellQCFilters(colData(x), 
                                 sub.fields="subsets_Mito_percent")$discard
  x
})

post_qc_dimensions <- lapply(sce_list, function(x) {
  print(dim(x))
})

names(post_qc_dimensions) <- sanitized_sample_names

# The other option is to simply mark the low-quality cells as such and retain them 
# in the downstream analysis. The aim here is to allow clusters of low-quality cells to
# form, and then to identify and ignore such clusters during interpretation of the results.
# This approach avoids discarding cell types that have poor values for the QC metrics,
# deferring the decision on whether a cluster of such cells represents a genuine biological
# state.

saveRDS(sce_list,
     file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/transcriptome/01_quality_controlled/sce_list.rds")

sce_list <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/transcriptome/01_quality_controlled/sce_list.rds")


colnames(sce_list[[1]])

colnames(sce_list[[2]])


#check common names between colnames(sce_list[[1]]) and colnames(sce_list[[2]])

common_names <- intersect(colnames(sce_list[[1]]), colnames(sce_list[[2]]))





                          
                          
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
