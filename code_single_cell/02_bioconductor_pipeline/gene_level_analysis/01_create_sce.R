## 3.4.1 from tabular formats 
library(readr)
library(scuttle)
library(ggplot2)
library(biomaRt)

input_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/all_samples/OUT"
gene_counts_path <- file.path(input_dir, "OUT.gene_grouped_counts.tsv")

# Load the data using read_tsv
gene_counts <- readSparseCounts(gene_counts_path)

#remove version number from gene names
rownames(gene_counts) <- gsub("\\.\\d+$", "", rownames(gene_counts))

#only keep protein coding genes 

# Connect to Ensembl
mart <- biomaRt::useEnsembl(biomart = "genes", 
                            dataset = "hsapiens_gene_ensembl")

# Retrieve gene IDs for protein-coding genes
protein_coding_genes <- biomaRt::getBM(
  attributes = c("ensembl_gene_id"),
  filters = "biotype",
  values = "protein_coding",
  mart = mart
)

# Convert to vector
protein_coding_gene_list <- protein_coding_genes$ensembl_gene_id

# Filter the counts
gene_counts <- gene_counts[rownames(gene_counts) %in% protein_coding_gene_list, ]


### create sce ####
sce <- SingleCellExperiment(assays = list(counts = gene_counts))

# samples 
sample_names  <- sub("_(?!.*_).*", "", colnames(gene_counts), perl = TRUE)
colData(sce)$sample <- sample_names
sample_names |> unique()
assay(sce) |> dim()

colData(sce)$day <- ifelse(grepl("D200", 
                                 colData(sce)$sample), 
                           "D200", "D100") 


colData(sce)$sample |> table()
# 10x_D100-EP1 10x_D200-EP1-1 10x_D200-EP1-2 
# 2420            993            815 

# D100 D200 
# 2420 1808 

# > assay(sce) |> dim()
# [1] 20087  4228

#add gene symbols to rowData
gene_symbols <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(sce),
  mart = mart
)

rownames(gene_symbols) <- gene_symbols$ensembl_gene_id
rowData(sce)$symbol <- gene_symbols[rownames(sce), "hgnc_symbol"]

  
# save sce
saveRDS(sce, file = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/gene_level/01_preqc_sce.rds")








