library(readr)

input_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/all_samples/OUT"

gene_counts_path <- file.path(input_dir, "OUT.gene_grouped_counts.tsv")
gene_tpm_path <- file.path(input_dir, "OUT.gene_grouped_tpm.tsv")
transcript_counts_path <- file.path(input_dir, "OUT.transcript_grouped_counts.tsv")
transcript_tpm_path <- file.path(input_dir, "OUT.transcript_grouped_tpm.tsv")

# Load the data using read_tsv
gene_counts <- read_tsv(gene_counts_path, col_names = TRUE)

# Set the row names using the first column
rownames <- gene_counts[[1]]

#remove version number from gene names
rownames <- sub("\\.\\d+$", "", rownames)

# Remove the first column (which is now just the original data, not needed)
gene_counts <- gene_counts[, -1]

# Convert to matrix
gene_counts <- as.matrix(gene_counts)

rownames(gene_counts) <- rownames
# Verify the row names and the matrix
head(gene_counts)


gene_tpm <- read_tsv(gene_tpm_path, col_names = TRUE)
rownames <- gene_tpm[[1]]
rownames <- sub("\\.\\d+$", "", rownames)
gene_tpm <- gene_tpm[, -1]
gene_tpm <- as.matrix(gene_tpm)
rownames(gene_tpm) <- rownames


# 
# transcript_counts <- fread(transcript_counts_path, header = TRUE)
# rownames(transcript_counts) <- transcript_counts[[1]]
# transcript_counts <- as.matrix(transcript_counts[, -1])
# 
# transcript_tpm <- fread(transcript_tpm_path, header = TRUE)
# rownames(transcript_tpm) <- transcript_tpm[[1]]
# transcript_tpm <- as.matrix(transcript_tpm[, -1])


# Create the SCE object
library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(
    counts = gene_counts,
    tpm = gene_tpm
    # transcript_counts = transcript_counts,
    # transcript_tpm = transcript_tpm
  )
)


# library(dplyr)
# sample_barcode <- readRDS("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/sample_barcodes.rds")
# sample_barcode$sample <- sub("_whitelist$", "", sample_barcode$sample)
# sample_barcode$sample <- sub("X", "x", sample_barcode$sample)
# sample_barcode$sample_base <- sub("_[^_]+$", "", sample_barcode$sample)
# sample_barcode$sample |> unique()
# sample_barcode <- sample_barcode |> dplyr::select(barcode, sample_base, day) |> distinct()
#  
# # > sample_barcode |> nrow()
# # [1] 4228
# 
# # Select relevant columns and remove duplicates, keeping only the first occurrence for each barcode
# sample_barcode_unique <- sample_barcode |>
#   dplyr::select(barcode, sample_base, day) |>
#   distinct(barcode, .keep_all = TRUE)

# Display the resulting dataset
# sample_barcode_unique |> nrow()
# 4185

sample_names <- sub("_.*$", "", colnames(sce) ) 
unique(sample_names)

## get unique barcodes from colnames gene_counts 


# Ensure that barcodes in sample_barcode_unique match colnames of the SCE object
sample_barcode_unique <- sample_barcode_unique[match(colnames(sce), sample_barcode_unique$barcode), ]
sample_barcode_unique <- sample_barcode_unique |>
  dplyr::mutate(sample_base = sub("^(.*-.*?)(_[^_]*).*$", "\\1", barcode))


sample_barcode_unique$sample_base |> unique()

# Add the sample information to colData of the SCE object
colData(sce) <- cbind(colData(sce), sample_barcode_unique)

# View the updated colData
colData(sce)


library(biomaRt)
# Connect to Ensembl database
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl")

gene_ids <- rownames(rowData(sce)) 


# Query the Biomart to get gene symbols
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = gene_ids,
                      mart = ensembl)

# Merge gene symbols into your rowData (using Ensembl gene IDs as row names)
rowData(sce)$gene_symbol <- gene_symbols$hgnc_symbol[match(rownames(rowData(sce)), gene_symbols$ensembl_gene_id)]

# all genes that start with "MT"
mt_genes <- grep("^MT-", rowData(sce)$gene_symbol, value = TRUE)

#get rownames of rowData that are  MT genes
mt_genes <- rownames(rowData(sce))[rowData(sce)$gene_symbol %in% mt_genes]

#### QC #####
library(scuttle)
# df <- perCellQCMetrics(sce, list(Mito = mt_genes))
# reasons <- perCellQCFilters(df)
# colSums(as.matrix(reasons))


sce <- addPerCellQCMetrics(sce , subsets=list(Mito=mt_genes))
df <- perCellQCMetrics(sce, subsets=list(Mito=mt_genes))
reasons <- perCellQCFilters(df, 
                            sub.fields=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))


# We then identify cells that are outliers for the various QC metrics, based on
# the median absolute deviation (MAD) from the median value of
# each metric across all cells. By default, we consider a value to be an outlier
# if it is more than 3 MADs from the median in the “problematic” direction.
# low_lib_size low_n_features        discard
# 161            408            408

sce$discard <- reasons$discard



######## diagnostic plot ########
library(scater)

pdf("/users/sparthib/retina_lrs/single_cell_plots/sce_qc.pdf")
gridExtra::grid.arrange(
  plotColData(sce, x="sample_base", y="sum", colour_by="discard") +  
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="sample_base", y="detected", colour_by="discard") + 
    scale_y_log10() +  ggtitle("Detected features"),
  ncol=1
)

# Close the PDF device to save the plots
dev.off()

pdf("/users/sparthib/retina_lrs/single_cell_plots/mito_high_lib_size.pdf")
plotColData(sce, x="sum", y="subsets_Mito_percent", colour_by="discard")
dev.off()


# The other option is to simply mark the low-quality cells as 
# such and retain them in the downstream analysis. 
# The aim here is to allow clusters of low-quality cells to form, 
# and then to identify and ignore such clusters during interpretation of the results. 
# This approach avoids discarding cell types that have poor values for the QC metrics, 
# deferring the decision on whether a cluster of such cells represents a genuine biological 
# state.


saveRDS(sce, "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/quantification_alternatives/01_IsoQuant/sce_files/01_post_qc_sce.rds")







