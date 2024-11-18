library(tidyr)
library(dplyr)
library(corrplot)
library(reshape2)
# get all quant files from /dcs04/hicks/data/sparthib/06_quantification/oarfish

file.list <- list.files("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/oarfish",
                        pattern = ".quant", full.names = TRUE)

# read in all quant files
all.files<- lapply(file.list, read.table, header = TRUE, sep = "\t")

# get the sample names
sample.names <- gsub(".quant", "", basename(file.list))

# add sample names to the dataframes
for(i in 1:length(all.files)){
  all.files[[i]]$sample <- sample.names[i]
}

# combine all quant files into one dataframe
all.files <- do.call(rbind, all.files)

all.files <- all.files |> select(sample, tname, num_reads)

#fill NA  num_reads with 0
all.files$num_reads[is.na(all.files$num_reads)] <- 0

nrow(all.files)
filtered_data <- all.files |> 
  group_by(tname) |> 
  filter(var(num_reads) != 0) |> 
  ungroup()

nrow(filtered_data)
wide_data <- filtered_data %>%
  pivot_wider(names_from = sample, values_from = num_reads, values_fill = 0)

colnames(wide_data) <- gsub("[- ]+", "_", colnames(wide_data))

# View updated column names
print(colnames(wide_data))


isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT"

# Function to load data
load_expression_data <- function(data_type) {
  counts <- read.table(file.path(isoquant_dir, paste0("OUT.", data_type, "_grouped_counts.tsv")), header = TRUE)
  tpm <- read.table(file.path(isoquant_dir, paste0("OUT.", data_type, "_grouped_tpm.tsv")), header = TRUE)
  colnames(counts) <- c("id", samples)
  colnames(tpm) <- c("id", samples)
  list(counts = counts, tpm = tpm)
}


# Function to filter rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}

# Function to prepare TPM data for PCA
prepare_tpm <- function(tpm, sample_indices) {
  selected_tpm <- tpm[, sample_indices]
  rownames(selected_tpm) <- tpm[, 1]
  remove_zero_var_rows(selected_tpm)
}

library(tibble)
samples <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
             "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
             "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")
isoform_data <- load_expression_data("transcript")
RO_isoform_tpm <- prepare_tpm(isoform_data$tpm, 2:8)
#remove isoform id version 
RO_isoform_tpm <- RO_isoform_tpm |> rownames_to_column(var = "isoform_id") |>
  mutate(isoform_id = gsub("\\..*", "", isoform_id)) |> column_to_rownames(var = "isoform_id")

colnames(RO_isoform_tpm)
colnames(RO_isoform_tpm) <- gsub("[- ]+", "_", colnames(RO_isoform_tpm))



wide_data <- wide_data |> column_to_rownames(var = "tname")
#remove version from isoform id
wide_data <- wide_data |> rownames_to_column(var = "isoform_id") |>
  mutate(isoform_id = gsub("\\..*", "", isoform_id)) |> column_to_rownames(var = "isoform_id")

matching_rows <- intersect(rownames(RO_isoform_tpm), rownames(wide_data))
matching_cols <- intersect(colnames(RO_isoform_tpm), colnames(wide_data))




##### DTE #######

library(edgeR)

oarfish_counts <- wide_data

oarfish_counts <- oarfish_counts |> select(c("EP1_BRN3B_RO","EP1_WT_hRO_2","EP1_WT_ROs_D45",
                         "H9_BRN3B_hRO_2","H9_BRN3B_RO","H9_CRX_hRO_2","H9_CRX_ROs_D45"))

targets  <- data.frame(Sample = c("EP1_BRN3B_RO","EP1_WT_hRO_2","EP1_WT_ROs_D45",
                                  "H9_BRN3B_hRO_2","H9_BRN3B_RO","H9_CRX_hRO_2","H9_CRX_ROs_D45") ,
                       Group = c("RO_D200", "RO_D100", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45"),
                       Replicate = c(1, 1, 1, 2, 2, 3 ,2),
                       stringsAsFactors = FALSE)


y <- DGEList(counts = oarfish_counts,
             samples = targets$Sample,
             group = targets$Group,
             genes = rownames(oarfish_counts))

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)



design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion



fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(D200_vs_D100 =  RO_D200 - RO_D100, 
                       D200_vs_D45 = RO_D100 - RO_D45,
                       D100_vs_D45 = RO_D100 - RO_D45,
                       levels=design)

for (i in 1:ncol(contr)){
  
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  # summary(qlf)
  is.de <- decideTests(qlf, p.value=0.05)
  summary(is.de)
  tt <- topTags(qlf,n = Inf)
  nrow(tt) #10261
  head(tt)
  
  tt$table$isoform_id <- tt$table$genes
  
  library(biomaRt)
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  
  annotLookup <- getBM(
    mart=mart,
    attributes=c( "ensembl_gene_id",
                  "ensembl_transcript_id",
                  "external_gene_name",
                  "transcript_biotype"),
    filter="ensembl_transcript_id",
    values=tt$table$isoform_id,
    uniqueRows=TRUE)
  
  colnames(annotLookup) <- c("gene_id", "isoform_id", "gene_name", "transcript_biotype")
  tt$table <- merge(tt$table, annotLookup,
                    by="isoform_id", all.x=TRUE)
  
  tt$table <- tt$table[order(tt$table$FDR),]
  tt$table$condition_1 <- names(contr[,i])[contr[,i]== 1]
  tt$table$condition_2 <- names(contr[,i])[contr[,i]== -1]
  
  file <- paste0("/users/sparthib/retina_lrs/processed_data/", colnames(contr)[i], "_oarfish_DTEs.tsv")
  readr::write_tsv(tt$table, file)
  
  
}

library(readr)
library(dplyr)
oarfish_D200_vs_D100 <- read_tsv("/users/sparthib/retina_lrs/processed_data/D200_vs_D100_oarfish_DTEs.tsv")
oarfish_D200_vs_D45 <- read_tsv("/users/sparthib/retina_lrs/processed_data/D200_vs_D45_oarfish_DTEs.tsv")
oarfish_D100_vs_D45 <- read_tsv("/users/sparthib/retina_lrs/processed_data/D100_vs_D45_oarfish_DTEs.tsv")


rbind(oarfish_D200_vs_D100, oarfish_D200_vs_D45, oarfish_D100_vs_D45) |> write_tsv("/users/sparthib/retina_lrs/processed_data/oarfish_DTEs.tsv")


oarfish_DTEs <- read_tsv("/users/sparthib/retina_lrs/processed_data/oarfish_DTEs.tsv")


oarfish_DTEs <- oarfish_DTEs |> dplyr::select(c("isoform_id", "logFC", "FDR", "condition_1", "condition_2"))


result <- oarfish_DTEs |>
  group_by(isoform_id) |>
  filter(FDR == min(FDR) & abs(logFC) == max(abs(logFC))) |>
  ungroup()   
  

result  <- result |> filter(FDR < 0.05) |> filter(abs(logFC) >= 1) 

#get isoforms with lowest FDR and highest logFC

sorted_result <- result |> 
  arrange(FDR, desc(abs(logFC)))  |> head(500)



isoquant_input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs"
switch_file <- read_tsv(file.path(isoquant_input_data_dir, "DGE_DTU_DTE.tsv"))

significant_DTEs <- switch_file |> dplyr::group_by(isoform_id ) |>
  filter(
    DTE_qval == min(DTE_qval) & 
      abs(DTE_log2FC) == max(abs(DTE_log2FC))
  ) |> filter(DTE_qval < 0.05 & abs(DTE_log2FC) >= 0.1) |>
  arrange(DTE_qval) |> 
  dplyr::select(isoform_id, gene_name) |>
  distinct() |> head(500)

# Get the isoforms that are common between the two methods

common_isoforms <- intersect(sorted_result$isoform_id, significant_DTEs$isoform_id)

isoquant_subset <- RO_isoform_tpm[matching_rows, matching_cols, drop = FALSE]
oarfish_subset <- wide_data[matching_rows, matching_cols, drop = FALSE]

isoquant_subset <- isoquant_subset[common_isoforms, ]
oarfish_subset <- oarfish_subset[common_isoforms, ]
# Calculate correlation matrix
correlation_matrix <- cor(isoquant_subset , oarfish_subset, use = "pairwise.complete.obs")

correlation_df <- melt(correlation_matrix, na.rm = TRUE)
colnames(correlation_df) <- c("Matrix1_Sample", "Matrix2_Sample", "Correlation")


# Define the desired order for rows and columns
sample_order <- c("EP1_WT_ROs_D45", "H9_CRX_ROs_D45", "H9_CRX_hRO_2", 
                  "EP1_WT_hRO_2", "H9_BRN3B_hRO_2", "H9_BRN3B_RO", "EP1_BRN3B_RO")

correlation_df$Matrix1_Sample <- factor(correlation_df$Matrix1_Sample, levels = sample_order)
correlation_df$Matrix2_Sample <- factor(correlation_df$Matrix2_Sample, levels = sample_order)



# Plot
library(ggplot2)
pdf("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/isoquant_oarfish_correlation.pdf")
p <- ggplot(correlation_df, aes(x = Matrix1_Sample, y = Matrix2_Sample, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) + 
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(-1, 1), na.value = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Correlation between Isoquant and Oarfish Significant DTE values",
       x = "Oarfish Values", y = "Isoquant Values", fill = "Correlation")
print(p)
dev.off()
