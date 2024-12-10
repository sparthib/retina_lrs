library("ggplot2")
library("DGEobj.utils")
library("patchwork")
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library("ggVennDiagram")
library(readxl)
library(stringr)
library(EnhancedVolcano)
library(edgeR)

### import functions from helper.R
source("/users/sparthib/retina_lrs/code/04_dtu_dge_dte/helper.R")

# Define directories and common variables
method <- "bambu"
comparison <- "FT_vs_RGC"
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "plots")

input_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                       method,comparison) 
input_gene_counts <- "gene_counts.RDS"
input_isoform_counts <- "isoform_counts.RDS"
input_isoform_cpm <- "isoform_cpm.RDS"

gene_counts <- readRDS(file.path(input_dir, input_gene_counts))
isoform_counts <- readRDS(file.path(input_dir, input_isoform_counts))
isoform_cpm <- readRDS(file.path(input_dir, input_isoform_cpm))


# colnames(gene_counts)[1] <- "gene_id"
# 
# dge <- DGEList(counts=gene_counts)
# 
# ## calculate norm. factors
# dge <- calcNormFactors(dge)
# 
# ## get normalized counts
# gene_cpm <- cpm(dge)
# saveRDS(gene_cpm, file.path(input_dir, "gene_cpm.RDS"))
gene_cpm <- readRDS(file.path(input_dir, "gene_cpm.RDS"))

gene_counts <- remove_zero_var_rows(gene_counts)
isoform_counts <- remove_zero_var_rows(isoform_counts)
isoform_cpm <- remove_zero_var_rows(isoform_cpm)
gene_cpm <- remove_zero_var_rows(gene_cpm)


samples <- colnames(gene_counts)
groups <- c("FT", "FT", "RGC", "RGC")


# PCA and plotting function
pca_plots_dir <- file.path(plots_dir, "pca")
if (!dir.exists(pca_plots_dir)) {
  dir.create(pca_plots_dir, recursive = TRUE)
}

# Prepare data and plot for isoforms
pca_plots_dir <- file.path(plots_dir, "pca")
plot_pca(isoform_cpm, samples, groups, "isoform", pca_plots_dir)
plot_pca(gene_cpm, samples, groups, "gene", pca_plots_dir)



####### Upset Plot ########

input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu",
                            method, comparison)
DGE_DTU_DTE <- read_tsv(file.path(input_data_dir,
                                  "DGE_DTE_DTU.tsv"))


gene_overlaps = DGE_DTU_DTE |> dplyr::select( gene_id, isoform_id,
                                                  condition_1, condition_2,DGE,DTU ,DTE )

gene_overlaps <- gene_overlaps |>
  mutate(condition = case_when(
    condition_1 == "FT" & condition_2 == "RGC" ~ "FT_vs_RGC",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))

gene_overlaps <- gene_overlaps |> mutate(DTU_FT_vs_RGC = ifelse(DTU == TRUE & condition == "FT_vs_RGC", TRUE, FALSE),
                                         DGE_FT_vs_RGC = ifelse(DGE == TRUE & condition == "FT_vs_RGC", TRUE, FALSE),
                                         DTE_FT_vs_RGC = ifelse(DTE == TRUE & condition == "FT_vs_RGC", TRUE, FALSE))

gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 
gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_FT_vs_RGC = any(DTU_FT_vs_RGC), 
            DGE_FT_vs_RGC = any(DGE_FT_vs_RGC), 
            DTE_FT_vs_RGC = any(DTE_FT_vs_RGC))

rownames(gene_overlaps) <- gene_overlaps$gene_id


# Prepare data for VennDiagram
venn_data <- list(
  DTU_FT_vs_RGC = rownames(gene_overlaps)[gene_overlaps$DTU_FT_vs_RGC],
  DGE_FT_vs_RGC  = rownames(gene_overlaps)[gene_overlaps$DGE_FT_vs_RGC],
  DTE_FT_vs_RGC = rownames(gene_overlaps)[gene_overlaps$DTE_FT_vs_RGC]
)


library("UpSetR")

upset_path <- file.path(plots_dir, "DGE_DTU_ROs_upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()


plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "FT" &
                                     condition_2 == "RGC"), "FT_vs_RGC",
             plots_dir)

######## RetNet Plots #######

RetNet_gene_list <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")

genes_and_diseases <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")

significant_DTUs <- DGE_DTU_DTE |> filter(DTU)

intersect(significant_DTUs$gene_name, RetNet_gene_list$Symbol)

# Initialize an empty list to store results
results_list <- list()
for (gene in intersect(significant_DTUs$gene_name, RetNet_gene_list$Symbol)) {
  # Filter the data
  result <- genes_and_diseases |>
    dplyr::filter(str_detect(mapped_and_identified_genes, gene))
  
  # Append results to the list
  results_list[[gene]] <- result$disease_category
}

# Convert the list to a dataframe
results_df <- data.frame(
  Gene = names(results_list),
  Disease_Category = I(results_list),
  stringsAsFactors = FALSE
)

##### Volcano plots ###### 

DTE_table <- read_tsv(file.path(input_data_dir, "DTE_table.tsv"))

#remove DTE_table isoform_id version number
DTE_table$isoform_id <- gsub("\\..*", "", DTE_table$isoform_id)

library(biomaRt)

us_mart <- useEnsembl(biomart = "ensembl", mirror = "asia")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)


results <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = DTE_table$isoform_id,
  mart = mart 
)
results <- results |> distinct()

DTE_table <- merge(DTE_table, results, by.x = "isoform_id", by.y = "ensembl_transcript_id", all.x = TRUE)


DTE_plots_dir <- file.path(plots_dir, "DTE")
if (!dir.exists(DTE_plots_dir)) {
  dir.create(DTE_plots_dir, recursive = TRUE)
}

DTE_table$external_gene_name[is.na(DTE_table$external_gene_name)] <- "unknown"


pdf(paste0(DTE_plots_dir, "FT_vs_RGC_DTEs_volcano.pdf"))
EnhancedVolcano(DTE_table,
                lab = DTE_table$external_gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-16,
                FCcutoff = 1.5)
dev.off()


####### DGE Volcano #########
# Define directories and read in the DGE table
DGE_table <- read_tsv(file.path(input_data_dir, "DGE_table.tsv"))
# Remove version number from gene_id
DGE_table <- DGE_table |>
  mutate(gene_id = gsub("\\..*", "", gene_id))

# Get gene names from biomaRt
library(biomaRt)
results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = DGE_table$gene_id,
  mart = mart
)
results <- results |> distinct()

# Merge gene names into DGE_table
DGE_table <- merge(DGE_table, results, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# Set up directory for plots
DGE_plots_dir <- file.path("DGE")
if (!dir.exists(DGE_plots_dir)) {
  dir.create(DGE_plots_dir, recursive = TRUE)
}

# Replace NA values in external_gene_name with "unknown"
DGE_table$external_gene_name[is.na(DGE_table$external_gene_name)] <- "unknown"

# Generate volcano plots for different comparisons

pdf(paste0(DGE_plots_dir, "FT_vs_RGC_volcano.pdf"))
EnhancedVolcano(DGE_table,
                lab = DGE_table$external_gene_name,
                x = 'logFC',
                y = 'PValue')
dev.off()


####### p value distribution######
# Set up directory for plots


pdf(paste0(DTE_plots_dir , "FT_vs_RGC_DTEs_pval_dist.pdf"))
ggplot(DTE_table, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for FT vs RGC DTEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()


pdf(paste0(DGE_plots_dir , "FT_vs_RGC_DGEs_pval_dist.pdf"))
ggplot(DGE_table, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for FT vs RGC DGEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()



generate_volcano_plots <- function(input_data_dir, plots_dir, comparison, table_type, conditions) {
  # Set up the file path and read the appropriate table
  table_file <- file.path(input_data_dir, paste0(table_type, "_table.tsv"))
  analysis_table <- read_tsv(table_file)
  
  # Remove version number and harmonize identifiers
  if (table_type == "DTE") {
    analysis_table$isoform_id <- ifelse(
      grepl("^ENST", analysis_table$isoform_id),  # Check if isoform_id starts with "ENST"
      gsub("\\..*", "", analysis_table$isoform_id),  # Remove everything after the first dot
      analysis_table$isoform_id  # Keep other isoform_id values unchanged
    )
  }
  
  # Harmonize gene IDs for DGE
  if (table_type == "DGE") {
    analysis_table <- analysis_table |> 
      mutate(gene_id = gsub("\\..*", "", gene_id))
  }
  
  # Use biomaRt to map IDs to gene names
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  id_column <- ifelse(table_type == "DTE", "isoform_id", "gene_id")
  filters <- ifelse(table_type == "DTE", "ensembl_transcript_id", "ensembl_gene_id")
  attributes <- switch(
    table_type,
    "DTE" = c("ensembl_transcript_id", "external_gene_name"),
    "DGE" = c("ensembl_gene_id", "external_gene_name")
  )
  results <- getBM(
    attributes = attributes,
    filters = filters,
    values = analysis_table[[id_column]],
    mart = ensembl
  )
  
  results <- results |> distinct()
  colnames(results) <-  c(id_column, "external_gene_name")
  
  # Merge gene names into the table
  analysis_table <- merge(analysis_table, results, by.x = id_column, 
                          by.y = id_column, all.x = TRUE)
  analysis_table$external_gene_name[is.na(analysis_table$external_gene_name)] <- "unknown"
  
  # Set up directory for plots
  plots_dir <- file.path(input_data_dir, "plots", "volcano")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # Generate volcano plots for specified conditions
  for (condition_pair in conditions) {
    condition_1 <- condition_pair[[1]]
    condition_2 <- condition_pair[[2]]
    filtered_table <- analysis_table |> filter(condition_1 == !!condition_1 & 
                                                 condition_2 == !!condition_2)
    
    output_file <- file.path(plots_dir, paste0(condition_1, "_vs_", condition_2, "_", table_type, "_volcano.pdf"))
    pdf(output_file)
    EnhancedVolcano(
      filtered_table,
      lab = filtered_table$external_gene_name,
      x = "logFC",
      y = "PValue",
      pCutoff = 10e-16,
      FCcutoff = 1.5
    )
    dev.off()
  }
}

# Example usage:
method <- "bambu"
comparison <- "ROs"

input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison) 

plots_dir <- file.path(input_data_dir, "plots", "volcano")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Define condition pairs for DTE and DGE
conditions <- list(
  c("A_RO_D200", "C_RO_D45"),
  c("B_RO_D100", "C_RO_D45"),
  c("A_RO_D200", "B_RO_D100")
)

# Generate plots for DTE
generate_volcano_plots(input_data_dir, plots_dir, 
                       comparison, "DTE", conditions)

# Generate plots for DGE
generate_volcano_plots(input_data_dir, plots_dir,
                       comparison, "DGE", conditions)



