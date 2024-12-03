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
source("/users/sparthib/retina_lrs/code/04_dtu_dge_dte/bambu/helper.R")


# Define directories and common variables

plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots"

input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs"
input_gene_counts <- "gene_counts.RDS"
input_isoform_counts <- "isoform_counts.RDS"
input_isoform_cpm <- "isoform_cpm.RDS"

gene_counts <- readRDS(file.path(input_dir, input_gene_counts))
isoform_counts <- readRDS(file.path(input_dir, input_isoform_counts))
isoform_cpm <- readRDS(file.path(input_dir, input_isoform_cpm))


colnames(gene_counts)[1] <- "gene_id"

dge <- DGEList(counts=gene_counts)

## calculate norm. factors
dge <- calcNormFactors(dge)

## get normalized counts
gene_cpm <- cpm(dge)
saveRDS(gene_cpm, file.path(input_dir, "gene_cpm.RDS"))

gene_counts <- remove_zero_var_rows(gene_counts)
isoform_counts <- remove_zero_var_rows(isoform_counts)
isoform_cpm <- remove_zero_var_rows(isoform_cpm)
gene_cpm <- remove_zero_var_rows(gene_cpm)


samples <- colnames(gene_counts)
groups <- c("RO_D45", "RO_D45", "RO_D100","RO_D100", 
                        "RO_D100", "RO_D200", "RO_D200")


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

# save venn diagram as pdf 
input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs"
DGE_DTU_DTE <- read_tsv(file.path(input_data_dir, "DGE_DTU_DTE.tsv"))


gene_overlaps = new_DGE_DTE_DTU |> dplyr::select( gene_id, isoform_id,
                                                  condition_1, condition_2,DGE, DTU , DTE )

gene_overlaps <- gene_overlaps |>
  mutate(condition = case_when(
    condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45" ~ "RO_D100_vs_RO_D45",
    condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45" ~ "RO_D200_vs_RO_D45",
    condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100" ~ "RO_D200_vs_RO_D100",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))

gene_overlaps <- gene_overlaps |> mutate(DTU_RO_D100_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D200_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D200_vs_RO_D100 = ifelse(DTU == TRUE & condition == "RO_D200_vs_RO_D100", TRUE, FALSE),
                                         DGE_RO_D200_vs_RO_D100 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D100", TRUE, FALSE),
                                         DGE_RO_D200_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DGE_RO_D100_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         DTE_RO_D100_vs_RO_D45 = ifelse(DTE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         DTE_RO_D200_vs_RO_D45 = ifelse(DTE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DTE_RO_D200_vs_RO_D100 = ifelse(DTE == TRUE & condition == "RO_D200_vs_RO_D100", TRUE, FALSE))

gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_RO_D100_vs_RO_D45 = any(DTU_RO_D100_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D45 = any(DTU_RO_D200_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D100 = any(DTU_RO_D200_vs_RO_D100), 
            DGE_RO_D200_vs_RO_D100 = any(DGE_RO_D200_vs_RO_D100), 
            DGE_RO_D200_vs_RO_D45 = any(DGE_RO_D200_vs_RO_D45), 
            DGE_RO_D100_vs_RO_D45 = any(DGE_RO_D100_vs_RO_D45),
            DTE_RO_D100_vs_RO_D45 = any(DTE_RO_D100_vs_RO_D45),
            DTE_RO_D200_vs_RO_D45 = any(DTE_RO_D200_vs_RO_D45),
            DTE_RO_D200_vs_RO_D100 = any(DTE_RO_D200_vs_RO_D100)) 

rownames(gene_overlaps) <- gene_overlaps$gene_id


# Prepare data for VennDiagram
venn_data <- list(
  DTU_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D45],
  DTU_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D45],
  DTU_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D100],
  DGE_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D100],
  DGE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D45],
  DGE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D45],
  DTE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D100_vs_RO_D45],
  DTE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D200_vs_RO_D45],
  DTE_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D200_vs_RO_D100]
)


library("UpSetR")

output_plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots"

upset_path <- file.path(output_plots_dir, "DGE_DTU_ROs_upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()




plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "B_RO_D100" &
                                     condition_2 == "C_RO_D45"), "RO_D100_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
                                     condition_2 == "C_RO_D45"), "RO_D200_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
                                     condition_2 == "B_RO_D100"), "RO_D100_vs_RO_D200")



######## RetNet Plots #######

RetNet_gene_list <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                               sheet = "genes_and_locations")

genes_and_diseases <- read_excel(file.path("/users", "sparthib", "retina_lrs" ,"raw_data", "RetNet.xlsx"),
                                 sheet = "diseases_and_genes")
colnames(genes_and_diseases) <- c("disease_category", "mapped_loci",
                                  "mapped_and_identified_genes")


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

DTE_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/DTE"
DTE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/DTE_table.tsv")
DTE_table <- DTE_table |> mutate(genes = gsub("\\..*", "", genes))


#remove DTE_table isoform_id version number
DTE_table$isoform_id <- gsub("\\..*", "", DTE_table$isoform_id)

#get gene names from biomart
isoform_id <- DTE_table$isoform_id
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

results <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = DTE_table$isoform_id,
  mart = ensembl
)
results <- results |> distinct()

DTE_table <- merge(DTE_table, results, by.x = "isoform_id", by.y = "ensembl_transcript_id", all.x = TRUE)


DTE_plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/plots/DTE"
if (!dir.exists(DTE_plots_dir)) {
  dir.create(DTE_plots_dir, recursive = TRUE)
}

DTE_table$external_gene_name[is.na(DTE_table$external_gene_name)] <- "unknown"

D100_vs_D45_DTEs <- DTE_table |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")

pdf(paste0(DTE_plots_dir, "D100_vs_D45_DTEs_volcano.pdf"))
EnhancedVolcano(D100_vs_D45_DTEs,
                lab = D100_vs_D45_DTEs$external_gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-16,
                FCcutoff = 1.5)
dev.off()

D200_vs_D45_DTEs <- DTE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")
pdf(paste0(DTE_plots_dir, "D200_vs_D45_DTEs_volcano.pdf"))
EnhancedVolcano(D200_vs_D45_DTEs,
                lab = D200_vs_D45_DTEs$external_gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-16,
                FCcutoff = 1.5)
dev.off()

D200_vs_D100_DTEs <- DTE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")
pdf(paste0(DTE_plots_dir, "D200_vs_D100_DTEs_volcano.pdf"))
EnhancedVolcano(D200_vs_D100_DTEs,
                lab = D200_vs_D45_DTEs$external_gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-16,
                FCcutoff = 1.5)
dev.off()




####### DGE Volcano #########
# Define directories and read in the DGE table
DGE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/DGE_table.tsv")

# Remove version number from gene_id
DGE_table <- DGE_table |>
  mutate(gene_id = gsub("\\..*", "", gene_id))

# Get gene names from biomaRt
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = DGE_table$gene_id,
  mart = ensembl
)
results <- results |> distinct()

# Merge gene names into DGE_table
DGE_table <- merge(DGE_table, results, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# Set up directory for plots
DGE_plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots/DGE"
if (!dir.exists(DGE_plots_dir)) {
  dir.create(DGE_plots_dir, recursive = TRUE)
}

# Replace NA values in external_gene_name with "unknown"
DGE_table$external_gene_name[is.na(DGE_table$external_gene_name)] <- "unknown"

# Generate volcano plots for different comparisons
D100_vs_D45_DGEs <- DGE_table |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")

pdf(paste0(DGE_plots_dir, "D100_vs_D45_DGEs_volcano.pdf"))
EnhancedVolcano(D100_vs_D45_DGEs,
                lab = D100_vs_D45_DGEs$external_gene_name,
                x = 'logFC',
                y = 'PValue')
dev.off()

D200_vs_D45_DGEs <- DGE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")
pdf(paste0(DGE_plots_dir, "D200_vs_D45_DGEs_volcano.pdf"))
EnhancedVolcano(D200_vs_D45_DGEs,
                lab = D200_vs_D45_DGEs$external_gene_name,
                x = 'logFC',
                y = 'PValue')
dev.off()

D200_vs_D100_DGEs <- DGE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")
pdf(paste0(DGE_plots_dir, "D200_vs_D100_DGEs_volcano.pdf"))
EnhancedVolcano(D200_vs_D100_DGEs,
                lab = D200_vs_D100_DGEs$external_gene_name,
                x = 'logFC',
                y = 'PValue',
                pCutoff = 10e-16,
                FCcutoff = 1.5)
dev.off()

####### p value distribution######

# Set up directory for plots
pval_plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots/dte"
if (!dir.exists(pval_plots_dir)) {
  dir.create(pval_plots_dir, recursive = TRUE)
}

# Generate p value distribution plots for different comparisons
pdf(paste0(pval_plots_dir, "D100_vs_D45_DTEs_pval_dist.pdf"))
ggplot(D100_vs_D45_DTEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D100 vs D45 DTEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()

pdf(paste0(pval_plots_dir, "D200_vs_D45_DTEs_pval_dist.pdf"))
ggplot(D200_vs_D45_DTEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D200 vs D45 DTEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()

pdf(paste0(pval_plots_dir, "D200_vs_D100_DTEs_pval_dist.pdf"))
ggplot(D200_vs_D100_DTEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D200 vs D100 DTEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()


#pvalue plots for DGE
D100_vs_D45_DGEs <- DGE_table |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")
dge_plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots/dge"
pdf(paste0(dge_plots_dir, "D100_vs_D45_DGEs_pval_dist.pdf"))
ggplot(D100_vs_D45_DGEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D100 vs D45 DGEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()

D200_vs_D45_DGEs <- DGE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")
pdf(paste0(dge_plots_dir, "D200_vs_D45_DGEs_pval_dist.pdf"))
ggplot(D200_vs_D45_DGEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D200 vs D45 DGEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()

D200_vs_D100_DGEs <- DGE_table |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")
pdf(paste0(dge_plots_dir, "D200_vs_D100_DGEs_pval_dist.pdf"))
ggplot(D200_vs_D100_DGEs, aes(x = PValue)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "P Value Distribution for D200 vs D100 DGEs", x = "P Value", y = "Count") +
  theme_minimal()
dev.off()





