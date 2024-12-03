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


# Define directories and common variables

plots_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/plots"

input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs"
input_gene_counts <- "gene_counts_ROs.RDS"
input_isoform_counts <- "isoform_counts_ROs.RDS"
input_isoform_cpm <- "isoform_cpm_ROs.RDS"

gene_counts <- readRDS(file.path(input_dir, input_gene_counts))
isoform_counts <- readRDS(file.path(input_dir, input_isoform_counts))
isoform_cpm <- readRDS(file.path(input_dir, input_isoform_cpm))


colnames(gene_counts)[1] <- "gene_id"

dge <- DGEList(counts=gene_counts)

## calculate norm. factors
dge <- calcNormFactors(dge)

## get normalized counts
gene_cpm <- cpm(dge)

gene_counts <- remove_zero_var_rows(gene_counts)
isoform_counts <- remove_zero_var_rows(isoform_counts)
isoform_cpm <- remove_zero_var_rows(isoform_cpm)
gene_cpm <- remove_zero_var_rows(gene_cpm)


samples <- colnames(gene_counts)
groups <- c("RO_D45", "RO_D45", "RO_D100","RO_D100", 
                        "RO_D100", "RO_D200", "RO_D200")


# Function to filter rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
}



# PCA and plotting function
pca_plots_dir <- file.path(plots_dir, "pca")
if (!dir.exists(pca_plots_dir)) {
  dir.create(pca_plots_dir, recursive = TRUE)
}
plot_pca <- function(tpm, samples, groups, output_name, output_dir) {
  pc <- prcomp(t(log2(tpm + 1)), scale = TRUE)
  pcr <- data.frame(pc$x[, 1:2], row.names = samples)  # PC scores
  
  # PCA plot
  p <- ggplot(pcr, aes(PC1, PC2, color = groups)) +  
    geom_point(size = 2) +
    theme_bw() +
    ggtitle(paste("PCA on", output_name, "Expression")) +
    xlab(paste("PC1 (", round(pc$sdev[1]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +
    ylab(paste("PC2 (", round(pc$sdev[2]^2 / sum(pc$sdev^2) * 100, 2), "%)")) +
    geom_label(aes(label = rownames(pcr)), size = 2, fill = "white", alpha = 0.7)
  
  # Scree plot
  var_explained <- pc$sdev[1:10]^2 / sum(pc$sdev^2)
  var_explained_df <- data.frame(
    Principal_Component = factor(1:length(var_explained)),
    Variance_Explained = var_explained
  )
  
  q <- ggplot(var_explained_df, aes(x = Principal_Component, y = Variance_Explained)) + 
    geom_bar(stat = "identity", fill = "skyblue") +
    xlab("Principal Component") +  
    ylab("Variance Explained") +  
    ggtitle("Scree Plot") +
    ylim(0, 1) +
    theme_minimal()
  
  # Save plots
  pdf(file.path(output_dir, paste0(output_name, "_level_pca.pdf")), width = 12, height = 6)
  print(p + q)
  dev.off()
}

# Prepare data and plot for isoforms
pca_plots_dir <- file.path(plots_dir, "pca")
plot_pca(isoform_cpm, samples, groups, "isoform", pca_plots_dir)
plot_pca(gene_cpm, samples, groups, "gene", pca_plots_dir)

######## HEATMAPS ########

heatmap_plots_dir <- file.path(plots_dir, "heatmaps")
if (!dir.exists(heatmap_plots_dir)) {
  dir.create(heatmap_plots_dir, recursive = TRUE)
}

input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs"

isoformFeatures <- read_tsv(file.path(input_data_dir, "isoformFeatures.tsv"))
isoformFeatures$isoform_switch_q_value




plot_heatmap <- function(input_data_dir, quant_name, compare, tmm, groups, output_plots_dir, type) {
  # Read the DEXSeqSwitchList file
  switch_file <- read_tsv(file.path(input_data_dir, "isoformFeatures.tsv"))
  if(type == "dtu"){
    significant_DTUs <- switch_file |> dplyr::group_by(isoform_id ) |>
      filter(isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1) |>
      arrange(isoform_switch_q_value) |> 
      dplyr::select(isoform_id, gene_name) |>
      distinct() |> 
      head(n = 200)
    
    # Keep only tmm rows that are in DTU_isoforms based on its rownames
    tpm <- RO_isoform_tpm
    tpm <- as.data.frame(tpm)
    tpm$isoform_id <- rownames(tpm)
    #remove version number 
    tpm$isoform_id <- gsub("\\..*", "", tpm$isoform_id)
    significant_DTUs$isoform_id <- gsub("\\..*", "", significant_DTUs$isoform_id)
    TPM_significant_isoforms <- tpm |>
      inner_join(significant_DTUs, by = "isoform_id") }
  
  else if(type == "dte"){ 
    
    }
  
  
  
  
  
  
  # Create a new column with gene name and isoform name
  TPM_significant_isoforms <- TPM_significant_isoforms |> 
    mutate(gene_isoform = paste(gene_name, isoform_id, sep = "_")) |>
    column_to_rownames(var = "gene_isoform") |>
    dplyr::select(-isoform_id, -gene_name)
  
  
  # Convert to matrix and scale
  TPM_matrix <- as.matrix( TPM_significant_isoforms )
  TPM_matrix <- t(scale(t(TPM_matrix), center = TRUE, scale = TRUE))
  
  # Color function for heatmap
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Heatmap annotation
  compare = "ROs"
  groups = RO_isoquant_groups
  ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left")
  if (compare == "FT_vs_RGC"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("FT" = "lightgreen", "RGC" = "brown")
                            ))
  }else if(compare == "ROs"){
    ha <- HeatmapAnnotation(type = groups, annotation_name_side = "left",
                            col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen")
                            ))
  }
  
  
  # Generate and draw heatmap
  pdf(file.path(heatmap_plots_dir, paste0(compare, "_DTU_heatmap.pdf")))
  ht_list <- Heatmap(TPM_matrix , name = paste0("scaled TPM expression of top", nrow(TPM_matrix), " DTU isoforms"),  row_km = 5,
                     col = col_fun, top_annotation = ha, show_row_names = TRUE,
                     show_column_names = TRUE, row_title = "isoforms",
                     row_names_gp = gpar(fontsize = 2.2),
                     column_names_gp = gpar(fontsize = 5),
                     show_row_dend = TRUE, show_column_dend = TRUE)
  draw(ht_list)
  dev.off()
  
}

plot_heatmap(input_data_dir, "isoquant", "ROs", ROs_isoquant_tpm, ROs_isoquant_groups, heatmap_plots_dir)


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



plot_barplot <- function(df, condition){ 
  
  dge_overlaps = df |> group_by(gene_id) |> dplyr::select( -DTE) |> 
    summarise(DGE=any(DGE), DTU=any(DTU))  
  
  dte_overlaps = df |> group_by(gene_id) |> dplyr::select(-DGE) |> 
    summarise(DTE = any(DTE), DTU=any(DTU)) 
  
  dge_contingency_table <- as_tibble(xtabs(~ DGE + DTU, data = dge_overlaps))
  dge_contingency_table <- plyr::ddply(dge_contingency_table, ~DGE, transform, Prop = n / sum(n))
  
  
  dge_bar <- ggplot(dge_contingency_table, aes(x = DGE, y = n, fill = DTU)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "DGE and DTU genes", x = "DGE", y = "Count") + 
    scale_y_continuous(labels = scales::comma) +  # Display y-axis as percentage
    theme_minimal() +
    geom_text(aes(label = n), 
              position = position_stack(vjust = 0.5), 
              size = 2)
  
  dte_contingency_table <- xtabs(~ DTE + DTU, data = dte_overlaps)
  dte_contingency_table <- plyr::ddply(as_tibble(dte_contingency_table), ~DTE, transform, Prop = n / sum(n))
  
  dte_bar <- ggplot(dte_contingency_table, aes(x = DTE, y = n, fill = DTU)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "DTE and DTU genes", x = "DTE", y = "Count") + 
    scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
    theme_minimal() +
    geom_text(aes(label = n), 
              position = position_stack(vjust = 0.5), 
              size = 2)
  
  library(patchwork)
  p <- dge_bar + dte_bar + plot_annotation(title = condition)
  ggsave(path = output_plots_dir, 
         device = "pdf", 
         plot = p, 
         filename = paste0(condition, "df_barplot.pdf"))
  
}

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

######## GO Analysis #########





