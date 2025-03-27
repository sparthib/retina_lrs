####### Upset Plot ########
library(readr)
library(dplyr)
library(UpSetR)
library(tidyr)
library(ggplot2)
# save venn diagram as pdf 
method = "bambu"
comparison = "ROs"
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison, "protein_coding")
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison,"protein_coding", "plots", "upset")
if (!dir.exists(plots_dir)){
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
}

DGE_DTU_DTE <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))

gene_overlaps = DGE_DTU_DTE |> dplyr::select( gene_id, isoform_id,
                                              condition_1, condition_2,DGE, 
                                              DTU , DTE )

gene_overlaps <- gene_overlaps |>
  mutate(condition = case_when(
    condition_1 == "Stage_1" & condition_2 == "Stage_2" ~ "Stage_1_vs_Stage_2",
    condition_1 == "Stage_1" & condition_2 == "Stage_3" ~ "Stage_1_vs_Stage_3",
    condition_1 == "Stage_2" & condition_2 == "Stage_3" ~ "Stage_2_vs_Stage_3",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))

nrow(gene_overlaps)

gene_overlaps <- gene_overlaps |> mutate(DTU_Stage_1_vs_Stage_2 = ifelse(DTU == TRUE & condition == "Stage_1_vs_Stage_2", TRUE, FALSE),
                                         DTU_Stage_1_vs_Stage_3 = ifelse(DTU == TRUE & condition == "Stage_1_vs_Stage_3", TRUE, FALSE),
                                         DTU_Stage_2_vs_Stage_3 = ifelse(DTU == TRUE & condition == "Stage_2_vs_Stage_3", TRUE, FALSE),
                                         DGE_Stage_2_vs_Stage_3 = ifelse(DGE == TRUE & condition == "Stage_2_vs_Stage_3", TRUE, FALSE),
                                         DGE_Stage_1_vs_Stage_3 = ifelse(DGE == TRUE & condition == "Stage_1_vs_Stage_3", TRUE, FALSE),
                                         DGE_Stage_1_vs_Stage_2 = ifelse(DGE == TRUE & condition == "Stage_1_vs_Stage_2", TRUE, FALSE),
                                         DTE_Stage_1_vs_Stage_2 = ifelse(DTE == TRUE & condition == "Stage_1_vs_Stage_2", TRUE, FALSE),
                                         DTE_Stage_1_vs_Stage_3 = ifelse(DTE == TRUE & condition == "Stage_1_vs_Stage_3", TRUE, FALSE),
                                         DTE_Stage_2_vs_Stage_3 = ifelse(DTE == TRUE & condition == "Stage_2_vs_Stage_3", TRUE, FALSE))


gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_Stage_1_vs_Stage_2 = any(DTU_Stage_1_vs_Stage_2), 
            DTU_Stage_1_vs_Stage_3 = any(DTU_Stage_1_vs_Stage_3), 
            DTU_Stage_2_vs_Stage_3 = any(DTU_Stage_2_vs_Stage_3), 
            DGE_Stage_2_vs_Stage_3 = any(DGE_Stage_2_vs_Stage_3),
            DGE_Stage_1_vs_Stage_3 = any(DGE_Stage_1_vs_Stage_3),
            DGE_Stage_1_vs_Stage_2 = any(DGE_Stage_1_vs_Stage_2),
            DTE_Stage_1_vs_Stage_2 = any(DTE_Stage_1_vs_Stage_2),
            DTE_Stage_1_vs_Stage_3 = any(DTE_Stage_1_vs_Stage_3),
            DTE_Stage_2_vs_Stage_3 = any(DTE_Stage_2_vs_Stage_3)
            ) 

readr::write_tsv(gene_overlaps, file.path(plots_dir, "gene_overlaps.tsv"))

rownames(gene_overlaps) <- gene_overlaps$gene_id

nrow(gene_overlaps)

table(gene_overlaps$DTU_Stage_1_vs_Stage_2,  
      gene_overlaps$DGE_Stage_1_vs_Stage_2)
table(gene_overlaps$DTU_Stage_1_vs_Stage_3,  
      gene_overlaps$DGE_Stage_2_vs_Stage_3)
table(gene_overlaps$DTU_Stage_2_vs_Stage_3,
      gene_overlaps$DGE_Stage_1_vs_Stage_3)
 

table(gene_overlaps$DTE_Stage_1_vs_Stage_2,  
      gene_overlaps$DGE_Stage_1_vs_Stage_2)
table(gene_overlaps$DTE_Stage_1_vs_Stage_3,
      gene_overlaps$DGE_Stage_1_vs_Stage_3)
table(gene_overlaps$DTE_Stage_2_vs_Stage_3,
      gene_overlaps$DGE_Stage_2_vs_Stage_3)



# Load required libraries

# Prepare data for VennDiagram
venn_data <- list(
  DTU_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_2],
  DTU_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_3],
  DTU_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_2_vs_Stage_3] ,#,
  DGE_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_2_vs_Stage_3],
  DGE_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_1_vs_Stage_3],
  DGE_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_1_vs_Stage_2],
  DTE_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_1_vs_Stage_2],
  DTE_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_1_vs_Stage_3],
  DTE_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_2_vs_Stage_3]
)

upset_path <- file.path(plots_dir, "Genes_upset_full.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()




# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "Stage_2" &
#                                      condition_2 == "Stage_1"), "Stage_1_vs_Stage_2",
#              plots_dir)
# 
# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "Stage_3" &
#                                      condition_2 == "Stage_1"), "Stage_1_vs_Stage_3",
#              plots_dir)
# 
# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "Stage_3" &
#                                      condition_2 == "Stage_2"), "RO_D100_vs_RO_D200",
#              plots_dir)

# Create the data frame
DTE_bar <- data.frame(
  Comparison = c("Stage1_vs_Stage2", "Stage1_vs_Stage3", "Stage2_vs_Stage3"),
  FALSE_FALSE = c(8377, 6592, 8202),  # non-DGE & non-DTE
  FALSE_TRUE = c(526, 869, 706),      # non-DGE & DTE
  TRUE_FALSE = c(275, 809, 252),      # DGE & non-DTE
  TRUE_TRUE = c(223, 1131, 241)       # DGE & DTE
)

# Reshape data into long format
data_long <- DTE_bar %>%
  pivot_longer(cols = -Comparison, names_to = "Category", values_to = "Count") %>%
  mutate(DGE = ifelse(Category %in% c("TRUE_FALSE", "TRUE_TRUE"), "DGE", "Non-DGE"),
         DTE = ifelse(Category %in% c("FALSE_TRUE", "TRUE_TRUE"), "DTE", "Non-DTE"))

# Plot stacked bar chart
stacked_DTE_path <- file.path(plots_dir, "stacked_DTE.pdf")
pdf(stacked_DTE_path, width = 10, height = 6)
ggplot(data_long, aes(x = DGE, y = Count, fill = DTE)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Comparison) +
  labs(title = "Stacked Bar Plot of DGE and DTE Overlaps",
       x = "DGE Status", y = "Count") +
  theme_minimal()
dev.off()


# Create the data frame for DTU vs DGE
DTU_bar <- data.frame(
  Comparison = c("Stage1_vs_Stage2", "Stage1_vs_Stage3", "Stage2_vs_Stage3"),
  FALSE_FALSE = c(7206, 6989, 6400),  # non-DGE & non-DTU
  FALSE_TRUE = c(628, 759, 1701),     # non-DTU & DGE
  TRUE_FALSE = c(1446, 1465, 1001),   # non-DGE & DTU
  TRUE_TRUE = c(121, 188, 299)        # DTU & DGE
)

# Reshape data into long format
data_dtu_long <- DTU_bar %>%
  pivot_longer(cols = -Comparison, names_to = "Category", values_to = "Count") %>%
  mutate(DGE = ifelse(Category %in% c("FALSE_TRUE", "TRUE_TRUE") , "DGE", "Non-DGE"),
         DTU = ifelse(Category %in% c("TRUE_FALSE", "TRUE_TRUE"), "DTU", "Non-DTU"))

# Plot stacked bar chart
stacked_DTU_path <- file.path(plots_dir, "stacked_DTU.pdf")
pdf(stacked_DTU_path, width = 10, height = 6)
ggplot(data_dtu_long, aes(x = DGE, y = Count, fill = DTU)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Comparison) +
  labs(title = "Stacked Bar Plot of DGE and DTU Overlaps",
       x = "DGE Status", y = "Count") +
  theme_minimal()
dev.off()


