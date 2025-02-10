####### Upset Plot ########
library(readr)
library(dplyr)
library(UpSetR)

# save venn diagram as pdf 
method = "bambu"
comparison = "ROs"
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison)
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "plots", "upset")
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

rownames(gene_overlaps) <- gene_overlaps$gene_id

table(gene_overlaps$DTU_Stage_1_vs_Stage_2,  
      gene_overlaps$DGE_Stage_1_vs_Stage_2)
table(gene_overlaps$DTU_Stage_1_vs_Stage_3,  
      gene_overlaps$DGE_Stage_2_vs_Stage_3)
table(gene_overlaps$DTU_Stage_2_vs_Stage_3,
      gene_overlaps$DGE_Stage_1_vs_Stage_3)
 

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
