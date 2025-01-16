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

##get gene_id that doesn't start with ENSG 
DGE_DTU_DTE |> filter(!grepl("^ENSG", gene_id))

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
                                         # DGE_RO_D200_vs_RO_D100 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D100", TRUE, FALSE),
                                         # DGE_RO_D200_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         # DGE_RO_D100_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         # DTE_RO_D100_vs_RO_D45 = ifelse(DTE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         # DTE_RO_D200_vs_RO_D45 = ifelse(DTE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         # DTE_RO_D200_vs_RO_D100 = ifelse(DTE == TRUE & condition == "RO_D200_vs_RO_D100", TRUE, FALSE))
)

gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_RO_D100_vs_RO_D45 = any(DTU_RO_D100_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D45 = any(DTU_RO_D200_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D100 = any(DTU_RO_D200_vs_RO_D100), 
            # DGE_RO_D200_vs_RO_D100 = any(DGE_RO_D200_vs_RO_D100), 
            # DGE_RO_D200_vs_RO_D45 = any(DGE_RO_D200_vs_RO_D45), 
            # DGE_RO_D100_vs_RO_D45 = any(DGE_RO_D100_vs_RO_D45),
            # DTE_RO_D100_vs_RO_D45 = any(DTE_RO_D100_vs_RO_D45),
            # DTE_RO_D200_vs_RO_D45 = any(DTE_RO_D200_vs_RO_D45),
            # DTE_RO_D200_vs_RO_D100 = any(DTE_RO_D200_vs_RO_D100)
            ) 

rownames(gene_overlaps) <- gene_overlaps$gene_id


# Prepare data for VennDiagram
venn_data <- list(
  DTU_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D45],
  DTU_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D45],
  DTU_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D100] #,
  # DGE_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D100],
  # DGE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D45],
  # DGE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D45],
  # DTE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D100_vs_RO_D45],
  # DTE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D200_vs_RO_D45],
  # DTE_RO_D200_vs_RO_D100 = rownames(gene_overlaps)[gene_overlaps$DTE_RO_D200_vs_RO_D100]
)

upset_path <- file.path(plots_dir, "DTU_genes_upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()




# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "B_RO_D100" &
#                                      condition_2 == "C_RO_D45"), "RO_D100_vs_RO_D45",
#              plots_dir)
# 
# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
#                                      condition_2 == "C_RO_D45"), "RO_D200_vs_RO_D45",
#              plots_dir)
# 
# plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" &
#                                      condition_2 == "B_RO_D100"), "RO_D100_vs_RO_D200",
#              plots_dir)
