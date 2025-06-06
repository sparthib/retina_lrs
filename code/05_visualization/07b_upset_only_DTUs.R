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
                                              DTU , DTE)

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
                                         DTU_Stage_2_vs_Stage_3 = ifelse(DTU == TRUE & condition == "Stage_2_vs_Stage_3", TRUE, FALSE))


gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_Stage_1_vs_Stage_2 = any(DTU_Stage_1_vs_Stage_2), 
            DTU_Stage_1_vs_Stage_3 = any(DTU_Stage_1_vs_Stage_3), 
            DTU_Stage_2_vs_Stage_3 = any(DTU_Stage_2_vs_Stage_3)) 

#### the DTU columns have NAs in them, due to no expression of gene in the condition,
#### so we need to replace them with FALSE
gene_overlaps <- gene_overlaps |> 
  mutate(
    DTU_Stage_1_vs_Stage_2 = ifelse(is.na(DTU_Stage_1_vs_Stage_2), FALSE, DTU_Stage_1_vs_Stage_2),
    DTU_Stage_1_vs_Stage_3 = ifelse(is.na(DTU_Stage_1_vs_Stage_3), FALSE, DTU_Stage_1_vs_Stage_3),
    DTU_Stage_2_vs_Stage_3 = ifelse(is.na(DTU_Stage_2_vs_Stage_3), FALSE, DTU_Stage_2_vs_Stage_3)
  )

sum(is.na(gene_overlaps))

readr::write_tsv(gene_overlaps, file.path(plots_dir, "gene_overlaps_DTU_genes.tsv"))

rownames(gene_overlaps) <- gene_overlaps$gene_id 

nrow(gene_overlaps)

# Load required libraries

# Prepare data for VennDiagram
venn_data <- list(
  DTU_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_2],
  DTU_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_3],
  DTU_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_2_vs_Stage_3]
)

upset_path <- file.path(plots_dir, "DTU_Genes_Upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()


###### DTU isoform overlap ######



isoform_overlaps = DGE_DTU_DTE |> dplyr::select( isoform_id,
                                              condition_1, condition_2,DGE, 
                                              DTU , DTE)


isoform_overlaps <- isoform_overlaps |>
  mutate(condition = case_when(
    condition_1 == "Stage_1" & condition_2 == "Stage_2" ~ "Stage_1_vs_Stage_2",
    condition_1 == "Stage_1" & condition_2 == "Stage_3" ~ "Stage_1_vs_Stage_3",
    condition_1 == "Stage_2" & condition_2 == "Stage_3" ~ "Stage_2_vs_Stage_3",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))

nrow(isoform_overlaps)

isoform_overlaps <- isoform_overlaps |> mutate(DTU_Stage_1_vs_Stage_2 = ifelse(DTU == TRUE & condition == "Stage_1_vs_Stage_2", TRUE, FALSE),
                                         DTU_Stage_1_vs_Stage_3 = ifelse(DTU == TRUE & condition == "Stage_1_vs_Stage_3", TRUE, FALSE),
                                         DTU_Stage_2_vs_Stage_3 = ifelse(DTU == TRUE & condition == "Stage_2_vs_Stage_3", TRUE, FALSE))


isoform_overlaps  <- isoform_overlaps |> distinct()
isoform_overlaps <- isoform_overlaps |> dplyr::select(-c( DGE, DTU, DTE,  condition )) 

isoform_overlaps = isoform_overlaps |> group_by(isoform_id) |> 
  summarise(DTU_Stage_1_vs_Stage_2 = any(DTU_Stage_1_vs_Stage_2), 
            DTU_Stage_1_vs_Stage_3 = any(DTU_Stage_1_vs_Stage_3), 
            DTU_Stage_2_vs_Stage_3 = any(DTU_Stage_2_vs_Stage_3)) 

isoform_overlaps <- isoform_overlaps |> 
  mutate(
    DTU_Stage_1_vs_Stage_2 = ifelse(is.na(DTU_Stage_1_vs_Stage_2), FALSE, DTU_Stage_1_vs_Stage_2),
    DTU_Stage_1_vs_Stage_3 = ifelse(is.na(DTU_Stage_1_vs_Stage_3), FALSE, DTU_Stage_1_vs_Stage_3),
    DTU_Stage_2_vs_Stage_3 = ifelse(is.na(DTU_Stage_2_vs_Stage_3), FALSE, DTU_Stage_2_vs_Stage_3)
  )


sum(is.na(isoform_overlaps))

readr::write_tsv(isoform_overlaps, file.path(plots_dir, "isoform_overlaps_DTU_genes.tsv"))

rownames(isoform_overlaps) <- isoform_overlaps$isoform_id

nrow(isoform_overlaps)


# Prepare data for VennDiagram
venn_data <- list(
  DTU_Stage_1_vs_Stage_2 = rownames(isoform_overlaps)[isoform_overlaps$DTU_Stage_1_vs_Stage_2],
  DTU_Stage_1_vs_Stage_3 = rownames(isoform_overlaps)[isoform_overlaps$DTU_Stage_1_vs_Stage_3],
  DTU_Stage_2_vs_Stage_3 = rownames(isoform_overlaps)[isoform_overlaps$DTU_Stage_2_vs_Stage_3]
)

upset_path <- file.path(plots_dir, "DTU_Isoforms_Upset.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange",
           set_size.show = FALSE)
print(p)
dev.off()