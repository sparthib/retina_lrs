####### Upset Plot ########
library(readr)
library(dplyr)
library(UpSetR)
library(tidyr)
library(ggplot2)
# save venn diagram as pdf 
method = "bambu"
comparison = "ROs"

code_dir <- Sys.getenv("retina_lrs_code")

input_data_dir <- file.path(code_dir, "processed_data/dtu/",
                            method, comparison, "protein_coding")
plots_dir <- file.path(code_dir, "processed_data/dtu/",
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

#### the DTU columns have NAs in them, due to no expression of gene in the condition,
#### so we need to replace them with FALSE
gene_overlaps <- gene_overlaps |> 
  mutate(
    DTU_Stage_1_vs_Stage_2 = ifelse(is.na(DTU_Stage_1_vs_Stage_2), FALSE, DTU_Stage_1_vs_Stage_2),
    DTU_Stage_1_vs_Stage_3 = ifelse(is.na(DTU_Stage_1_vs_Stage_3), FALSE, DTU_Stage_1_vs_Stage_3),
    DTU_Stage_2_vs_Stage_3 = ifelse(is.na(DTU_Stage_2_vs_Stage_3), FALSE, DTU_Stage_2_vs_Stage_3)
  )

sum(is.na(gene_overlaps))

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
  DTU_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_2_vs_Stage_3],
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



#### DGE based on bambu gene counts matrix ####
#### this includes genes that don't have specific isoforms ####



DGE_DTU_DTE <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))

DGE_table <- read_tsv(file.path(input_data_dir, "DGE_table.tsv"))

joined_df <- DGE_DTU_DTE |> 
  left_join(DGE_table, by = "gene_id") |> 
  dplyr::select(gene_id, isoform_id, condition_1, condition_2, DGE, DTU, DTE,
                gene_name, gene_biotype) |> 
  mutate(condition = case_when(
    condition_1 == "Stage_1" & condition_2 == "Stage_2" ~ "Stage_1_vs_Stage_2",
    condition_1 == "Stage_1" & condition_2 == "Stage_3" ~ "Stage_1_vs_Stage_3",
    condition_1 == "Stage_2" & condition_2 == "Stage_3" ~ "Stage_2_vs_Stage_3",
    TRUE ~ NA_character_
  ))

## number of genes in DGE_table 
unique_DGE_ids <- DGE_table |> pull(gene_id) |> unique() 
# 16971

genes_with_isoforms <- DGE_DTU_DTE |> pull(gene_id) |> unique() 
# 15683

genes_not_with_isoforms <- setdiff(unique_DGE_ids, genes_with_isoforms) 
# 1288 + 15683 = 16971
DGE_table <- DGE_table |>
  mutate(DGE = FDR < 0.05 & abs(logFC) >= 1)

## only keep genes in set
DGE_table <- DGE_table |> 
  filter(gene_id %in% genes_not_with_isoforms) 

DGE_1_vs_2 <- DGE_table |> filter(condition_1 =="Stage_1",
                                               condition_2 == "Stage_2", 
                                              DGE == "TRUE") |> pull(gene_id)

DGE_1_vs_3 <- DGE_table |> filter(condition_1 =="Stage_1",
                                               condition_2 == "Stage_3", 
                                              DGE == "TRUE") |> pull(gene_id)
# FALSE  TRUE 
# 792   496 

DGE_2_vs_3 <- DGE_table |> filter(condition_1 =="Stage_2",
                                               condition_2 == "Stage_3", 
                                              DGE == "TRUE") |> pull(gene_id)

# DGE
# FALSE  TRUE 
# 1177   111 

#nrow gene_overlaps
gene_overlaps <- read_tsv( file.path(plots_dir, "gene_overlaps.tsv"))

rownames(gene_overlaps) <- gene_overlaps$gene_id 


# [1] "gene_id"                "DTU_Stage_1_vs_Stage_2" "DTU_Stage_1_vs_Stage_3"
# [4] "DTU_Stage_2_vs_Stage_3" "DGE_Stage_2_vs_Stage_3" "DGE_Stage_1_vs_Stage_3"
# [7] "DGE_Stage_1_vs_Stage_2" "DTE_Stage_1_vs_Stage_2" "DTE_Stage_1_vs_Stage_3"
# [10] "DTE_Stage_2_vs_Stage_3"

# initiate a dataframe with these colnames

gene_only_gene_overlaps <- data.frame(gene_id = genes_not_with_isoforms,
                          DTU_Stage_1_vs_Stage_2 = FALSE,
                          DTU_Stage_1_vs_Stage_3 = FALSE,
                          DTU_Stage_2_vs_Stage_3 = FALSE,
                          DGE_Stage_2_vs_Stage_3 = FALSE,
                          DGE_Stage_1_vs_Stage_3 = FALSE,
                          DGE_Stage_1_vs_Stage_2 = FALSE,
                          DTE_Stage_1_vs_Stage_2 = FALSE,
                          DTE_Stage_1_vs_Stage_3 = FALSE,
                          DTE_Stage_2_vs_Stage_3 = FALSE)

gene_only_gene_overlaps <- gene_only_gene_overlaps |>
  mutate(DGE_Stage_1_vs_Stage_2 = gene_id %in% DGE_1_vs_2)

gene_only_gene_overlaps$DGE_Stage_1_vs_Stage_2 |> table()

gene_only_gene_overlaps <- gene_only_gene_overlaps |>
  mutate(DGE_Stage_1_vs_Stage_3 = gene_id %in% DGE_1_vs_3)
gene_only_gene_overlaps$DGE_Stage_1_vs_Stage_3 |> table()

gene_only_gene_overlaps <- gene_only_gene_overlaps |>
  mutate(DGE_Stage_2_vs_Stage_3 = gene_id %in% DGE_2_vs_3)
gene_only_gene_overlaps$DGE_Stage_2_vs_Stage_3 |> table()


### rbind to the gene_overlaps
gene_overlaps <- rbind(gene_overlaps, gene_only_gene_overlaps)

gene_overlaps$only_found_at_gene_level <- gene_overlaps$gene_id %in% genes_not_with_isoforms


### this is the final gene_overlaps table

readr::write_tsv(gene_overlaps, file.path(plots_dir,
                                          "gene_overlaps_including_genes_with_unidentified_isoforms.tsv"))

venn_data <- list(
  DTU_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_2],
  DTU_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_1_vs_Stage_3],
  DTU_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTU_Stage_2_vs_Stage_3],
  DGE_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_2_vs_Stage_3],
  DGE_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_1_vs_Stage_3],
  DGE_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DGE_Stage_1_vs_Stage_2],
  DTE_Stage_1_vs_Stage_2 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_1_vs_Stage_2],
  DTE_Stage_1_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_1_vs_Stage_3],
  DTE_Stage_2_vs_Stage_3 = rownames(gene_overlaps)[gene_overlaps$DTE_Stage_2_vs_Stage_3]
)

upset_path <- file.path(plots_dir, "Genes_upset_including_genes_with_unidentified_isoforms.pdf")
pdf(upset_path, width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
           main.bar.color = "steelblue",
           matrix.color = "darkorange")
print(p)
dev.off()

