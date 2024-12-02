
##### isoformFeatures 

input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs"

isoformFeatures <- read_tsv(file.path(input_data_dir, "isoformFeatures.tsv"))


DTE_table <- read_tsv(file.path(input_data_dir, "DTE_table.tsv"))

DGE_table <- read_tsv(file.path(input_data_dir, "DGE_table.tsv"))

colnames(isoformFeatures)

isoformFeatures$gene_id <- gsub("\\..*", "", isoformFeatures$gene_id)
isoformFeatures$isoform_id <- gsub("\\..*", "", isoformFeatures$isoform_id)

colnames(DTE_table)

DTE_table$gene_id <- gsub("\\..*", "", DTE_table$gene_id)
DTE_table$isoform_id <- gsub("\\..*", "", DTE_table$isoform_id)

DTE_table <- DTE_table |> dplyr::select(c( "isoform_id", "logFC", "logCPM", "F", "PValue", "FDR",  "condition_1", "condition_2"))
colnames(DTE_table) <- c("isoform_id", "DTE_log2FC", "DTE_logCPM", "DTE_F", "DTE_pval", "DTE_qval", "condition_1", "condition_2")


colnames(DGE_table)

DGE_table$gene_id <- gsub("\\..*", "", DGE_table$gene_id)

DGE_table <- DGE_table |> dplyr::select(c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR", "condition_1", "condition_2"))
colnames(DGE_table) <- c("gene_id", "DGE_log2FC", "DGE_logCPM", "DGE_F", "DGE_pval", "DGE_qval", "condition_1", "condition_2")

isoformFeatures <- isoformFeatures |> distinct()
DTE_table <- DTE_table |> distinct()

new_DGE_DTE_DTU <- left_join(isoformFeatures, DTE_table, 
            by = c("isoform_id", "condition_1", "condition_2"))

new_DGE_DTE_DTU <- left_join(new_DGE_DTE_DTU, DGE_table, 
            by = c("gene_id", "condition_1", "condition_2"))


# Add DTU column based on isoform_switch_q_value and dIF
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTU = isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1)

# Add DTE column based on DTE_qval and DTE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DTE = DTE_qval < 0.05 & abs(DTE_log2FC) >= 1)

# Add DGE column based on DGE_qval and DGE_log2FC
new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> 
  mutate(DGE = DGE_qval < 0.05 & abs(DGE_log2FC) >= 1)


new_DGE_DTE_DTU |> dplyr::select( gene_id, 
                              condition_1, condition_2,DGE, DTU )

table(new_DGE_DTE_DTU$DTU, useNA = "always")
table(new_DGE_DTE_DTU$DGE, useNA = "always")
table(new_DGE_DTE_DTU$DTE, useNA = "always")



#convert NA to FALSE
new_DGE_DTE_DTU$DTU[is.na(new_DGE_DTE_DTU$DTU)] <- FALSE
new_DGE_DTE_DTU$DGE[is.na(new_DGE_DTE_DTU$DGE)] <- FALSE
new_DGE_DTE_DTU$DTE[is.na(new_DGE_DTE_DTU$DTE)] <- FALSE



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


colnames(new_DGE_DTE_DTU)

new_DGE_DTE_DTU <- new_DGE_DTE_DTU |> dplyr::select(
  c( "gene_id", "isoform_id",   "DTU", "dIF", "IF1", "IF2",
  "gene_value_1"  , "gene_value_2"  , "condition_1", "condition_2", 
  "DGE_log2FC",  "DGE_logCPM",   "DGE_pval","DGE_qval", "DGE", 
  "DTE_log2FC"  , "DTE_logCPM"  , "DTE_pval"  , "DTE_qval"  , "DTE"
))

write_tsv(new_DGE_DTE_DTU, file.path(input_data_dir, "DTU_left_join_DGE_DTE_DTU.tsv")
          

