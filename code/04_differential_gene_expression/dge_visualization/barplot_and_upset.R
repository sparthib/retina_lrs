# save venn diagram as pdf 
library("ggVennDiagram")


DGE_DTU_DTE <- readr::read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE.tsv")

gene_overlaps = DGE_DTU_DTE |> dplyr::select( gene_id, condition_1, condition_2,DGE, DTU )

gene_overlaps <- gene_overlaps %>%
  mutate(condition = case_when(
    condition_1 == "RO_D100" & condition_2 == "RO_D45" ~ "RO_D100_vs_RO_D45",
    condition_1 == "RO_D200" & condition_2 == "RO_D45" ~ "RO_D200_vs_RO_D45",
    condition_1 == "RO_D100" & condition_2 == "RO_D200" ~ "RO_D100_vs_RO_D200",
    TRUE ~ NA_character_  # This line handles any other cases that don't match the above
  ))

#create new columns 


gene_overlaps <- gene_overlaps |> mutate(DTU_RO_D100_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D200_vs_RO_D45 = ifelse(DTU == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DTU_RO_D100_vs_RO_D200 = ifelse(DTU == TRUE & condition == "RO_D100_vs_RO_D200", TRUE, FALSE),
                                         DGE_RO_D100_vs_RO_D200 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D200", TRUE, FALSE),
                                         DGE_RO_D200_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D200_vs_RO_D45", TRUE, FALSE),
                                         DGE_RO_D100_vs_RO_D45 = ifelse(DGE == TRUE & condition == "RO_D100_vs_RO_D45", TRUE, FALSE))
                                      
gene_overlaps  <- gene_overlaps |> distinct()
gene_overlaps <- gene_overlaps |> dplyr::select(-c( DGE, DTU,  condition )) 

gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise(DTU_RO_D100_vs_RO_D45 = any(DTU_RO_D100_vs_RO_D45), 
            DTU_RO_D200_vs_RO_D45 = any(DTU_RO_D200_vs_RO_D45), 
            DTU_RO_D100_vs_RO_D200 = any(DTU_RO_D100_vs_RO_D200), 
            DGE_RO_D100_vs_RO_D200 = any(DGE_RO_D100_vs_RO_D200), 
            DGE_RO_D200_vs_RO_D45 = any(DGE_RO_D200_vs_RO_D45), 
            DGE_RO_D100_vs_RO_D45 = any(DGE_RO_D100_vs_RO_D45)) 

rownames(gene_overlaps) <- gene_overlaps$gene_id


# Prepare data for VennDiagram
venn_data <- list(
  DTU_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D45],
  DTU_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D200_vs_RO_D45],
  DTU_RO_D100_vs_RO_D200 = rownames(gene_overlaps)[gene_overlaps$DTU_RO_D100_vs_RO_D200],
  DGE_RO_D100_vs_RO_D200 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D200],
  DGE_RO_D200_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D200_vs_RO_D45],
  DGE_RO_D100_vs_RO_D45 = rownames(gene_overlaps)[gene_overlaps$DGE_RO_D100_vs_RO_D45]
)


library("UpSetR")
pdf("./processed_data/dtu/DTU_gandall/bambu/ROs/plots/DGE_DTU_ROs_upset.pdf",
    width = 10, height = 6)
p <- upset(fromList(venn_data), nsets = 15,order.by = "freq", 
      main.bar.color = "steelblue",
      matrix.color = "darkorange")
print(p)
dev.off()



############### Bar plots, DTU, DGE, DTE ###############

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
  ggsave(path = "./processed_data/dtu/DTU_gandall/bambu/ROs/plots/",
         device = "pdf", plot = p, filename = paste0(condition, "df_barplot.pdf"))
  
}

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "RO_D100" &
                                     condition_2 == "RO_D45"), "RO_D100_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "RO_D200" &
                                     condition_2 == "RO_D45"), "RO_D200_vs_RO_D45")

plot_barplot(DGE_DTU_DTE |> filter(condition_1 == "RO_D100" &
                                     condition_2 == "RO_D200"), "RO_D100_vs_RO_D200")



