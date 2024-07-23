# save venn diagram as pdf 
library("ggVennDiagram")

gene_overlaps = DGE_DTU_DTE |> dplyr::select( gene_id, condition_1, condition_2,DGE, DTU)

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

gene_overlaps <- gene_overlaps |> dplyr::select(-c(gene_id, DGE, DTU, condition )) 


# Prepare data for VennDiagram
venn_data <- list(
  DTU_RO_D100_vs_RO_D45 = which(gene_overlaps$DTU_RO_D100_vs_RO_D45),
  DTU_RO_D200_vs_RO_D45 = which(gene_overlaps$DTU_RO_D200_vs_RO_D45),
  DTU_RO_D100_vs_RO_D200 = which(gene_overlaps$DTU_RO_D100_vs_RO_D200),
  DGE_RO_D100_vs_RO_D200 = which(gene_overlaps$DGE_RO_D100_vs_RO_D200),
  DGE_RO_D200_vs_RO_D45 = which(gene_overlaps$DGE_RO_D200_vs_RO_D45),
  DGE_RO_D100_vs_RO_D45 = which(gene_overlaps$DGE_RO_D100_vs_RO_D45)
)
data_unique <- venn_data %>%
  distinct()

# Convert logical values to numeric
data_binary <- data_unique %>%
  select(-gene_id) %>%
  mutate(across(everything(), ~ as.integer(.))) %>%
  as.data.frame()

# Add gene_id as row names
rownames(data_binary) <- data_unique$gene_id

# Check the transformed data
head(data_binary)


library(UpSetR)

upset(binary_data_df[,2:7], 
      main.bar.color = "steelblue", 
      matrix.color = "darkorange", 
      keep.order = TRUE)


pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE_venn.pdf")
fig <- ggVennDiagram(list(DTU_RO_D100_vs_RO_D45 = which(gene_overlaps$DTU_RO_D100_vs_RO_D45) ,
)) + 
  scale_fill_gradient(low="grey",high = "red")
print(fig)
dev.off()



# Plot the Venn diagram
ggVennDiagram(venn_data) + 
  scale_fill_gradient(low = "blue", high = "red") +
  theme(legend.position = "bottom")


# # Create the 3-way contingency table
# contingency_table <- xtabs(~ DGE + DTU + condition, data = gene_overlaps)
# 
# # Convert the contingency table to a data frame
# contingency_df <- as.data.frame(contingency_table)
# 
# # Pivot the data frame longer
# contingency_long <- contingency_df %>%
#   pivot_longer(cols = -c(DGE, DTU, condition), names_to = "variable", values_to = "value")
# 
# # Display the pivoted data
# contingency_long


gene_overlaps = gene_overlaps |> group_by(gene_id) |> 
  summarise( DGE=any(DGE), DTU=any(DTU), condition)  


pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/RO_D100_vs_RO_D45/DGE_DTU_DTE_venn.pdf")
fig <- ggVennDiagram(list(DTU = which(gene_overlaps$DTU), 
                          DGE = which(gene_overlaps$DGE),
                          DTE = which(gene_overlaps$DTE), 
)) + 
  scale_fill_gradient(low="grey",high = "red")
print(fig)
dev.off()





