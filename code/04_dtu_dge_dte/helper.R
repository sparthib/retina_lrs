
# Function to filter rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, function(x) min(x) != max(x)), ]
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


plot_barplot <- function(df, condition, output_plots_dir){ 
  
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


###### GO Analysis ######

# Function to load and filter DGE data
load_dge_data <- function(file_path, cond1, cond2) {
  read_tsv(file_path) %>%
    filter(condition_1 == cond1, condition_2 == cond2)
}


get_gene_list <- function(dge, fdr = 0.05, log2fc_cutoff = 1) {
  filtered <- dge %>%
    select(DGE_log2FC, DGE_qval, gene_id) %>%
    distinct() %>%
    filter(DGE_qval < fdr, abs(DGE_log2FC) >= log2fc_cutoff)
  
  gene_list <- setNames(
    sort(filtered$DGE_log2FC, decreasing = TRUE),
    gsub("\\..*", "", filtered$gene_id)
  )
  gene_list
}

# Function to perform over-representation analysis and save results
run_ora <- function(genelist, ont, output_data_dir) {
  ego <- enrichGO(
    gene = names(genelist),
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ont,
    pAdjustMethod = "fdr",
    minGSSize = 100,
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.01,
    readable = TRUE
  )
  
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }
  
  write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("DGE_ora_", ont, ".tsv")))
}


# Function to generate plots for GO analysis
run_ora_plots <- function(genelist, ont, output_plot_dir) {
  ego <- enrichGO(
    gene = names(genelist),
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ont,
    pAdjustMethod = "fdr",
    minGSSize = 100,
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.01,
    readable = TRUE
  )
  
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  pdf(file.path(output_plot_dir, paste0("DGE_ora_dotplot_", ont, ".pdf")))
  print(dotplot(ego, showCategory = 15))
  dev.off()
}




