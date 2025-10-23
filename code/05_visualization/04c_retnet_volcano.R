library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)
# install.packages("DGEobj.utils")
library("DGEobj.utils")
library(grid)


analysis_type <- "ROs"
quant_method <- "bambu"

retnet_df <- read_tsv( file = "/users/sparthib/retina_lrs/processed_data/dtu/retnet_disease_genes.tsv" ) 

volcano_plot <- function(data, x_col, y_col, gene_label_col = "gene_name",
                         cutoff_log_qval =  1.30103, cutoff_logFC = 1, condition_1, condition_2, table_type) {
  # Classify points as significant or not
  data <- data |>
    mutate(
      significant = case_when(
        !!sym(y_col) > cutoff_log_qval & !!sym(x_col) >= cutoff_logFC ~ "Upregulation in condition 2",
        !!sym(y_col) > cutoff_log_qval & !!sym(x_col) <= -cutoff_logFC ~ "Upregulation in condition 1",
        TRUE ~ "Not Significant"
      )
    )
  
  color_map <- list(
    "FT|RGC"         = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown",  "Upregulation in condition 1" = "skyblue"),
    "RGC|Stage_1"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "orange",  "Upregulation in condition 1" = "brown"),
    "RGC|Stage_2"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen",    "Upregulation in condition 1" = "brown"),
    "RGC|Stage_3"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "brown"),
    "Stage_1|Stage_2"= c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen",    "Upregulation in condition 1" = "orange"),
    "Stage_1|Stage_3"= c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "orange"),
    "Stage_2|Stage_3"= c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "seagreen")
  )
  
  
  # Create key
  pair_key <- paste(condition_1, condition_2, sep = "|")
  
  # Assign colors
  colors <- color_map[[pair_key]]
  
  
  # Axis labels
  x_label <- switch(table_type,
                    "DTE" = "log2 Fold Change",
                    "DGE" = "log2 Fold Change",
                    "DTU" = "dIF",
                    x_col)  # fallback
  
  y_label <- "-log10(FDR)"
  
  # Create the plot
  plot <- ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col), color = significant)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-cutoff_logFC, cutoff_logFC), linetype = "dashed", color = "pink") +
    geom_hline(yintercept = cutoff_log_qval, linetype = "dashed", color = "pink") +
    labs(
      title = paste("Volcano Plot for", table_type, ":", method, condition_1, "vs", condition_2),
      x = x_label,
      y = y_label,
      color = "Significance"
    ) +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Add gene labels for significant points
  if (table_type %in% c("DTU", "DTE")){
    top_labels <- data |>
      filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) |>
      filter(!!sym(y_col) > cutoff_log_qval) |>
      filter(abs(!!sym(x_col)) >= cutoff_logFC) |>
      mutate(logFC_sign = sign(!!sym(x_col))) |>
      group_by(logFC_sign) |>
      arrange(desc(!!sym(y_col))) |>
      slice_head(n = 20) |>
      ungroup()
  }
  
  if (table_type == "DGE"){
    top_labels <- data |> dplyr::select(!!sym(x_col), !!sym(y_col), !!sym(gene_label_col)) |>
      distinct() |>
      filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) |>
      filter(!!sym(y_col) > cutoff_log_qval) |>
      filter(abs(!!sym(x_col)) >= cutoff_logFC) |>
      mutate(logFC_sign = sign(!!sym(x_col))) |>
      filter(logFC_sign != 0) |>  # Remove zero values
      group_by(logFC_sign) |>
      arrange(desc(!!sym(y_col))) |>
      slice_head(n = 10) |>  # 10 from each group = 20 total (10 positive, 10 negative)
      ungroup() 
  }
  
  plot <- plot +
    geom_text_repel(
      data = top_labels,
      aes(label = !!sym(gene_label_col)),
      size = 2,
      color = "black",
      max.overlaps = Inf)
  plot
}


generate_volcano_plots <- function(input_base_dir, comparison, table_type, conditions) {
  
  input_data_dir <- file.path(input_base_dir, method, comparison, "protein_coding") 
  conditions <- if (comparison == "ROs") {
    list(c("Stage_1", "Stage_2"), c("Stage_1", "Stage_3"), c("Stage_2", "Stage_3"))
  } else if (comparison == "FT_vs_RGC") {
    list(c("FT", "RGC"))
  } else if (comparison == "RO_vs_RGC") { 
    list(c("RGC", "Stage_1"), c("RGC", "Stage_2"), c("RGC", "Stage_3"))
  }
  
  DGE_DTE_DTU <- file.path(input_data_dir, "DGE_DTE_DTU.tsv")
  
  analysis_table <- read_tsv(DGE_DTE_DTU)
  analysis_table <- analysis_table |> filter(gene_name %in% retnet_df$gene_name)
  
  # Create plots directory
  plots_dir <- file.path(input_data_dir, "plots", "retnet")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Generate plots for each condition pair
  for (condition_pair in conditions) {
    condition_1 <- condition_pair[[1]]
    condition_2 <- condition_pair[[2]]
    
    filtered_table <- analysis_table |> 
      filter(condition_1 == !!condition_1 & condition_2 == !!condition_2)
    
    output_file <- file.path(plots_dir, paste0(condition_1, "_vs_", condition_2, "_", table_type, "_volcano.pdf"))
    
    if (table_type == "DTE"){ 
      filtered_table$log10DTE_qval <- -log10(filtered_table$DTE_qval)
      pdf(output_file)
      
      p <- volcano_plot(filtered_table, x_col = "DTE_log2FC", y_col = "log10DTE_qval", gene_label_col = "gene_name",
                        cutoff_log_qval = 1.30103, cutoff_logFC = 1, condition_1 = condition_1, condition_2 = condition_2,
                        table_type = table_type)
      print(p)
      dev.off()
      
    } else if (table_type == "DGE") {
      filtered_table$log10DGE_qval <- -log10(filtered_table$DGE_qval)
      filtered_table <- filtered_table |> dplyr::select("gene_name", "DGE_log2FC", "log10DGE_qval") |> distinct()
      pdf(output_file)
      
      p <- volcano_plot(filtered_table, x_col = "DGE_log2FC", y_col = "log10DGE_qval", gene_label_col = "gene_name",
                        cutoff_log_qval = 1.30103, cutoff_logFC = 1, condition_1 = condition_1, condition_2 = condition_2,
                        table_type = table_type)
      print(p)
      dev.off()
    } else if (table_type == "DTU") {
      filtered_table$log10DTU_qval <- -log10(filtered_table$DTU_qval)
      pdf(output_file)
      
      p <- volcano_plot(filtered_table, x_col = "dIF", y_col = "log10DTU_qval", gene_label_col = "gene_name",
                        cutoff_log_qval = 1.30103, cutoff_logFC = 0.1, condition_1 = condition_1, condition_2 = condition_2,
                        table_type = table_type)
      print(p)
      dev.off()
    }
    
  }
}


# Generate plots for each method and comparison
methods <- c("bambu")
comparisons <- c("ROs", "FT_vs_RGC", "RO_vs_RGC")
input_base_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/"

for (method in methods) {
  for (comparison in comparisons) {
    conditions <- if (comparison == "ROs") {
      list(c("Stage_1", "Stage_2"), c("Stage_1", "Stage_3"), c("Stage_2", "Stage_3"))
    } else if( comparison == "FT_vs_RGC") { 
      list(c("FT", "RGC"))
    } else if (comparison == "RO_vs_RGC") { 
      list(c("RGC", "Stage_1"), c("RGC", "Stage_2"), c("RGC", "Stage_3"))
    }
    
    generate_volcano_plots(
      input_base_dir = input_base_dir,
      comparison = comparison,
      table_type = "DTE",
      conditions = conditions
    )
    
    generate_volcano_plots(
      input_base_dir = input_base_dir,
      comparison = comparison,
      table_type = "DTU",
      conditions = conditions
    )
    
    generate_volcano_plots(
      input_base_dir = input_base_dir,
      comparison = comparison,
      table_type = "DGE",
      conditions = conditions
    )
  }
}