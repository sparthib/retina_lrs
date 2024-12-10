library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("/users/sparthib/retina_lrs/code/05_visualization/helper.R")
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison) 

volcano_plot <- function(data, x_col = "logFC", y_col = "FDR", gene_label_col = "external_gene_name",
                         cutoff_qval = 0.05, cutoff_logFC = 1, condition_1, condition_2) {
  # Classify points as significant or not
  data <- data |> 
    mutate(
      significant = ifelse(
        !!sym(y_col) < cutoff_qval & abs(!!sym(x_col)) >= cutoff_logFC,
        "Significant",
        "Not Significant"
      )
    )
  
  # Define color scheme
  colors <- c("Not Significant" = "gray", "Significant" = "darkgreen")
  
  # Create the plot
  plot <- ggplot(data, aes(x = !!sym(x_col), y = -log10(!!sym(y_col)), color = significant)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-cutoff_logFC, cutoff_logFC), linetype = "dashed", color = "pink") +
    geom_hline(yintercept = -log10(cutoff_qval), linetype = "dashed", color = "pink") +
    labs(
      title = paste("Volcano Plot for", table_type, ":", method, condition_1, "vs", condition_2),
      x = "Log2 Fold Change",
      y = "-log10 FDR",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Add gene labels for significant points
  plot <- plot +
    geom_text_repel(
      data = data |> 
        arrange(!!sym(y_col)) |> 
        slice_head(n = 20),
      aes(label = !!sym(gene_label_col)), 
      size = 2,
      color = "black",        # Set annotation text color to black
      max.overlaps = Inf
    )
  plot
}

generate_volcano_plots <- function(input_data_dir, comparison, table_type, conditions) {
  # Load analysis table
  table_file <- file.path(input_data_dir, paste0(table_type, "_table.tsv"))
  analysis_table <- read_tsv(table_file)
  
  # Process identifiers based on table type
  if (table_type == "DTE") {
    analysis_table <- analysis_table |> 
      mutate(isoform_id = ifelse(
        grepl("^ENST", isoform_id),
        gsub("\\..*", "", isoform_id),
        isoform_id
      ))
  }
  
  if (table_type == "DGE") {
    analysis_table <- analysis_table |> 
      mutate(gene_id = gsub("\\..*", "", gene_id))
  }
  
  # Map IDs to gene names using biomaRt
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  id_column <- if (table_type == "DTE") "isoform_id" else "gene_id"
  filters <- if (table_type == "DTE") "ensembl_transcript_id" else "ensembl_gene_id"
  attributes <- if (table_type == "DTE") c("ensembl_transcript_id", "external_gene_name") else c("ensembl_gene_id", "external_gene_name")
  
  results <- getBM(
    attributes = attributes,
    filters = filters,
    values = analysis_table[[id_column]],
    mart = ensembl
  ) |> distinct()
  
  colnames(results) <- c(id_column, "external_gene_name")
  analysis_table <- merge(analysis_table, results, by = id_column, all.x = TRUE) |>
    mutate(external_gene_name = replace_na(external_gene_name, "unknown"))
  
  # Create plots directory
  plots_dir <- file.path(input_data_dir, "plots", "volcano")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Generate plots for each condition pair
  for (condition_pair in conditions) {
    condition_1 <- condition_pair[[1]]
    condition_2 <- condition_pair[[2]]
    
    filtered_table <- analysis_table |> 
      filter(condition_1 == !!condition_1 & condition_2 == !!condition_2)
    
    output_file <- file.path(plots_dir, paste0(condition_1, "_vs_", condition_2, "_", table_type, "_volcano.pdf"))
    pdf(output_file)
    p <- volcano_plot(filtered_table, x_col = "logFC", y_col = "FDR", gene_label_col = "external_gene_name",
                      cutoff_qval = 0.05, cutoff_logFC = 1, condition_1 = condition_1, condition_2 = condition_2)
    print(p)
    dev.off()
  }
}

# Generate plots for each method and comparison
methods <- c("bambu", "Isoquant")
comparisons <- c("ROs", "FT_vs_RGC")
input_base_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/"

for (method in methods) {
  for (comparison in comparisons) {
    conditions <- if (comparison == "ROs") {
      list(c("B_RO_D100", "C_RO_D45"), c("A_RO_D200", "C_RO_D45"), c("A_RO_D200", "B_RO_D100"))
    } else {
      list(c("FT", "RGC"))
    }
    
    generate_volcano_plots(
      input_data_dir = file.path(input_base_dir, method, comparison),
      comparison = comparison,
      table_type = "DTE",
      conditions = conditions
    )
    
    generate_volcano_plots(
      input_data_dir = file.path(input_base_dir, method, comparison),
      comparison = comparison,
      table_type = "DGE",
      conditions = conditions
    )
  }
}


