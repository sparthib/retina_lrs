library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)

volcano_plot <- function(data, x_col, y_col, gene_label_col = "gene_name",
                         cutoff_log_qval = 1.30103, cutoff_logFC = 1, condition_1, condition_2, 
                         table_type, is_allele_comparison = FALSE) {
  
  # Determine if this is an allele comparison (H1 vs H2)
  if (missing(is_allele_comparison)) {
    # Only H1 vs H2 comparisons should be treated as allele comparisons
    is_allele_comparison <- (grepl("H[12]_Stage[123].*vs.*H[12]_Stage[123]", paste(condition_1, condition_2, sep = "_vs_")) |
                               grepl("H[12]_FT.*vs.*H[12]_FT", paste(condition_1, condition_2, sep = "_vs_")) |
                               grepl("H[12]_RGC.*vs.*H[12]_RGC", paste(condition_1, condition_2, sep = "_vs_"))) &
      grepl("H1.*vs.*H2|H2.*vs.*H1", paste(condition_1, condition_2, sep = "_vs_"))
  }
  
  # Classify points as significant or not
  if (is_allele_comparison) {
    # Simple classification for allele comparisons - just significant or not
    data <- data |>
      mutate(
        significant = case_when(
          !!sym(y_col) > cutoff_log_qval & abs(!!sym(x_col)) >= cutoff_logFC ~ "Significant",
          TRUE ~ "Not Significant"
        )
      )
    
    # Simple color scheme for allele comparisons
    colors <- c("Not Significant" = "gray", "Significant" = "red")
    
  } else {
    # Original classification for non-allele comparisons
    data <- data |>
      mutate(
        significant = case_when(
          !!sym(y_col) > cutoff_log_qval & !!sym(x_col) >= cutoff_logFC ~ "Upregulation in condition 2",
          !!sym(y_col) > cutoff_log_qval & !!sym(x_col) <= -cutoff_logFC ~ "Upregulation in condition 1",
          TRUE ~ "Not Significant"
        )
      )
    
    # Original color mapping for different comparisons
    color_map <- list(
      "FT|RGC"         = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown",  "Upregulation in condition 1" = "skyblue"),
      "RGC|Stage_1"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "orange",  "Upregulation in condition 1" = "brown"),
      "RGC|Stage_2"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen",    "Upregulation in condition 1" = "brown"),
      "RGC|Stage_3"    = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "brown"),
      "Stage_1|Stage_2"= c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen",    "Upregulation in condition 1" = "orange"),
      "Stage_1|Stage_3"= c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "orange"),
      "Stage_2|Stage_3"= c("Not Significant" = "gray", "Upregulation in condition 2" = "purple",    "Upregulation in condition 1" = "seagreen"),
      "Stage1|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "orange"),
      "Stage2|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "seagreen"),
      "Stage3|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "purple"),
      # Add within-haplotype stage comparisons (these should have directional colors)
      "H1_Stage1|H1_Stage2" = c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen", "Upregulation in condition 1" = "orange"),
      "H1_Stage2|H1_Stage3" = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple", "Upregulation in condition 1" = "seagreen"),
      "H1_Stage1|H1_Stage3" = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple", "Upregulation in condition 1" = "orange"),
      "H2_Stage1|H2_Stage2" = c("Not Significant" = "gray", "Upregulation in condition 2" = "seagreen", "Upregulation in condition 1" = "orange"),
      "H2_Stage2|H2_Stage3" = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple", "Upregulation in condition 1" = "seagreen"),
      "H2_Stage1|H2_Stage3" = c("Not Significant" = "gray", "Upregulation in condition 2" = "purple", "Upregulation in condition 1" = "orange"),
      "H1_Stage1|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "orange"),
      "H1_Stage2|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "seagreen"),
      "H1_Stage3|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "purple"),
      # Add FT vs RGC comparisons within haplotypes
      "H1_FT|H1_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "skyblue"),
      "H2_FT|H2_RGC" = c("Not Significant" = "gray", "Upregulation in condition 2" = "brown", "Upregulation in condition 1" = "skyblue")
    )
    
    # Create key for color mapping
    pair_key <- paste(condition_1, condition_2, sep = "|")
    colors <- color_map[[pair_key]]
    
    # Fallback to default colors if not found
    if (is.null(colors)) {
      colors <- c("Not Significant" = "gray", "Upregulation in condition 2" = "blue", "Upregulation in condition 1" = "red")
    }
  }
  
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
      title = paste("Volcano Plot for", table_type, ":", paste(condition_1, "vs", condition_2)),
      x = x_label,
      y = y_label,
      color = "Significance"
    ) +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Add gene labels for significant points
  if (is_allele_comparison) {
    # For allele comparisons, label top 20 by absolute fold change and significance
    top_labels <- data |>
      filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) |>
      filter(!!sym(y_col) >= cutoff_log_qval & abs(!!sym(x_col)) >= cutoff_logFC) |>
      arrange(desc(abs(!!sym(x_col))), !!sym(y_col)) |>
      slice_head(n = 20)
    
  } else {
    # Original labeling logic for non-allele comparisons
    if (table_type %in% c("DTU", "DTE")){
      top_labels <- data |>
        filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) |>
        filter(!!sym(y_col) >= cutoff_log_qval) |>
        mutate(logFC_sign = sign(!!sym(x_col))) |>
        group_by(logFC_sign) |>
        arrange(desc(!!sym(y_col))) |>
        slice_head(n = 20) |>
        ungroup()
    }
    
    if (table_type == "DGE"){
      top_labels <- data |> 
        dplyr::select(!!sym(x_col), !!sym(y_col), !!sym(gene_label_col)) |>
        distinct() |>
        filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col))) |>
        filter(!!sym(y_col) >= cutoff_log_qval) |>
        mutate(logFC_sign = sign(!!sym(x_col))) |>
        filter(logFC_sign != 0) |>
        group_by(logFC_sign) |>
        arrange(desc(!!sym(y_col))) |>
        slice_head(n = 10) |>
        ungroup() 
    }
  }
  
  plot <- plot +
    geom_text_repel(
      data = top_labels,
      aes(label = !!sym(gene_label_col)),
      size = 2,
      color = "black",
      max.overlaps = Inf
    )
  
  return(plot)
}

# Function to generate volcano plots for allele-specific analysis
generate_allele_volcano_plots <- function(input_dir, contrast_names, table_type = "DGE") {
  
  # Create plots directory
  plots_dir <- file.path(input_dir, "plots", "volcano")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  
  # Generate plots for each contrast
  for (contrast_name in contrast_names) {
    # Read the data file
    input_file <- file.path(input_dir, paste0(contrast_name, "_DGEs.tsv"))
    
    if (!file.exists(input_file)) {
      warning(paste("File not found:", input_file))
      next
    }
    
    data <- read_tsv(input_file)
    
    # Extract condition names from contrast name
    conditions <- strsplit(contrast_name, "_vs_")[[1]]
    condition_1 <- conditions[1]
    condition_2 <- conditions[2]
    
    # Determine if this is an allele comparison
    is_allele_comp <- grepl("H[12]_Stage[123].*vs.*H[12]_Stage[123]", contrast_name) &
      grepl("H1.*vs.*H2|H2.*vs.*H1", contrast_name) |
      grepl("H[12]_FT.*vs.*H[12]_FT", contrast_name) &
      grepl("H1.*vs.*H2|H2.*vs.*H1", contrast_name) |
      grepl("H[12]_RGC.*vs.*H[12]_RGC", contrast_name) &
      grepl("H1.*vs.*H2|H2.*vs.*H1", contrast_name)
    
    # Set appropriate cutoffs
    cutoff_logFC <- if (table_type == "DTU") 0.1 else 1
    
    # Create the plot
    output_file <- file.path(plots_dir, paste0(contrast_name, "_", table_type, "_volcano.pdf"))
    
    pdf(output_file)
    
    p <- volcano_plot(
      data = data,
      x_col = "logFC", 
      y_col = "neg_log10_FDR",
      gene_label_col = "gene_name",
      cutoff_log_qval = 1.30103,  # -log10(0.05)
      cutoff_logFC = cutoff_logFC,
      condition_1 = condition_1,
      condition_2 = condition_2,
      table_type = table_type,
      is_allele_comparison = is_allele_comp
    )
    
    print(p)
    dev.off()
    
    cat("Generated plot:", output_file, "\n")
  }
}

# FT vs RGC contrast names


# Generate volcano plots for ROs
# generate_allele_volcano_plots(
#   input_dir = ros_dge_output_dir,
#   contrast_names = ros_contrast_names,
#   table_type = "DGE"
# )

# Generate volcano plots for FT vs RGC
# generate_allele_volcano_plots(
#   input_dir = "/users/sparthib/retina_lrs/processed_data/ASE/DGE/FT_vs_RGC/filtered_counts",
#   contrast_names = ft_rgc_contrast_names,
#   table_type = "DGE"
# )
