library(readr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
# install.packages("DGEobj.utils")
library("DGEobj.utils")
library(grid)


###### ADD ENSEMBL ID TO THE DATA AND SAVE ######
# head(splicing_factors)
# 
# #convert gene symbol to ENSEMBL ID
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ids <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                filters = "hgnc_symbol", values = splicing_factors$`Gene Symbol`,
#                mart = ensembl)
# 
# 
# # Merge the data
# splicing_factors <- merge(splicing_factors, gene_ids, by.x = "Gene Symbol", by.y = "hgnc_symbol")
# 
# # Save the data
# write_csv(splicing_factors, here("raw_data", "GeneCards-Pathway-Splicing.csv"))


# Load required helper functions
source("/users/sparthib/retina_lrs/code/04_dtu_dge_dte/bambu/helper.R")
# Handle rows with zero variance
remove_zero_var_rows <- function(mat) {
  mat[apply(mat, 1, var) != 0, , drop = FALSE]
}

load_gene_counts_matrix <- function(analysis_type, quant_method, counts_matrix_dir, splicing_factors_path) {
  # Validate inputs
  if (!analysis_type %in% c("FT_vs_RGC", "ROs")) {
    stop("Invalid analysis_type. Choose 'FT_vs_RGC' or 'ROs'.")
  }
  if (!quant_method %in% c("bambu", "isoquant")) {
    stop("Invalid quant_method. Choose 'bambu' or 'isoquant'.")
  }
  
  # Load splicing factors
  splicing_factors <- read_csv(splicing_factors_path) |>
    select(ensembl_gene_id, `Gene Symbol`) |>
    rename(gene_id = ensembl_gene_id, gene_name = `Gene Symbol`) |>
    distinct(gene_id, .keep_all = TRUE)
  
  # Define file paths based on inputs
  isoform_file <- file.path(counts_matrix_dir, quant_method, analysis_type, "isoform_cpm.RDS")
  gene_file <- file.path(counts_matrix_dir, quant_method, analysis_type, "gene_cpm.RDS")
  
  # Load isoform and gene TPM matrices
  isoform_tpm <- readRDS(isoform_file)
  gene_tpm <- readRDS(gene_file)
  
  # Process gene TPM
  rownames(gene_tpm) <- gsub("\\..*", "", rownames(gene_tpm)) # Remove version numbers
  gene_tpm <- gene_tpm[rownames(gene_tpm) %in% splicing_factors$gene_id, ] # Filter to splicing factors
  

  gene_tpm <- remove_zero_var_rows(gene_tpm)
  
  # Define groups based on analysis_type
  groups <- if (analysis_type == "ROs") {
    c("RO_D45", "RO_D45", "RO_D100", "RO_D100", "RO_D100", "RO_D200", "RO_D200")
  } else {
    c("FT", "FT", "RGC", "RGC")
  }
  
  # Define output directory based on analysis_type and quant_method
  output_plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu", 
                                quant_method, 
                                analysis_type, 
                                "plots", 
                                "splicing_factor_analysis")
  
  # Create the output directory if it doesn't exist
  dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(list(isoform_tpm = isoform_tpm, gene_tpm = gene_tpm, splicing_factors = splicing_factors, groups = groups, output_plots_dir = output_plots_dir))
}





# Function to plot heatmap
plot_heatmap <- function(tpm, groups, compare, output_plots_dir, splicing_factors_df) {
  # Prepare TPM data for significant genes
  tpm <- tpm |>
    as.data.frame() |>
    rownames_to_column(var = "gene_id") |>
    inner_join(splicing_factors_df, by = "gene_id") |>
    mutate(gene = paste(gene_name, gene_id, sep = "_")) |>
    column_to_rownames(var = "gene") |>
    select(-gene_id, -gene_name)
  
  write.tsv(rownames(tpm), file.path(output_plots_dir, paste0(compare, "_significant_genes.tsv")))
  # Convert to matrix and scale
  tpm_matrix <- tpm |>
    as.matrix() |>
    t() |>
    scale(center = TRUE, scale = TRUE) |>
    t()
  
  # Define heatmap color function
  col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red"))
  
  # Define annotation based on comparison
  ha <- switch(
    compare,
    "FT_vs_RGC" = HeatmapAnnotation(
      type = groups,
      annotation_name_side = "left",
      col = list(type = c("FT" = "lightgreen", "RGC" = "brown")),
      annotation_name_gp = gpar(fontsize = 0.75)
    ),
    "ROs" = HeatmapAnnotation(
      type = groups,
      annotation_name_side = "left",
      col = list(type = c("RO_D200" = "purple", "RO_D45" = "orange", "RO_D100" = "seagreen")),
      annotation_name_gp = gpar(fontsize = 0.75)
    )
  )
  
  # Create and save heatmap as PDF
  pdf(file.path(output_plots_dir, paste0(compare, "_splicing_factors_heatmap.pdf")))
  ht_list <- Heatmap(
    tpm_matrix,
    name = "Scaled TPM Expression of Splicing Factors",
    row_km = 5,
    col = col_fun,
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_title = "Isoforms",
    row_names_gp = gpar(fontsize = 1),
    column_names_gp = gpar(fontsize = 5),
    show_row_dend = TRUE,
    show_column_dend = TRUE
  )
  draw(ht_list)
  dev.off()
}

# Plot heatmap
groups <- c("RO_D45", "RO_D45", "RO_D100", "RO_D100", "RO_D100", "RO_D200", "RO_D200")
plot_heatmap(gene_tpm, groups, "ROs", output_plots_dir, splicing_factors_df)



counts_matrix_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/"
splicing_factors_path <- "/users/sparthib/retina_lrs/raw_data/GeneCards-Pathway-Splicing.csv"

# Load data for ROs using bambu
data <- load_gene_counts_matrix(
  analysis_type = "ROs", 
  quant_method = "bambu", 
  counts_matrix_dir = counts_matrix_dir, 
  splicing_factors_path = splicing_factors_path
)

# Access isoform and gene TPM matrices, groups, and output directory
isoform_tpm <- data$isoform_tpm
gene_tpm <- data$gene_tpm
groups <- data$groups
output_plots_dir <- data$output_plots_dir

# Output the dynamic output_plots_dir
print(output_plots_dir)


# Plot heatmap wrapper
plot_all_heatmaps <- function() {
  # Define combinations
  methods <- c("bambu", "isoquant")
  comparisons <- c("ROs") # Add "FT_vs_RGC" if needed
  
  # Iterate through all combinations
  for (method in methods) {
    for (comparison in comparisons) {
      # Load data
      data <- load_gene_counts_matrix(comparison, method, counts_matrix_dir, splicing_factors_path)
      gene_tpm <- data$gene_tpm
      groups <- data$groups
      output_dir <- data$output_plots_dir
  
      
      # Run plot_heatmap
      tryCatch({
        plot_heatmap(
          dge_dir = NULL, # Not used in this example
          quant_name = method,
          compare = comparison,
          tpm = gene_tpm,
          groups = groups,
          output_plots_dir = output_dir
        )
        message(sprintf("Successfully plotted heatmap for %s - %s", method, comparison))
      }, error = function(e) {
        message(sprintf("Failed to plot heatmap for %s - %s: %s", method, comparison, e$message))
      })
    }
  }
}

plot_all_heatmaps()



