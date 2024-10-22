# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Define paths for data loading
base_data_dir <- "/dcs04/hicks/data/sparthib"
# retina_data_path <- file.path(base_data_dir, "retina_lrs/08_coverage/GAlignments/all_samples_transcript_info.rds")
sgnex_data_path <- file.path(base_data_dir, "sg_nex_data/coverage/transcript_info/all_samples_transcript_info.rds")

# Load data (choose appropriate dataset by commenting/uncommenting)
# all_sample_transcript_info <- readRDS(retina_data_path)
all_sample_transcript_info <- readRDS(sgnex_data_path)

# Set number of bins
nbins <- 100

# Compute total and average coverage
all_sample_transcript_info <- all_sample_transcript_info %>%
  mutate(
    total_coverage = rowSums(across(starts_with("coverage"))),
    average_coverage = total_coverage / nbins,
    gc_content_bin = cut(percentage_gene_gc_content, breaks = seq(0, 100, 10))
  )

# Create grouped datasets for plotting
gc_bin_vs_average_coverage <- all_sample_transcript_info %>%
  select(gc_content_bin, average_coverage) %>%
  group_by(gc_content_bin)

transcript_length_average_coverage <- all_sample_transcript_info %>%
  select(length_bin, average_coverage) %>%
  group_by(length_bin)

transcript_biotypes_average_coverage <- all_sample_transcript_info %>%
  select(transcript_biotype, average_coverage) %>%
  group_by(transcript_biotype)

# Define output directory and create it if necessary
output_plot_dir <- "/users/sparthib/retina_lrs/plots/coverage/high_level_2D/"
dir.create(output_plot_dir, showWarnings = FALSE)

# Function to generate and save plots
# Function to generate and save plots with point count annotation
save_plot <- function(data, x, y, title, x_label, y_label, filename, width = 8, height = 6, angle = 0) {
  # Calculate the number of data points per group for annotation
  data_count <- data %>%
    group_by(!!sym(x)) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Generate the boxplot with count annotation
  p <- ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
    geom_boxplot() +
    geom_text(
      data = data_count,
      aes(x = !!sym(x), y = max(data[[y]], na.rm = TRUE), label = count),
      vjust = -0.5, size = 3
    ) +
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal()
  
  # Apply custom text angle if needed
  if (angle > 0) {
    p <- p + theme(axis.text.x = element_text(size = 7, angle = angle, hjust = 1))
  }
  
  # Save the plot
  ggsave(filename, plot = p, device = "pdf", width = width, height = height)
}


# Plot and save for GC Content Bin
save_plot(
  data = gc_bin_vs_average_coverage,
  x = "gc_content_bin",
  y = "average_coverage",
  title = "Cross Transcript Coverage vs. GC Content Bin",
  x_label = "GC Content Bin",
  y_label = "Cross Transcript Coverage",
  filename = paste0(output_plot_dir, "sgnex_cross_transcript_coverage_vs_gc_content_bin.pdf")
)

# Plot and save for Transcript Length Bin
save_plot(
  data = transcript_length_average_coverage,
  x = "length_bin",
  y = "average_coverage",
  title = "Cross Transcript Coverage vs. Transcript Length Bin",
  x_label = "Transcript Length Bin",
  y_label = "Cross Transcript Coverage",
  filename = paste0(output_plot_dir, "sgnex_cross_transcript_coverage_vs_transcript_length_bin.pdf")
)

# Plot and save for Transcript Biotype with angled x-axis labels
save_plot(
  data = transcript_biotypes_average_coverage,
  x = "transcript_biotype",
  y = "average_coverage",
  title = "Cross Transcript Coverage vs. Transcript Biotype",
  x_label = "Transcript Biotype",
  y_label = "Cross Transcript Coverage",
  filename = paste0(output_plot_dir, "sgnex_cross_transcript_coverage_vs_transcript_biotype.pdf"),
  angle = 45
)


