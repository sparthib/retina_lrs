library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)
library(readr)
library(ggplot2)

plot_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/vcf_stats/H9_EP1/variants_per_read"

tsv_files <- list.files(plot_output_dir, pattern = "tsv$", full.names = TRUE)

# Read all TSV files into a list of data frames
df_list <- lapply(tsv_files, read_tsv)

#add sample name to each data frame
df_list <- lapply(seq_along(df_list), function(i) {
  df <- df_list[[i]]
  sample_name <- tools::file_path_sans_ext(basename(tsv_files[i]))
  df$sample <- sample_name
  return(df)
})

# Combine all data frames into one
combined_df <- do.call(rbind, df_list)

#remove the 'read_id' column
combined_df$read_id <- NULL


#plot the variants_overlapped.counts_per_read and  variants_overlapped.Freq
# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)

# Assuming your data is in a dataframe called 'df'
# If you need to read it from a file, uncomment and modify the line below:
# df <- read.csv("your_file.csv")

# Clean up sample names (remove file extensions if present)
combined_df$sample_clean <- gsub("_variants.*", "", combined_df$sample)


combined_df <- combined_df |>
  dplyr::group_by(sample_clean) |>
  dplyr::mutate(Freq_percentage = (variants_overlapped.Freq / sum(variants_overlapped.Freq)) * 100) |>
  dplyr::ungroup()

write_tsv(combined_df, file.path(plot_output_dir, "variants_frequency_percentage.tsv"), col_names = TRUE)

combined_df <- combined_df |> 
  dplyr::filter(Freq_percentage > 0.1) 


# Main plot: Faceted by sample with both line and points
pdf(file.path(plot_output_dir, "variants_frequency_faceted_by_sample.pdf"), width = 12, height = 8)
ggplot(combined_df, aes(x = variants_overlapped.counts_per_read, 
               y = Freq_percentage)) +
  geom_point(color = "darkblue", size = 1, alpha = 0.8) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_x_continuous() +
  facet_wrap(~sample_clean, scales = "free_y", ncol = 2) +
  labs(title = "Variants Frequency Percentage vs. Counts Per Read",
       subtitle = "Faceted by Sample",
       x = "Counts Per Read",
       y = "Percentage of reads") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank())

dev.off()

