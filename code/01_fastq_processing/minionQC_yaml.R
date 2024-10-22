library(dplyr)
library(ggplot2)
library(yaml)
library(tidyr)
library(patchwork)


input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/02_MinIONQC"

#get all yaml files in the directory and subdirectories
yaml_files <- list.files(input_dir, pattern = "\\.yaml$", recursive = TRUE, full.names = TRUE)

#read yaml files
yaml_data <- lapply(yaml_files, yaml::yaml.load_file)

#unlist all nested lists
yaml_data <- lapply(yaml_data, function(x) unlist(x, recursive = TRUE))

#convert list to data frame
yaml_data <- do.call(rbind, yaml_data)
yaml_data <- as_tibble(yaml_data)

#convert to numberic values except for the first and last columns
yaml_data <- yaml_data %>%
  mutate(across(-c(1, ncol(yaml_data)), as.numeric))

sample_names <- c( "DG-WT-hRGCs"  , "EP1-BRN3B-ROs", "EP1-WT_hRO_2" ,  "EP1-WT_ROs_D45" ,"H9-BRN3B_hRO_2" ,"H9-BRN3B-ROs",
                  "H9-CRX_hRO_2" ,  "H9-CRX_ROs_D45", "H9-FT_1" ,     "H9-FT_2"   ,     "H9-hRGC_1"  ,    
                  "H9-hRGC_2" , "hRGCs", "YZ-15T_hRGC","YZ-3T_hRGC"   )

yaml_data$sample <- sample_names


columns_of_interest <- c("All reads.mean.length", "All reads.median.q", "All reads.N50.length", "All reads.total.reads")
# Pivot data to long format
yaml_long <- yaml_data %>%
  pivot_longer(cols = all_of(columns_of_interest), 
               names_to = "Metric", 
               values_to = "Value")


custom_palette <- c(
  "darkgreen", "#377EB8", "#4DAF4A", "seagreen", "#FFFF33", 
  "#A65628", "#984EA3", "#999999", "#D95F02", "pink",
  "#4575B4", "#91BFDB", "#313695", "#A500D8", "violet",
  "#D73027"
)

# Assuming `yaml_data$sample` is already set to `sample_names`
yaml_data$sample <- factor(yaml_data$sample, levels = sample_names)

# Columns to plot
columns_of_interest <- c("All reads.mean.length", "All reads.median.q", "All reads.N50.length", "All reads.total.reads")

# Create a list to hold all plots
plots <- list()

# Loop through each metric and create a plot
for (i in seq_along(columns_of_interest)) {
  metric <- columns_of_interest[i]
  
  # Create the boxplot for the current metric with all samples together
  p <- ggplot(yaml_data, aes(x = factor(1), y = !!sym(metric))) + # Color by sample
    geom_boxplot(outlier.shape = NA) + # Boxplot without outliers
    geom_jitter(aes(color = sample), width = 0.2, size = 2, alpha = 0.9) + # Jittered points for each sample
    labs(title = metric, x = "", y = "Value") + # Title and y-axis label
    theme(axis.text.x = element_blank(), # Remove x-axis text
          axis.ticks.x = element_blank()) + # Remove x-axis ticks
    scale_color_manual(values = custom_palette) # Use the custom palette
  
  # Remove the legend for all but the first plot
  if (i < 4) {
    p <- p + theme(legend.position = "none") # Remove legend
  }
  
  # Add the plot to the list
  plots[[metric]] <- p
}

# Combine all plots into one
combined_plot <- wrap_plots(plots, ncol = 2) # You can change `ncol` to arrange in multiple columns

# Save the combined plot to a PDF
ggsave("/users/sparthib/retina_lrs/plots/02_MinIONQC/combined_boxplots.pdf", plot = combined_plot, width = 10, height = 12)



