library(dplyr)
library(ggplot2)
library(yaml)
library(tidyr)
library(patchwork)


input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/02_MinIONQC"

plot_metrics <- function(df, file){ 
  plots <- list()
  for (i in seq_along(columns_of_interest)) {
    metric <- columns_of_interest[i]
    
    # Create the boxplot for the current metric with all samples together
    p <- ggplot(df, aes(x = factor(1), y = !!sym(metric))) + # Color by sample
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
  ggsave(file , plot = combined_plot, width = 10, height = 12)
  
}



#get all yaml files in the directory and subdirectories
yaml_files <- list.files(input_dir, pattern = "\\.yaml$", 
                         recursive = TRUE, full.names = TRUE)

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

sample_names <- c( "DG-WT-hRGCs", "EP1-BRN3B-ROs", "EP1-WT_hRO_2" ,  "EP1-WT_ROs_D45" ,"H9-BRN3B_hRO_2" ,"H9-BRN3B-ROs",
                  "H9-CRX_hRO_2" ,  "H9-CRX_ROs_D45", "H9-FT_1" ,     "H9-FT_2"   ,     "H9-hRGC_1"  ,    
                  "H9-hRGC_2" , "hRGCs", "YZ-15T_hRGC","YZ-3T_hRGC"   )

yaml_data$sample <- sample_names

RO_samples <- c("EP1-BRN3B-ROs", "EP1-WT_hRO_2" ,  "EP1-WT_ROs_D45" ,"H9-BRN3B_hRO_2" ,"H9-BRN3B-ROs",
                "H9-CRX_hRO_2" ,  "H9-CRX_ROs_D45")

FT_RGC_samples <- c("H9-FT_1" ,  "H9-FT_2" ,"H9-hRGC_1"  ,    
                    "H9-hRGC_2")

yaml_data <- yaml_data |>
  rename(
    "Mean Length of all reads" = "All reads.mean.length",
    "Median q value of all reads" = "All reads.median.q",
    "Median N50 of all reads" = "All reads.N50.length",
    "Total number of reads" = "All reads.total.reads"
  )

columns_of_interest <- c("Mean Length of all reads", "Median q value of all reads", 
                         "Median N50 of all reads", "Total number of reads")

# Pivot data to long format
yaml_long <- yaml_data |>
  pivot_longer(cols = all_of(columns_of_interest), 
               names_to = "Metric", 
               values_to = "Value")

custom_palette <- c("#000000","#E69F00" ,"#56B4E9", "#009E73" ,"#F0E442", "#0072B2",
                    "#CC79A7", "#D55E00"  , "#999999")

# Assuming `yaml_data$sample` is already set to `sample_names`
yaml_data$sample <- factor(yaml_data$sample, levels = sample_names)

RO_yaml_data <- yaml_data %>%
  filter(sample %in% RO_samples)
FT_RGC_yaml_data <- yaml_data %>%
  filter(sample %in% FT_RGC_samples)

RO_yaml_data <- RO_yaml_data |> select(c("sample","Mean Length of all reads", "Median q value of all reads", 
                         "Median N50 of all reads", "Total number of reads" ))

FT_RGC_yaml_data <- FT_RGC_yaml_data |> select(c("sample","Mean Length of all reads", "Median q value of all reads", 
                         "Median N50 of all reads", "Total number of reads" ))


write.csv(RO_yaml_data, "/users/sparthib/retina_lrs/plots/02_MinIONQC/RO_stats.csv", row.names = FALSE)
write.csv(FT_RGC_yaml_data, "/users/sparthib/retina_lrs/plots/02_MinIONQC/FT_RGC_stats.csv", row.names = FALSE)

df <- RO_yaml_data
file <- "/users/sparthib/retina_lrs/plots/02_MinIONQC/ROs_combined_boxplots.pdf"
# Create a list to hold all plots


# Loop through each metric and create a plot

plot_metrics(RO_yaml_data, file = "/users/sparthib/retina_lrs/plots/02_MinIONQC/ROs_combined_boxplots.pdf")
plot_metrics(FT_RGC_yaml_data, file = "/users/sparthib/retina_lrs/plots/02_MinIONQC/FT_RGC_combined_boxplots.pdf")


###### Q > 10 ########

yaml_data <- yaml_data |>
  rename(
    "Mean Length of reads Q >= 10" = "Q>=10.mean.length",
    "Median q value of reads Q >= 10" = "Q>=10.median.q",
    "Median N50 of reads Q >= 10" = "Q>=10.N50.length",
    "Total number of reads Q >= 10" = "Q>=10.total.reads"
  )

columns_of_interest <- c("Mean Length of reads Q >= 10", "Median q value of reads Q >= 10", 
                         "Median N50 of reads Q >= 10", "Total number of reads Q >= 10")



RO_yaml_data <- yaml_data %>%
  filter(sample %in% RO_samples)
FT_RGC_yaml_data <- yaml_data %>%
  filter(sample %in% FT_RGC_samples)

RO_yaml_data <- RO_yaml_data |> select(c("sample","Mean Length of reads Q >= 10", "Median q value of reads Q >= 10", 
                                         "Median N50 of reads Q >= 10", "Total number of reads Q >= 10" ))
FT_RGC_yaml_data <- FT_RGC_yaml_data |> select(c("sample","Mean Length of reads Q >= 10", "Median q value of reads Q >= 10", 
                                                 "Median N50 of reads Q >= 10", "Total number of reads Q >= 10"))

write.csv(RO_yaml_data, "/users/sparthib/retina_lrs/plots/02_MinIONQC/RO_stats_Q10.csv", row.names = FALSE)
write.csv(FT_RGC_yaml_data, "/users/sparthib/retina_lrs/plots/02_MinIONQC/FT_RGC_stats_Q10.csv", row.names = FALSE)


plot_metrics(RO_yaml_data, file = "/users/sparthib/retina_lrs/plots/02_MinIONQC/ROs_combined_boxplots_Q10.pdf")
plot_metrics(FT_RGC_yaml_data, file = "/users/sparthib/retina_lrs/plots/02_MinIONQC/FT_RGC_combined_boxplots_Q10.pdf")


