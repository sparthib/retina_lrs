library(reshape2)
library(ggplot2)
# Define the directory where the CSV files are stored
data_dir <- "/users/sparthib/retina_lrs/processed_data/exon_exon"

# List all the CSV files in the directory
csv_files <- list.files(data_dir, pattern = "_junction_per_read.csv$", full.names = TRUE)

# Read each CSV file into a list
junction_data_list <- lapply(csv_files, read.csv)

# Optionally, name the list elements based on the filenames (without the path and extension)
names(junction_data_list) <- gsub(".*/|\\.csv$", "", csv_files)

# Now you can access each dataset by its name in the list
head(junction_data_list[["DG-WT-hRGC_junction_per_read"]])

# To combine all datasets into a single dataframe, you can use do.call and rbind
combined_junction_data <- do.call(rbind, junction_data_list)
colnames(combined_junction_data) <- c("sample", "junctions_0", "junction_1",
                                      "junctions_2", "junctions_3", "junctions_4",
                                      "junctions_5", "junctions_6", "junctions_7",
                                      "junctions_8", "junctions_9", "junctions_10")


combined_junction_data$sample <- rownames(combined_junction_data)
#save the combined data
write.csv(combined_junction_data, file = "/users/sparthib/retina_lrs/processed_data/exon_exon/combined_junction_data.csv", row.names = FALSE)

custom_palette <- c(
  "darkgreen", "#377EB8", "#4DAF4A", "seagreen", "#FFFF33", 
  "#A65628", "#984EA3", "#999999", "#D95F02", "pink",
  "#4575B4", "#91BFDB", "#313695", "#A500D8", "violet",
  "#D73027"
)

# plot the data
combined_junction_data_long <- melt(combined_junction_data, id.vars = "sample", variable.name = "junctions", value.name = "percentage")
#only plot junctions 1-10
combined_junction_data_long <- combined_junction_data_long[combined_junction_data_long$junctions != "junctions_0",]

pdf("/users/sparthib/retina_lrs/plots/exon_exon/combined_junction_data.pdf")
p <- ggplot(combined_junction_data_long, aes(x = junctions, y = percentage)) +
  geom_boxplot(outlier.shape = NA) +  # Box plot without outliers
  geom_jitter(aes(color = sample), width = 0.2, size = 1, alpha = 0.9) +  # Add sample points
  theme_minimal() +
  labs(title = "Exon junction data by sample", x = "Number of exon-exon junctions", y = "percentage unique genes expressed") +
  theme(axis.text.x = element_blank(), # Remove x-axis text
        axis.ticks.x = element_blank()) + # Remove x-axis ticks
  scale_color_manual(values = custom_palette) +# Use the custom palette
  scale_y_continuous(
    breaks = seq(0.01, 0.13, by = 0.01),  # Specify breaks from 0.01 to 0.13 by 0.01
    limits = c(0.01, 0.13),  # Set y-axis limits from 0.01 to 0.13
    labels = scales::number_format(accuracy = 0.01)  # Format labels to 2 decimal places
  )

print(p)
dev.off()
