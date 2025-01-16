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

custom_palette <- c("#000000","#E69F00" ,"#56B4E9", "#009E73" ,"#F0E442", "#0072B2",
                    "#CC79A7", "#D55E00"  , "#999999")

# plot the data
combined_junction_data_long <- melt(combined_junction_data, id.vars = "sample", variable.name = "junctions", value.name = "percentage")
#only plot junctions 1-10
combined_junction_data_long <- combined_junction_data_long[combined_junction_data_long$junctions != "junctions_0",]

combined_junction_data_long$sample <- gsub("_junction_per_read", "", combined_junction_data_long$sample)
combined_junction_data_long$junctions <- gsub("junction_", "", combined_junction_data_long$junctions)
combined_junction_data_long$junctions <- gsub("junctions_", "", combined_junction_data_long$junctions)
combined_junction_data_long$junctions <- factor(
  as.numeric(as.character(combined_junction_data_long$junctions)),
  levels = sort(unique(as.numeric(as.character(combined_junction_data_long$junctions))))
)

combined_junction_data_long$sample |> unique()
combined_junction_data_long$junctions |> unique()
combined_junction_data_long$sample <- as.character(combined_junction_data_long$sample)

combined_junction_data_long$sample |> unique()

RO_samples <- c("EP1-BRN3B-RO" ,  "EP1-WT_hRO_2" ,  "EP1-WT_ROs_D45",
                "H9-BRN3B_hRO_2" ,"H9-BRN3B-RO"   , "H9-CRX_hRO_2" ,  "H9-CRX_ROs_D45")

FT_RGC_samples <- c("H9-FT_1" ,  "H9-FT_2" ,"H9-hRGC_1" , "H9-hRGC_2")

RO_combined_junction_data_long <- combined_junction_data_long |>
  dplyr::filter(sample %in% RO_samples)

FT_RGC_combined_junction_data_long <- combined_junction_data_long |>
  dplyr::filter(sample %in% FT_RGC_samples)

df <- RO_combined_junction_data_long
file <- "/users/sparthib/retina_lrs/plots/exon_exon/RO_combined_boxplots.pdf"

pdf(file)
p <- ggplot(df, aes(x = junctions, y = percentage)) +
  geom_boxplot(outlier.shape = NA) +  # Box plot without outliers
  geom_jitter(aes(color = sample), width = 0.2, size = 1, alpha = 0.9) +  # Add sample points
  theme_minimal() +
  labs(
    title = "Exon junction data by sample", 
    x = "Number of exon-exon junctions", 
    y = "percentage unique genes expressed"
  ) +
  # theme(
  #   axis.text.x = element_blank(), # Remove x-axis text
  #   axis.ticks.x = element_blank() # Remove x-axis ticks
  # ) +
  scale_color_manual(values = custom_palette) +  # Use the custom palette
  scale_y_continuous(
    breaks = seq(0.01, 0.13, by = 0.01),  # Specify breaks from 0.01 to 0.13 by 0.01
    limits = c(0.01, 0.13),  # Set y-axis limits from 0.01 to 0.13
    labels = scales::number_format(accuracy = 0.01)  # Format labels to 2 decimal places
  )

print(p)
dev.off()

# Remove the lines that blank out the x-axis text and ticks
# theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
