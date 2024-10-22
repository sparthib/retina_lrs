library(tidyr)
library(dplyr)
# get all quant files from /dcs04/hicks/data/sparthib/06_quantification/oarfish

file.list <- list.files("/dcs04/hicks/data/sparthib/06_quantification/oarfish",
                        pattern = ".quant", full.names = TRUE)

# read in all quant files
all.files<- lapply(file.list, read.table, header = TRUE, sep = "\t")

# get the sample names
sample.names <- gsub(".quant", "", basename(file.list))

# add sample names to the dataframes
for(i in 1:length(all.files)){
  all.files[[i]]$sample <- sample.names[i]
}

# combine all quant files into one dataframe
all.files <- do.call(rbind, all.files)

all.files <- all.files |> select(sample, tname, num_reads)

#fill NA  num_reads with 0
all.files$num_reads[is.na(all.files$num_reads)] <- 0

nrow(all.files)
filtered_data <- all.files |> 
  group_by(tname) |> 
  filter(var(num_reads) != 0) |> 
  ungroup()

nrow(filtered_data)
wide_data <- filtered_data %>%
  pivot_wider(names_from = sample, values_from = num_reads, values_fill = 0)

# Remove the tname column for correlation calculation
wide_data_matrix <- wide_data %>%
  select(-tname) %>%
  as.matrix()

# Calculate pairwise correlation matrix
correlation_matrix <- cor(wide_data_matrix)

correlation_long <- melt(correlation_matrix)

# Define the desired sample order
sample_order <- c( "EP1-BRN3B-RO", "EP1-WT_ROs_D45", "H9-BRN3B-RO", 
                  "H9-CRX_ROs_D45", "H9-FT_1", "H9-FT_2", "H9-hRGC_1", 
                  "H9-hRGC_2", "hRGC", "DG-WT-hRGC", "YZ-15T_hRGC", "YZ-3T_hRGC")

# Convert Var1 and Var2 to factors with the specified levels
correlation_long$Var1 <- factor(correlation_long$Var1, levels = sample_order)
correlation_long$Var2 <- factor(correlation_long$Var2, levels = sample_order)

# View the correlation matrix
pdf("/users/sparthib/retina_lrs/plots/oarfish_correlation.pdf")
p <- ggplot(correlation_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "red", 
                      limits = c(0.85, 1), 
                      name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise Correlation Matrix Heatmap (0.85 to 1)", x = "Samples", y = "Samples")
print(p)
dev.off()
# Plot the correlation matrix



