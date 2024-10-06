# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define a function to extract total, primary mapped, and supplementary counts from flagstat output
extract_flagstat <- function(flagstat_lines) {
  total <- as.numeric(sub(" .*", "", grep("in total", flagstat_lines, value = TRUE)))
  primary_mapped <- as.numeric(sub(" .*", "", grep("primary mapped", flagstat_lines, value = TRUE)))
  supplementary <- as.numeric(sub(" .*", "", grep("supplementary", flagstat_lines, value = TRUE)))
  return(c(total, primary_mapped, supplementary))
}

# Read the flagstat file (which contains concatenated outputs)
flagstat_file <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/logs/all.flagstat.txt"  # Replace with actual path
flagstat_data <- readLines(flagstat_file)

# Split flagstat data by samples (each sample has a block of lines, preceded by "Sample: <sample name>")
sample_blocks <- split(flagstat_data, cumsum(grepl("^Sample:", flagstat_data)))

# Initialize an empty data frame for results
df <- data.frame(Sample = character(), Total_Reads = numeric(), Primary_Mapped = numeric(), Supplementary_Reads = numeric())

# Loop through each block (sample) and extract the data
for (block in sample_blocks) {
  # Extract the sample name from the first line (e.g., "Sample: DG-WT-hRGC_bam_flagstat.txt")
  sample_name <- sub("Sample: ", "", block[1])
  
  # Extract the relevant read counts
  counts <- extract_flagstat(block)
  
  # Append the results to the data frame
  df <- rbind(df, data.frame(Sample = sample_name, Total_Reads = counts[1], Primary_Mapped = counts[2], Supplementary_Reads = counts[3]))
}



# Calculate percentages
df <- df |> mutate(Primary_Mapped_Percent = (Primary_Mapped / Total_Reads) * 100,
         Supplementary_Reads_Percent = (Supplementary_Reads / Total_Reads) * 100)

#remove "_bam_flagstat.txt" from Sample names
df$Sample <- gsub("_bam_flagstat.txt", "", df$Sample)

# Prepare data for plotting in long format
df_long <- df |> select(Sample, Primary_Mapped_Percent, Supplementary_Reads_Percent) %>%
  pivot_longer(cols = c(Primary_Mapped_Percent, Supplementary_Reads_Percent), 
               names_to = "Category", values_to = "Percentage") %>%
  mutate(Category = recode(Category, 
                           Primary_Mapped_Percent = "Primary Mapped", 
                           Supplementary_Reads_Percent = "Supplementary"))

# Plot the bar graph with percentage annotations
output_dir <- "/users/sparthib/retina_lrs/plots/bam_qc"
pdf(file.path(output_dir, "sample_wise_primary_reads_percentage.pdf"), width = 10, height = 6)
p <- ggplot(df_long, aes(x = Sample, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = ifelse(Category == "Primary Mapped", 
                               paste0(df$Primary_Mapped), 
                               paste0(df$Supplementary_Reads))), 
            position = position_stack(vjust = 0.5), size = 2, color = "black") +
  labs(title = "Sample-wise Percentage of Total, Primary Mapped, and Supplementary Reads", 
       y = "Percentage", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
