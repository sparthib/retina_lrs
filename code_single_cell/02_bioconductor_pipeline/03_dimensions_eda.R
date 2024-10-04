
library(ggplot2)
library(dplyr)
library(tidyr)

# Input data
data <- data.frame(
  Sample = c("10x_D100-EP1_A1", "10x_D100-EP1_A2", "10x_D100-EP1_B1", "10x_D100-EP1_B2",
             "10x_D200-EP1-1_A1", "10x_D200-EP1-1_A2", "10x_D200-EP1-1_B1", "10x_D200-EP1-1_B2",
             "10x_D200-EP1-2_A1", "10x_D200-EP1-2_A2", "10x_D200-EP1-2_B1", "10x_D200-EP1-2_B2"),
  Genome = c(18117032, 699749, 32457291, 3874547, 31997316, 2519136, 35374512, 2136278, 
             19896090, 539154, 20918783, 654685),
  Transcriptome = c(9797100, 407313, 16332639, 2247286, 10797588, 1328475, 22593735, 
                    1106423, 6927127, 289590, 7094870, 349144)
)

# Convert data to long format for easier plotting
data_long <- data %>%
  pivot_longer(cols = c("Genome", "Transcriptome"), names_to = "Category", values_to = "Reads")

# Calculate percentages for each sample
data_long <- data_long %>%
  group_by(Category) %>%
  mutate(Percentage = Reads / sum(Reads) * 100)

# Plot the data using ggplot2
ggplot(data_long, aes(x = Sample, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Primary Mapped Reads: Genome vs Transcriptome",
       x = "Sample",
       y = "Number of Reads") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Category, scales = "free_y")



##### 

data <- data.frame(
  Sample = c("10x_D100_EP1_A1", "10x_D100_EP1_A2", "10x_D100_EP1_B1", "10x_D100_EP1_B2",
             "10x_D200_EP1_1_A1", "10x_D200_EP1_1_A2", "10x_D200_EP1_1_B1", "10x_D200_EP1_1_B2",
             "10x_D200_EP1_2_A1", "10x_D200_EP1_2_A2", "10x_D200_EP1_2_B1", "10x_D200_EP1_2_B2"),
  isoforms = c(19098, 2426, 25329, 8344, 19975, 5486, 25686, 4840, 15095, 1693, 15328, 1959),
  cells= c(2390, 2287, 2419, 2343, 987, 902, 992, 903, 802, 614, 809, 619)
)

# Convert data to long format for easier plotting
data_long <- data |>
  pivot_longer(cols = c("isoforms", "cells"), names_to = "Category", values_to = "Value")

# Plot the data using ggplot2, faceted by Category
ggplot(data_long, aes(x = Sample, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Value), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5) +  # Annotate with values
  labs(title = "Number of isoforms and cells",
       x = "Sample",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Category, scales = "free_y")