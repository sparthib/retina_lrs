library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

isoquant_ROs
isoquant_FT_RGC

isoquant_FT_RGC <- read_tsv("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/categories.tsv",
                            col_names = FALSE)

colnames(isoquant_FT_RGC) <- c("isoform", "chr", "strand", "X4",
                               "X5", "category", "gene_id", "isoform_id",
                               "X9", "X10", "subcategory")
                          


unique(isoquant_FT_RGC$category)

output_plots_dir <- "/users/sparthib/retina_lrs/plots/sqanti/isoquant/"
#make a stacked bar plot for all the categories and print to pdf 
#add percentage to the plot

pdf(paste0(output_plots_dir, "FT_RGC_categories.pdf"))
p <- isoquant_FT_RGC |>
  group_by(category) |>
  summarise(count = n()) |>
  mutate(percentage = count / sum(count) * 100) |>
  ggplot(aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity") +
  coord_polar("y") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "FT_RGC",
       x = "",
       y = "Proportion of category in observed novel isoforms") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5))
print(p)
dev.off()

isoquant_ROs <- read_tsv("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/categories.tsv")

colnames(isoquant_ROs) <- c("isoform", "chr", "strand", "X4",
                               "X5", "category", "gene_id", "isoform_id",
                               "X9", "X10", "subcategory")

unique(isoquant_ROs$category)

output_plots_dir <- "/users/sparthib/retina_lrs/plots/sqanti/isoquant/"
#make a stacked bar plot for all the categories and print to pdf
#add percentage to the plot

pdf(paste0(output_plots_dir, "ROs_categories.pdf"))
p <- isoquant_ROs |>
  group_by(category) |>
  summarise(count = n()) |>
  mutate(percentage = count / sum(count) * 100) |>
  ggplot(aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity") +
  coord_polar("y") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "ROs",
       x = "",
       y = "Proportion of category in observed novel isoforms") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5))
print(p)
dev.off()


categories <- paste(as.vector(isoquant_FT_RGC$subcategory), collapse = ";" )
categories <- strsplit(categories, ";")
categories <- unlist(categories)
categories <- unique(categories)
#remove "." from categories


for (value in categories) {
  isoquant_FT_RGC[[value]] <- ifelse(str_detect(isoquant_FT_RGC$subcategory, value), 1, 0)
}

ft_categories <- tibble(
  category = character(),
  percentage = numeric()
)
for (value in categories){ 
  percentage <- sum(isoquant_FT_RGC[[value]]) * 100 / nrow(isoquant_FT_RGC)
  ft_categories  <- ft_categories  |> add_row(category = value, percentage = percentage)
}

#make a bar plot for each of the categories with y axis as percentage
output_plots_dir <- "/users/sparthib/retina_lrs/plots/sqanti/isoquant/"
pdf(paste0(output_plots_dir, "FT_RGC_subcategories.pdf"))
p <- ft_categories |>
  ggplot(aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightblue") +  # Set bar color to light blue
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Adjust axis tick label size and rotation
        axis.text.y = element_text(size = 8),  # Adjust y-axis tick label size
        axis.title = element_text(size = 10)) +  # Adjust axis title size
  labs(title = "FT_RGC",
       x = "Subcategories",
       y = "Proportion of subcategory in observed novel isoforms") +
  geom_text(data = subset(ft_categories, percentage >= 15),
            aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2) +  # Adjust horizontal justification for better positioning
  theme(axis.line = element_line(colour = "black")) 
print(p)
dev.off()

categories <- paste(as.vector(isoquant_ROs$subcategory), collapse = ";" )
categories <- strsplit(categories, ";")
categories <- unlist(categories)
categories <- unique(categories)
#remove "." from categories

for (value in categories) {
  isoquant_ROs[[value]] <- ifelse(str_detect(isoquant_ROs$subcategory, value), 1, 0)
}

ro_categories <- tibble(
  category = character(),
  percentage = numeric()
)

for (value in categories){ 
  percentage <- sum(isoquant_ROs[[value]]) * 100 / nrow(isoquant_ROs)
  ro_categories  <- ro_categories  |> add_row(category = value, percentage = percentage)
}

#make a bar plot for each of the categories with y axis as percentage
output_plots_dir <- "/users/sparthib/retina_lrs/plots/sqanti/isoquant/"
pdf(paste0(output_plots_dir, "ROs_subcategories.pdf"))
p <- ro_categories |>
  ggplot(aes(x = category, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightblue") +  # Set bar color to light blue
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Adjust axis tick label size and rotation
        axis.text.y = element_text(size = 8),  # Adjust y-axis tick label size
        axis.title = element_text(size = 10)) +  # Adjust axis title size
  labs(title = "FT_RGC",
       x = "Subcategories",
       y = "Proportion of subcategory in observed novel isoforms") +
  geom_text(data = subset(ro_categories, percentage >= 15),
            aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2) +  # Adjust horizontal justification for better positioning
  theme(axis.line = element_line(colour = "black")) 
print(p)
dev.off()


