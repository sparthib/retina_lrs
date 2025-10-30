
library(readr)
library(dplyr)
library(ggplot2)

# Initialize an empty data frame to collect all data
chrom_counts <- data.frame()


comparisons <- c("H1_Stage1_vs_H2_Stage1", "H1_Stage2_vs_H2_Stage2",
                 "H1_Stage3_vs_H2_Stage3","H1_Stage1_vs_H1_Stage2",
                 "H1_Stage2_vs_H1_Stage3", "H1_Stage1_vs_H1_Stage3",
                 "H2_Stage1_vs_H2_Stage2", "H2_Stage2_vs_H2_Stage3",
                 "H2_Stage1_vs_H2_Stage3")

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/ROs/"

allele_comparisons <- c("H1_Stage1_vs_H2_Stage1", "H1_Stage2_vs_H2_Stage2",
                        "H1_Stage3_vs_H2_Stage3")

for (i in seq_len(length(allele_comparisons))) {
  read_file <- paste0(dge_output_dir, allele_comparisons[i], "_DGEs.tsv")
  tt <- read_tsv(read_file, show_col_types = FALSE)
  tt <- tt |> filter(significant == TRUE)
  
  # Count chromosome occurrences
  chrom_table <- table(tt$chromosome_name)
  
  # Convert to data frame and add comparison name
  chrom_df <- as.data.frame(chrom_table)
  colnames(chrom_df) <- c("chromosome", "count")
  chrom_df$comparison <- allele_comparisons[i]
  
  # Add to full data
  chrom_counts <- bind_rows(chrom_counts, chrom_df)
}

pdf(file.path(dge_output_dir,"chromosome_counts.pdf"), width = 10, height = 6)
ggplot(chrom_counts, aes(x = chromosome, y = count, fill = comparison)) +
  geom_col(position = "dodge") +
  labs(
    title = "Significant DE Genes per Chromosome",
    x = "Chromosome",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()