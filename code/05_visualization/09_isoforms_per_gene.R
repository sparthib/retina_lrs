####### Upset Plot ########
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rtracklayer)

method <- "bambu"
comparison <- "ROs"

# 1. Load counts matrix

long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                              method, comparison, "filtered_by_counts_and_biotype", "filtered_isoform_cpm.RDS")

long_read_counts <- readRDS(long_read_counts)

# long_read_counts <- remove_zero_var_rows(long_read_counts)

rownames(long_read_counts) <- gsub("\\.\\d+$", "",
                                   rownames(long_read_counts))

#2. Load GTF file 
gtf_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads/ROs_protein_coding_annotations.gtf"
gtf <- import.gff(gtf_file)
gtf_df <- as.data.frame(gtf)
rm(gtf)

gtf_df <- gtf_df |> 
  filter(type == "transcript") |> 
  select(gene_id, transcript_id) |> distinct()
gtf_df$transcript_id <- gsub("\\.\\d+$", "", gtf_df$transcript_id)

#3. Merge counts with gene_id and transcript_id

long_read_counts <- as.data.frame(long_read_counts) |>
  tibble::rownames_to_column("transcript_id")  |>
  pivot_longer(
    cols = -transcript_id,
    names_to = "sample",
    values_to = "expression"
  )

# remove rows where expression is 0
long_read_counts <- long_read_counts |> 
  filter(expression > 0)

counts_matrix_isoforms <- left_join(long_read_counts, gtf_df, 
                                    by = c("transcript_id" = "transcript_id"))

# check for NAs
na_count <- sum(is.na(counts_matrix_isoforms$gene_id))

isoform_per_gene <- counts_matrix_isoforms |> 
  group_by(gene_id, sample) |> count()

# group by sample and get the table for number of isoforms per gene
isoform_dist_per_sample <- isoform_per_gene |>
  ungroup() |>
  count(sample, n, name = "num_genes")

sample_stage_map <- tibble::tibble(
  sample = c("EP1_BRN3B_RO", "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
             "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2"),
  stage = c("Stage_3", "Stage_1", "Stage_2", 
            "Stage_3", "Stage_2", "Stage_1", "Stage_2")
)

isoform_dist_per_sample <- isoform_dist_per_sample |>
  left_join(sample_stage_map, by = "sample")

# Boxplot of isoforms per gene per sample

plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison,"protein_coding", "plots")

write_tsv(isoform_dist_per_sample, 
            file = file.path(plots_dir, "isoform_per_gene_per_sample.tsv"))

pdf(file.path(plots_dir, "isoform_per_gene_per_sample.pdf"))

ggplot(isoform_dist_per_sample, aes(x = factor(n), y = num_genes)) +
  geom_boxplot(fill = "gray90", color = "black", outlier.shape = NA) +  # single box per `n`
  geom_jitter(aes(color = stage), width = 0.2, alpha = 0.7, size = 2) +  # points by sample colored by stage
  scale_color_manual(
    values = c(
      "Stage_1" = "orange",
      "Stage_2" = "seagreen",
      "Stage_3" = "purple"
    )
  ) +
  labs(
    x = "Number of Isoforms per Gene (n)",
    y = "Number of Genes (per Sample)",
    title = "Distribution of Gene Counts per Isoform Number Across Samples"
  ) +
  theme_minimal()

dev.off()



#4. plot barplot
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison,"protein_coding", "plots")
pdf(file.path(plots_dir, "isoform_per_gene.pdf"))
isoform_per_gene |> 
  ggplot(aes(x = n)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Number of isoforms per gene",
       x = "Number of isoforms",
       y = "Number of genes")
dev.off()








