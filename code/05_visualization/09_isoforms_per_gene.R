####### Upset Plot ########
library(readr)
library(dplyr)
library(UpSetR)
library(tidyr)
library(ggplot2)
library(rtracklayer)

method <- "bambu"
comparison <- "ROs"

# 1. Load counts matrix
long_read_counts <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                              method, comparison, "filtered_by_counts_and_biotype", "isoform_cpm.RDS")

long_read_counts <- readRDS(long_read_counts)

# long_read_counts <- remove_zero_var_rows(long_read_counts)

rownames(long_read_counts) <- gsub("\\.\\d+$", "", rownames(long_read_counts))

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
library(dplyr)

counts_matrix_isoforms <- data.frame(isoform_id = rownames(long_read_counts)) |> 
  left_join(gtf_df, by = c("isoform_id" = "transcript_id"))

isoform_per_gene <- counts_matrix_isoforms |> 
  group_by(gene_id) |> count()

table(isoform_per_gene$n)

# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 5060 3464 2264 1476  879  500  317  195  132  101   64   44   36   19   14    4 
# 17   18   19   20   21   22   23   24   25   26   27   28   31   33   37 
# 12   11    2    5    3    3    2    3    1    1    1    1    1    1    1

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

