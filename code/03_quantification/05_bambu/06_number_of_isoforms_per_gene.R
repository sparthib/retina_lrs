library(dplyr)
library(tidyr)
library(readr)
library(rtracklayer)

gtf <- "extended_annotations_fa_contigs_only.gtf"
gtf <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads",
                 gtf)
gtf_gr <- import(gtf, format = "gtf")

# Convert to a data.frame
gtf_df <- as.data.frame(gtf_gr)
head(gtf_df)


gene_and_isoforms <- gtf_df |> select(c(gene_id, transcript_id)) |> distinct()
nrow(gene_and_isoforms)

RO_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered"
RO_counts <- readRDS(file.path(RO_counts_dir, "isoform_counts.RDS"))
FT_vs_RGC_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered"
FT_vs_RGC_counts <- readRDS(file.path(FT_vs_RGC_counts_dir, "isoform_counts.RDS"))


#remove rows that are 0s 
RO_counts <- RO_counts[ rowSums(RO_counts) != 0, ]
FT_vs_RGC_counts <- FT_vs_RGC_counts[ rowSums(FT_vs_RGC_counts) != 0, ]
RO_isoforms <- rownames(RO_counts)
length(RO_isoforms)
FT_vs_RGC_isoforms <- rownames(FT_vs_RGC_counts)
length(FT_vs_RGC_isoforms)

#only keep isoforms in gene_amd_isoforms that are in RO_counts
RO_gene_and_isoforms <- gene_and_isoforms[gene_and_isoforms$transcript_id %in% RO_isoforms,]
nrow(RO_gene_and_isoforms)
FT_vs_RGC_gene_and_isoforms <- gene_and_isoforms[gene_and_isoforms$transcript_id %in% FT_vs_RGC_isoforms,]
nrow(FT_vs_RGC_gene_and_isoforms)


#count the number of isoforms per gene
RO_gene_and_isoforms <- RO_gene_and_isoforms |> group_by(gene_id) |> summarise(n_isoforms = n())
FT_vs_RGC_gene_and_isoforms <- FT_vs_RGC_gene_and_isoforms |> group_by(gene_id) |> summarise(n_isoforms = n())


nrow(RO_gene_and_isoforms)
nrow(FT_vs_RGC_gene_and_isoforms)

summary(RO_gene_and_isoforms$n_isoforms)
summary(FT_vs_RGC_gene_and_isoforms$n_isoforms)


table(RO_gene_and_isoforms$n_isoforms)
table(FT_vs_RGC_gene_and_isoforms$n_isoforms)


# 48280 genes were expressed with mean number of isoforms per gene being 3.555 in ROs.
# 171647 isoforms were expressed in ROs.
# 37998 genes were expressed with mean number of isoforms per gene being 3.387 in FT_vs_RGC.
# 128690 isoforms were expressed in FT_vs_RGC.





