library(tidyr)
library(dplyr)
library(ggplot2)

##### load 

transcript_gc_df <- readRDS("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/transcript_meta_info.rds")
covg_info_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/GAlignments"
rds_file <- "all_samples_transcript_info.rds"

covg_info <- readRDS(file.path(covg_info_dir, rds_file))


covg_info$isoform |> head()


isoform_ids <- covg_info$isoform |> unique()
#get gene names from biomart

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query the database
annotLookup <- getBM(
  mart = mart,
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = isoform_ids,
  uniqueRows = TRUE
)

colnames(annotLookup) <- c("isoform", "gene_id", "gene_name")


#group by isoform and average covg, merge gene info

isoform_info <- covg_info |> 
  group_by(isoform) |> 
  summarise(across(contains("coverage"), mean, na.rm = TRUE))

isoform_info <- merge(isoform_info, annotLookup, by = "isoform", all.x = TRUE)

colnames(transcript_gc_df) <- c("isoform", 
                                "isoform_GC_percent", "isoform_length")

isoform_info <- merge(isoform_info , transcript_gc_df, 
                      by = "isoform", all.x = TRUE)
  
  
#only keep multi-isoform genes

#get number of multi-isoform genes 
gene_counts <- isoform_info |> 
  group_by(gene_id) |> 
  summarise(n_isoforms = n())

gene_counts |> filter(n_isoforms > 1) |> nrow()
# 11366

isoform_info <- isoform_info |> 
  filter(gene_id %in% gene_counts$gene_id[gene_counts$n_isoforms > 1])


#for each gene, plot 



# Reshape data from wide to long format
coverage_long <- isoform_info |>
  pivot_longer(
    cols = starts_with("coverage_"),
    names_to = "coverage_point",
    names_prefix = "coverage_",
    values_to = "coverage_value"
  ) |>
  mutate(coverage_point = as.numeric(coverage_point))


# Add isoform_GC_percent and isoform_length to the isoform labels for annotation
coverage_long <- coverage_long %>%
  mutate(
    isoform_label = paste(
      isoform, 
      "(GC:", isoform_GC_percent, 
      "%, Length:", isoform_length, "bp)"
    )
  )


# Plot and save each gene on a separate page
pdf("/dcs04/hicks/data/sparthib/retina_lrs/08_coverage/coverage_plots/gene_wise/coverage_plot.pdf")
unique_genes <- unique(coverage_long$gene_name)

for (gene in unique_genes) {
  # Filter data for the current gene
  gene_data <- coverage_long %>% filter(gene_name == gene)
  
  # Plot for the current gene
  p <- ggplot(gene_data, aes(x = coverage_point, y = coverage_value, color = isoform_label)) +
    geom_line() +
    labs(
      x = "Coverage Point",
      y = "Coverage Value",
      title = paste("Coverage Plot for Gene:", gene)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Isoform (GC %, Length)"))
  
  # Print the plot to the PDF
  print(p)
}

# Close the PDF
dev.off()