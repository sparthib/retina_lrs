library(biomaRt)
library(ggplot2)

us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)  

method <- "bambu"
comparison <- "ROs"

counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"
counts_dir <- file.path(counts_dir, method, comparison)

RO_isoforms <- readRDS(file.path(counts_dir, "filtered", "isoform_counts.RDS"))


# get number of rownames that start with Bambu
RO_num_bambu_isoforms <- sum(grepl("^Bambu", rownames(isoforms)))

# only keep rows that have more than 0 counts across the row
RO_isoforms <- RO_isoforms[rowSums(RO_isoforms) > 0, ]

nrow(RO_isoforms)  # Number of isoforms with counts
# 171647
# get biotypes for these isoforms 

rownames(RO_isoforms) <- gsub("\\..*", "", rownames(RO_isoforms))  # Remove Ensembl version numbers

RO_isoform_annotLookup <- getBM(
  mart = mart,
  attributes = c("ensembl_transcript_id", "external_gene_name",
                 "gene_biotype", "transcript_biotype"),
  filter = "ensembl_transcript_id",
  values =  rownames(RO_isoforms) ,
  uniqueRows = TRUE
)

# > nrow(RO_isoform_annotLookup)  # Number of unique isoforms with annotations
# [1] 170804
# > nrow(RO_isoforms)
# [1] 171647

# Use table() and convert to data frame
RO_isoform_table <- as.data.frame(table(RO_isoform_annotLookup$gene_biotype))
colnames(RO_isoform_table) <- c("biotype", "count")
## add NA row for biotypes that are not present
RO_isoform_table <- rbind(RO_isoform_table, 
                          data.frame(biotype = "NA", count = 843))

RO_isoform_table$biotype <- as.character(RO_isoform_table$biotype)

# Then group rare biotypes
RO_isoform_table <- RO_isoform_table |>
  mutate(prop = count / sum(count),
         biotype_grouped = ifelse(prop < 0.02, "Other", biotype))


# Aggregate rare biotypes (<2% of total) into "Other"
RO_isoform_table$prop <- RO_isoform_table$count / sum(RO_isoform_table$count)
RO_isoform_table$biotype_grouped <- ifelse(RO_isoform_table$prop < 0.02, 
                                           "Other", RO_isoform_table$biotype)

# Summarize counts per group
grouped_table <- RO_isoform_table %>%
  group_by(biotype_grouped) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100,
         percent_label = paste0(round(percentage, 1), "%"),
         ypos = cumsum(count) - count / 2)

# Make biotype_grouped a factor with levels ordered by count (for legend order)
grouped_table$biotype_grouped <- factor(grouped_table$biotype_grouped,
                                        levels = grouped_table$biotype_grouped)

grouped_table$label <- paste0(grouped_table$biotype_grouped, "\n", grouped_table$percent_label)

pdf("/users/sparthib/retina_lrs/plots/biotypes/RO_biotype_distribution_grouped_percent.pdf",
    width = 8, height = 8)

pie(grouped_table$percentage, 
    labels = grouped_table$label,
    main = "Biotype Distribution (Grouped <2%)",
    col = RColorBrewer::brewer.pal(n = nrow(grouped_table), name = "Set3"),
    cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2
)

dev.off()



comparison <- "FT_vs_RGC"
counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"
counts_dir <- file.path(counts_dir, method, comparison)

FT_vs_RGC_isoforms <- readRDS(file.path(counts_dir, "filtered", "isoform_counts.RDS"))
# get number of rownames that start with Bambu
FT_vs_RGC_num_bambu_isoforms <- sum(grepl("^Bambu", rownames(FT_vs_RGC_isoforms)))
# only keep rows that have more than 0 counts across the row
FT_vs_RGC_isoforms <- FT_vs_RGC_isoforms[rowSums(FT_vs_RGC_isoforms) > 0, ]


rownames(FT_vs_RGC_isoforms) <- gsub("\\..*", "",
                                     rownames(FT_vs_RGC_isoforms))  # Remove Ensembl version numbers

FT_vs_RGC_isoform_annotLookup <- getBM(
  mart = mart,
  attributes = c("ensembl_transcript_id", "external_gene_name",
                 "gene_biotype", "transcript_biotype"),
  filter = "ensembl_transcript_id",
  values =  rownames(FT_vs_RGC_isoforms) ,
  uniqueRows = TRUE
)

# > nrow(FT_vs_RGC_isoform_annotLookup)  # Number of unique isoforms with annotations
# [1] 127958
# > nrow(FT_vs_RGC_isoforms)
# [1] 128690


FT_vs_RGC_isoform_table <- as.data.frame(table(FT_vs_RGC_isoform_annotLookup$gene_biotype))
colnames(FT_vs_RGC_isoform_table) <- c("biotype", "count")
## add NA FT_vs_RGCw for biotypes that are not present
FT_vs_RGC_isoform_table <- rbind(FT_vs_RGC_isoform_table, 
                          data.frame(biotype = "NA", count = 732))

FT_vs_RGC_isoform_table$biotype <- as.character(FT_vs_RGC_isoform_table$biotype)

# Then gFT_vs_RGCup rare biotypes
FT_vs_RGC_isoform_table <- FT_vs_RGC_isoform_table |>
  mutate(prop = count / sum(count),
         biotype_grouped = ifelse(prop < 0.02, "Other", biotype))


# Aggregate rare biotypes (<2% of total) into "Other"
FT_vs_RGC_isoform_table$prop <- FT_vs_RGC_isoform_table$count / sum(FT_vs_RGC_isoform_table$count)
FT_vs_RGC_isoform_table$biotype_grouped <- ifelse(FT_vs_RGC_isoform_table$prop < 0.02, 
                                           "Other", FT_vs_RGC_isoform_table$biotype)

# Summarize counts per group
grouped_table <- FT_vs_RGC_isoform_table %>%
  group_by(biotype_grouped) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100,
         percent_label = paste0(round(percentage, 1), "%"),
         ypos = cumsum(count) - count / 2)

# Make biotype_grouped a factor with levels ordered by count (for legend order)
grouped_table$biotype_grouped <- factor(grouped_table$biotype_grouped,
                                        levels = grouped_table$biotype_grouped)

grouped_table$label <- paste0(grouped_table$biotype_grouped,
                              "\n", grouped_table$percent_label)

pdf("/users/sparthib/retina_lrs/plots/biotypes/FT_vs_RGC_biotype_distribution_grouped_percent.pdf",
    width = 8, height = 8)

pie(grouped_table$percentage, 
    labels = grouped_table$label,
    main = "Biotype Distribution (Grouped <2%)",
    col = RColorBrewer::brewer.pal(n = nrow(grouped_table), name = "Set3"),
    cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2
)

dev.off()
