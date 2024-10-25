library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
gtf_file <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/OUT/OUT.transcript_models.gtf"


# Import the GTF file as a GRanges object
gtf_gr <- import(gtf_file)

# Convert to a data frame if necessary
gtf_df <- as.data.frame(gtf_gr)

# View the first few rows
head(gtf_df)

# Count the number of unique genes
num_genes <- length(unique(gtf_df$gene_id))


#group the data by gene_id and count the number of isoforms per gene
n_iso_per_gene <- gtf_df |> filter(type == "transcript") |> group_by(gene_id) |> summarise(n_isoforms = n())
#number of genes in sample =  35173

nrow(n_iso_per_gene)


sum(n_iso_per_gene$n_isoforms)

# Calculate the average number of isoforms per gene
avg_iso_per_gene <- mean(n_iso_per_gene$n_isoforms)

#get range 
range(n_iso_per_gene$n_isoforms)

#get the gene id for max number of isoforms
max_iso_gene <- n_iso_per_gene[which.max(n_iso_per_gene$n_isoforms),]
# ENSG00000234741.10  

#get number of single-isoform genes
n_single_iso_genes <- sum(n_iso_per_gene$n_isoforms == 1)
# 18938


#### COLNAMES gtf_df ####


n_exon_per_gene <- gtf_df |> filter(type == "exon") |> group_by(gene_id) |> summarise(n_exons= n())

# Calculate the average number of exons per gene
avg_exons_per_gene <- mean(n_exon_per_gene$n_exons)
# 16.58957

range(n_exon_per_gene$n_exons)
# 492

# merge the two dataframes 
n_iso_exon_per_gene <- merge(n_iso_per_gene, n_exon_per_gene, by = "gene_id")

# Create dot plot using ggplot2
p <- ggplot(n_iso_exon_per_gene , aes(x = n_isoforms, y = n_exons)) +
  geom_point(color = "blue", size = 1) +
  labs(title = "Dot Plot: Number of Isoforms vs. Number of Exons",
       x = "Number of Isoforms",
       y = "Number of Exons") +
  theme_minimal()

# Save the plot as PDF
pdf("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/iso_vs_exon.pdf")
print(p)
dev.off()

# get the gene with max width 

transcripts_df <- gtf_df |> filter(type == "transcript") 

#get max width transcript from transcripts_df

transcripts_df |> filter(width == max(width))
# 
# transcripts_df |> filter(width == max(width))
# seqnames   start      end   width strand source       type score phase
# 1     chr9 8314246 10613002 2298757      - HAVANA transcript    NA    NA
# gene_id transcripts      gene_type gene_name   hgnc_id
# 1 ENSG00000153707.19        <NA> protein_coding     PTPRD HGNC:9668
# havana_gene     transcript_id similar_reference_id alternatives
# 1 OTTHUMG00000021005.6 ENST00000381196.9                 <NA>         <NA>
#   Canonical exons exon exon_id   tag transcript_type transcript_name
# 1      True    46 <NA>    <NA> basic  protein_coding       PTPRD-203
# transcript_support_level  ont    havana_transcript        protein_id
# 1                        5 <NA> OTTHUMT00000055395.3 ENSP00000370593.3



