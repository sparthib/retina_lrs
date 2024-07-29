library(biomaRt)
library(readr)
library(here)


# Load the data
splicing_factors <- read_csv(here("raw_data", "GeneCards-Pathway-Splicing.csv"))


###### ADD ENSEMBL ID TO THE DATA AND SAVE ######
# head(splicing_factors)
# 
# #convert gene symbol to ENSEMBL ID
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_ids <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
#                filters = "hgnc_symbol", values = splicing_factors$`Gene Symbol`,
#                mart = ensembl)
# 
# 
# # Merge the data
# splicing_factors <- merge(splicing_factors, gene_ids, by.x = "Gene Symbol", by.y = "hgnc_symbol")
# 
# # Save the data
# write_csv(splicing_factors, here("raw_data", "GeneCards-Pathway-Splicing.csv"))


##### ADD GENE COUNTS #####

# Load FT vs RGC gene counts matrix 

bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation"
counts <- read.table(file.path(bambu_dir, "counts_gene.txt"),
                     header = TRUE)


#remove "_primary_over_30_chr_only_sorted" in column names 
colnames(counts) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(counts))
colnames(counts)[1] <- "gene_id"
head(counts)

group <- factor(c("FT", "FT", "RGC", "RGC"))




