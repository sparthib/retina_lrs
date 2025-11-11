library(biomaRt)
library(readr)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

common_isoforms <- file.path(data_dir,"06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)
nrow(common_isoforms)


common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                                      !grepl("^BambuGene", isoquant_gene_id))
nrow(common_isoforms)

# remove version ID from ref_gene
common_isoforms$ref_gene <- gsub("\\..*", "", common_isoforms$ref_gene)

common_isoforms <- common_isoforms[1:496,]
tail(common_isoforms)

### get gene names from biomart

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = common_isoforms$ref_gene,
                    mart = mart)

# merge gene names with common_isoforms

common_isoforms <- merge(common_isoforms, gene_names, 
                         by.x = "ref_gene", 
                         by.y = "ensembl_gene_id", all.x = TRUE)

write_tsv(common_isoforms, 
          file = file.path(code_dir, "processed_data/dtu/bambu",
                           "common_isoforms_with_gene_names.tsv"),
          col_names = TRUE)



