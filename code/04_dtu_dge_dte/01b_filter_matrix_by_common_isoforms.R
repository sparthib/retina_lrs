

common_isoforms <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)

#if ref_gene column has values starting with BambuGene or isoquant_gene_id has 
#column starting with BambuGene

common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                             !grepl("^BambuGene", isoquant_gene_id))


colnames(common_isoforms)

filter_counts <- function( method, comparison ) {
    method <- "bambu"
    comparison <- "ROs"
    counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"
    counts_dir <- file.path(counts_dir, method, comparison)
  
    for (file in list.files(counts_dir, full.names=TRUE)) {
      
      if(basename(file) == "isoform_counts.RDS" | basename(file) == "isoform_cpm.RDS") {
        counts <- readRDS(file)
        print(paste("pre filter", nrow(counts)))
        counts <- counts[
            (grepl("^Bambu", rownames(counts)) & rownames(counts) %in% common_isoforms$ref_id) |
                (grepl("^transcript", rownames(counts)) & rownames(counts) %in% common_isoforms$isoquant_isoform_id) |
                grepl("^ENST", rownames(counts)),
        ]
        print(paste("post filter", nrow(counts)))
        output_file <- dir.create(file.path(counts_dir, "filtered"), showWarnings=FALSE, recursive=TRUE)
        output_file <- file.path(counts_dir, "filtered", basename(file))
        saveRDS(counts, output_file)
      }
    }
    
    for (file in list.files(counts_dir, full.names=TRUE)) {
      
      if(basename(file) == "gene_counts.RDS" | basename(file) == "gene_cpm.RDS") {
        counts <- readRDS(file)
        counts <- counts[grepl("^ENSG", rownames(counts)),]
        output_file <- file.path(counts_dir, "filtered", basename(file))
        saveRDS(counts, output_file)
      }
    }
}

filter_counts("bambu", "ROs")
filter_counts("bambu", "FT_vs_RGC")
filter_counts("Isoquant", "ROs")
filter_counts("Isoquant", "FT_vs_RGC")

### Verify that the filtering worked
# filtered_counts_dir <- file.path(counts_dir, "Isoquant", "ROs", "filtered")
# isoform_counts <- read.table(file.path(filtered_counts_dir, "isoform_counts.RDS"), header=TRUE, sep="\t")
# gene_counts <- read.table(file.path(filtered_counts_dir, "gene_counts.RDS"), header=TRUE, sep="\t")
# 
# nrow(isoform_counts)
# #check if there are any rownames with BambuGene or transcript
# sum(grepl("^Bambu", rownames(isoform_counts)))
# sum(grepl("^transcript", rownames(isoform_counts)))
# sum(grepl("^ENST", rownames(isoform_counts)))
# 
# 
# nrow(gene_counts)
# sum(grepl("^Bambu", rownames(gene_counts))
