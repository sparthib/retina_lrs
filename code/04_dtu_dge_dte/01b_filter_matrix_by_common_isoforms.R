

common_isoforms <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu",
                             "bambu_isoquant_refmap.txt")
common_isoforms <- read.table(common_isoforms, header=TRUE, sep="\t")
head(common_isoforms)

#if ref_gene column has values starting with BambuGene or isoquant_gene_id has 
#column starting with BambuGene

common_isoforms <- common_isoforms |> dplyr::filter(!grepl("^BambuGene", ref_gene) & 
                                             !grepl("^BambuGene", isoquant_gene_id))
nrow(common_isoforms)


colnames(common_isoforms)

counts_dir <-"/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"

bambu_RO_genes <- readRDS(file.path(counts_dir, "bambu", "ROs", "gene_counts_ROs.RDS"))
bambu_FT_vs_RGC_genes <- readRDS(file.path(counts_dir, "bambu", "FT_vs_RGC", "gene_counts_FT_vs_RGC.RDS"))

bambu_RO_isoforms <- readRDS(file.path(counts_dir, "bambu", "ROs", "isoform_counts_ROs.RDS"))
bambu_FT_vs_RGC_isoforms <- readRDS(file.path(counts_dir, "bambu", "FT_vs_RGC", "isoform_counts_FT_vs_RGC.RDS"))


filter_counts <- function( method, comparison ) {
    counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices"
    counts_dir <- file.path(counts_dir, method, comparison)
  
    for (file in list.files(counts_dir, full.names=TRUE)) {
      
      if (grepl("isoform", basename(file))) {
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
      
      if (grepl("gene", basename(file))) {
        counts <- readRDS(file)
        counts <- counts[grepl("^ENSG", rownames(counts)),]
        output_file <- file.path(counts_dir, "filtered", basename(file))
        saveRDS(counts, output_file)
      }
    }
}

filter_counts("bambu", "ROs")
filter_counts("bambu", "FT_vs_RGC")



 