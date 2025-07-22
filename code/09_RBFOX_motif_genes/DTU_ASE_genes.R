

library(Rsubread)

bam_file <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf/H9-hRGC_1_sorted_haplotagged.bam"
sample_name <- basename(bam_file)
sample_name <- gsub("\\.bam$", "", sample_name)

gtf_file <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
gene_counts <- featureCounts(files = bam_file,annot.ext= gtf_file,
              isGTFAnnotationFile = TRUE, GTF.featureType = "gene",
              isLongRead = TRUE)

gene_counts <- gene_counts$counts
colnames(gene_counts) <- sample_name

gene_counts <- as.data.frame(gene_counts) |>
  dplyr::filter(.data[[sample_name]] > 0)

gene_counts$gene_id <- rownames(gene_counts)
colnames(gene_counts) <- c("count", "gene_id")

gene_counts$gene_id <- gsub("\\..*","",gene_counts$gene_id)

require(biomaRt)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c( "external_gene_name", "ensembl_gene_id",
                "gene_biotype"),
  filter="ensembl_gene_id",
  values=gene_counts$gene_id,
  uniqueRows=TRUE)


# merge gene counts with gene annotations
gene_counts <- merge(gene_counts, annotLookup, by.x = "gene_id",
                     by.y = "ensembl_gene_id")

# order gene counts by count in descending order
gene_counts <- gene_counts[order(gene_counts$count, decreasing = TRUE), ]
readr::write_csv(gene_counts, file.path( "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf/",
                 paste0(sample_name, ".csv")))

###### DTU specific genes ######
gene_counts_all <- readr::read_csv(paste0("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf/", 
                                      sample_name,".csv"))
gene_counts_hp1 <- readr::read_csv(paste0("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf/", 
                                      "H9-hRGC_1_hp1_gene_counts.csv"))
gene_counts_hp2 <- readr::read_csv(paste0("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/vcf/",
                                      "H9-hRGC_1_hp2_gene_counts.csv"))


gene_counts_all_in_hp1 <- gene_counts_all |> dplyr::filter(gene_id %in% gene_counts_hp1$gene_id)
gene_counts_all_in_hp2 <- gene_counts_all |> dplyr::filter(gene_id %in% gene_counts_hp2$gene_id)

gene_counts_all_in_hp1 <- gene_counts_all_in_hp1 |> 
  dplyr::filter(external_gene_name == "JUND")

method <- "bambu"
comparison <- "FT_vs_RGC"
input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "protein_coding" )

FT_vs_RGC_DTU <- readr::read_tsv(file.path(input_data_dir,
                                "DGE_DTE_DTU.tsv"))

FT_vs_RGC_DTU <- FT_vs_RGC_DTU |> dplyr::filter(DTU == TRUE)

DTU_ASE_genes <- intersect(gene_counts$gene_id,
                       FT_vs_RGC_DTU$gene_id)
DTU_ASE_genes <- unique(DTU_ASE_genes)

length(DTU_ASE_genes)

DTU_ASE_genes <- gene_counts |> dplyr::filter(gene_id %in% DTU_ASE_genes)


write.csv(DTU_ASE_genes, file.path(input_data_dir,
                                "DTU_ASE_genes_HP2.csv"), row.names = FALSE)


#### 
DTU_ASE_genes <- readr::read_csv(file.path(input_data_dir,
                                "DTU_ASE_genes_HP2.csv"))




