library(IsoformSwitchAnalyzeR)
library(dplyr)
library(rtracklayer)
library(Biostrings)

# https://github.com/GoekeLab/bambu/issues/354
# https://github.com/GoekeLab/bambu/issues/440

#### LOAD GTF #### 
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"
# Define file paths
gtf_file <- paste0(bambu_dir, "/extended_annotations.gtf")

gtf <- import(gtf_file)
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)


##### Get isoform counts info ######
method <- "bambu"
comparison <- "ROs"

matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered_by_counts_and_biotype")

full_isoform_counts <- file.path(bambu_dir, "counts_transcript.txt")
full_isoform_counts <- read.table(full_isoform_counts, header = TRUE, sep = "\t")

full_isoform_counts$TXNAME <- gsub("\\..*", "", full_isoform_counts$TXNAME)
full_isoform_counts$GENEID <- gsub("\\..*", "", full_isoform_counts$GENEID)
full_isoform_counts <- full_isoform_counts[,1:9]
colnames(full_isoform_counts) <- c("isoform_id", "gene_id", 
                              "EP1_BRN3B_RO", "EP1_WT_hRO_2", 
                              "EP1_WT_ROs_D45", "H9_BRN3B_hRO_2", 
                              "H9_BRN3B_RO", "H9_CRX_hRO_2", "H9_CRX_ROs_D45")


gene_counts <- file.path(matrix_dir, "filtered_gene_counts.RDS")
gene_counts <- readRDS(gene_counts)
rownames(gene_counts) <- gsub("\\..*", "", rownames(gene_counts))

##### Subset full isoform counts to genes that are in gene_counts ####
isoform_counts <- full_isoform_counts[full_isoform_counts$gene_id %in% 
                                         rownames(gene_counts), ]

RO_group <- c("Stage_3", "Stage_2", "Stage_1", "Stage_2",
              "Stage_3", "Stage_2","Stage_1")
library(edgeR)
dge <- DGEList(counts = isoform_counts)
keep_10 <- filterByExpr(dge, group = RO_group, min.count = 10) ##default is 10.
keep_2 <- filterByExpr(dge, group = RO_group, min.count = 2)

dge_10 <- dge[keep_10, ]
dge_2 <- dge[keep_2, ]

dge_10$genes$gene_id |> unique() |> length()
dge_2$genes$gene_id |> unique() |> length()






####
rdata_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "protein_coding", "rds", "SwitchList.rds")
SwitchListFiltered <- readRDS(rdata_path)

isoformFeatures <- SwitchListFiltered$isoformFeatures
genes_in_isoform_matrix <- isoformFeatures$gene_id |> unique()
length(genes_in_isoform_matrix)

setdiff_genes <- setdiff(rownames(gene_counts), 
                         genes_in_isoform_matrix)

length(setdiff_genes)

## get number of exons per gene

exon_counts_only_in_setdiff <- gtf |> as.data.frame() |>
  filter(type == "exon") |>
  filter(gene_id %in% setdiff_genes) |>
  group_by(gene_id) |>
  summarise(exon_count = n()) |>
  as.data.frame()

all_exon_counts <- gtf |> as.data.frame() |>
  filter(type == "exon") |>
  filter(gene_id %in% rownames(gene_counts)) |>
  group_by(gene_id) |>
  summarise(exon_count = n()) |>
  as.data.frame()


isoform_counts_only_in_setdiff <- gtf |> as.data.frame() |>
  filter(type == "transcript") |>
  filter(gene_id %in% setdiff_genes) |>
  group_by(gene_id) |>
  summarise(isoform_count = n()) |>
  as.data.frame()

all_isoform_counts <- gtf |> as.data.frame() |>
  filter(type == "transcript") |>
  filter(gene_id %in% rownames(gene_counts)) |>
  group_by(gene_id) |>
  summarise(isoform_count = n()) |>
  as.data.frame()

median(exon_counts_only_in_setdiff$exon_count)
range(exon_counts_only_in_setdiff$exon_count)
median(isoform_counts_only_in_setdiff$isoform_count)
range(isoform_counts_only_in_setdiff$isoform_count)


median(all_exon_counts$exon_count)
range(all_exon_counts$exon_count)
median(all_isoform_counts$isoform_count)
range(all_isoform_counts$isoform_count)



##### 
fasta_path <- file.path(bambu_dir, "extended_annotations.fasta")
fasta <- Biostrings::readDNAStringSet(fasta_path)

# Look at transcript names
names(fasta) <- gsub("\\..*", "", names(fasta))
names(fasta) <- gsub(" .*", "", names(fasta))

###isoforms that were identified

isoforms_identified_fa <- fasta[which(names(fasta) %in% 
                                        (isoformFeatures$isoform_id |> unique() ))]

isoforms_identified_df <- data.frame(
  transcript_id = names(isoforms_identified_fa),
  length = width(isoforms_identified_fa)
)

median(isoforms_identified_df$length)


### isoforms not identified 

isoforms_not_identified <- all_isoform_counts <- gtf |> as.data.frame() |>
  filter(type == "transcript") 
isoforms_not_identified$transcript_id <- gsub("\\..*", "", isoforms_not_identified$transcript_id)

isoforms_not_identified <- isoforms_not_identified |> filter(gene_id %in% rownames(gene_counts)) |> select(transcript_id, gene_id) |> distinct() |>
  filter(!transcript_id %in% isoformFeatures$isoform_id)


isoforms_not_identified_fa <- fasta[which(names(fasta) %in% 
                                            isoforms_not_identified$transcript_id)]

isoforms_not_identified_df <- data.frame(
  transcript_id = names(isoforms_not_identified_fa),
  length = width(isoforms_not_identified_fa)
)


#### compare ranges #####
median(isoforms_not_identified_df$length)
range(isoforms_not_identified_df$length)

median(isoforms_identified_df$length)
range(isoforms_identified_df$length)
