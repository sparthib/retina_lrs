library(readr)
library(dplyr)
library(biomaRt)
library(edgeR)
library(ggplot2)

# get all the files called "counts_gene.txt" in all subdirs 
files <- list.files(path = "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/H9_EP1_bambu",
                    pattern = "counts_gene.txt",
                    recursive = TRUE, full.names = TRUE)

# get sample name from the path
get_sample_name <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  sample <- parts[length(parts) - 2]
  haplotype <- parts[length(parts) - 1]
  paste0(sample, "_", haplotype)
}

samples <- sapply(files, get_sample_name)
samples 
names(samples) <- NULL

# read all the counts files into a list of data frames
counts_list <- lapply(files, read_tsv, col_types = cols())
names(counts_list) <- samples

#change colnames to c("gene_id", "sample")
for(i in seq_along(counts_list)) {
  colnames(counts_list[[i]]) <- c("gene_id", names(counts_list)[i])
}

# check nrow for each data frame
nrows <- sapply(counts_list, nrow)

# all data frames should have the same gene_id order
all_same <- all(sapply(counts_list, function(df) all(df$gene_id == counts_list[[1]]$gene_id)))

#cbind all data frames by gene_id
if(all_same) {
  counts_matrix <- Reduce(function(x, y) merge(x, y, by = "gene_id"), counts_list)
} else {
  stop("Gene IDs do not match across all samples.")
}


# remove version number for gene_id column 
counts_matrix$gene_id <- gsub("\\..*$", "", counts_matrix$gene_id)


#remove zero value rows
counts_matrix <- counts_matrix[rowSums(counts_matrix[,-1]) > 0, ]

nrow(counts_matrix)

# only keep PTC genes 
mart <- useMart(biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl")

# Fetch gene symbols using getBM
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                     filters = "ensembl_gene_id",
                     values = counts_matrix$gene_id,
                     mart = mart)

#filter for protein coding genes
ptc_genes <- gene_info |> 
  filter(gene_biotype == "protein_coding") |> 
  pull(ensembl_gene_id)

counts_matrix_ptc <- counts_matrix |> 
  filter(gene_id %in% ptc_genes)

# nrow(counts_matrix_ptc)
# [1] 16966

#convert to cpm and do edgeR expression filtering
alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")


groups <- c("Stage_3", "Stage_3", 
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "Stage_2", "Stage_2",
            "Stage_3", "Stage_3",
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "FT", "FT",  
            "FT", "FT",
            "RGC", "RGC",
            "RGC", "RGC")

groups <- paste0(groups,"_", alleles)

colnames(counts_matrix_ptc)[-1] <- samples

min_counts <- 3 #since ASE is reduced to ~25% of all alignments
dge <- DGEList(counts = counts_matrix_ptc[,-1], samples = samples,
               group = groups, gene = counts_matrix_ptc$gene_id)
keep <- filterByExpr(dge, group = groups, min.count = min_counts)
dge <- dge[keep, ]
dge <- calcNormFactors(dge, method ="TMM")
# Convert to CPM (Counts Per Million)
filtered_cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
rownames(filtered_cpm_matrix ) <- dge$genes$genes
filtered_counts_matrix <- dge$counts 
rownames(filtered_counts_matrix) <- dge$genes$genes

dir.create("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/", recursive = TRUE, showWarnings = FALSE)

readr::write_tsv(as.data.frame(filtered_cpm_matrix) |> 
                         mutate(gene_id = rownames(filtered_cpm_matrix)),
                       file = "/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_cpm.tsv")
readr::write_tsv(as.data.frame(filtered_counts_matrix) |> 
                         mutate(gene_id = rownames(filtered_counts_matrix)),
                       file = "/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_counts.tsv")

