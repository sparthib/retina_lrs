library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(biomaRt)
library(readr)
library(dplyr)

samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2") 

gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gene_counts"

# load all counts matrix in the directory
files <- list.files(gene_counts_dir, pattern = "_counts.txt", full.names = TRUE)

# create a list of data frames

# Load all the count data into a list of data frames
counts_list <- lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1)  # Adjust if needed (e.g., use sep = "," for CSV)
})

# Combine all data frames into a single matrix by column binding them
counts_matrix <- do.call(cbind, counts_list)

## remove rows with 0 values

counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]

colnames(counts_matrix) <- gsub(".bam", "", colnames(counts_matrix))


groups <- c("H1", "H2", "H1", "H2",
            "H1", "H2", "H1", "H2",
            "H1", "H2", "H1", "H2",
            "H1", "H2", "H1", "H2")


library(edgeR)

y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             group = groups,
             genes = rownames(counts_matrix))

y <- normLibSizes(y)

design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(H2 - H1, levels=design)
qlf <- glmQLFTest(fit, contrast=contr)

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)

# > summary(is.de)
# -1*H1 1*H2
# Down          614
# NotSig      41063
# Up            291

# Extract the results
tt <- topTags(qlf,n = Inf)
nrow(tt) 
head(tt)
colnames(tt$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")


tt$table$gene_id <- gsub("\\..*", "", tt$table$gene_id)
tt$table <- tt$table[order(tt$table$FDR),]
tt$table$condition_1 <- names(contr[,1])[contr[,1] == -1]
tt$table$condition_2 <- names(contr[,1])[contr[,1] == 1]


mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( 
    "external_gene_name",
    "gene_biotype", "ensembl_gene_id"),
  filter="ensembl_gene_id",
  values=tt$table$gene_id,
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id")
annotLookup <- annotLookup |> dplyr::distinct()

# Merge the annotation with the results
merged_results <- merge(tt$table, annotLookup, by = "gene_id", all.x = TRUE)

merged_results_sig <- merged_results |> 
  dplyr::filter(FDR < 0.05 & abs(logFC) >= 1)

merged_results_sig$gene_biotype |> table()



# lncRNA                              miRNA
# 186                                 10
# misc_RNA               processed_pseudogene
# 14                                 57
# protein_coding                     rRNA
# 497                                  2
# rRNA_pseudogene                  snoRNA
# 3                                 35
# snRNA                                TEC
# 6                                 10
# transcribed_processed_pseudogene   transcribed_unitary_pseudogene
# 16                                  5
# transcribed_unprocessed_pseudogene  unprocessed_pseudogene
# 13                                  7



write_tsv(counts_matrix,
          file = file.path(gene_counts_dir, "H9_DNA_Seq_data_gene_counts.tsv"))
write_tsv(merged_results,
          file = file.path(gene_counts_dir, "H9_DNA_Seq_data_DE_results.tsv"))


#### visualization of ASE results (volcano plot)
library(ggplot2)
library(stringr)
library(ggrepel)
results <- merged_results
results$significant <- with(results, FDR < 0.05 & abs(logFC) > 1)

# remove pseudogenes and sort by FDR
results <- results %>% 
  filter(!str_detect(gene_biotype, "^_pseudogene")) %>%
  arrange(FDR) 
results$gene_biotype |> table()

# volcano plot
ggplot(results, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = results[results$significant, ], 
                  aes(label = gene_name), 
                  size = 3) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()

