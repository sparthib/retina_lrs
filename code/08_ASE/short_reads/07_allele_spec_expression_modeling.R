library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2") 

gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gene_counts"

# load all counts matrix in the directory
files <- list.files(gene_counts_dir,
                    pattern = "_counts.txt", 
                    full.names = TRUE)

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


#remove version number from row names
rownames(counts_matrix) <- gsub("\\..*", "", rownames(counts_matrix))

# only keep PTC genes

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( 
    "external_gene_name",
    "gene_biotype", "ensembl_gene_id"),
  filter="ensembl_gene_id",
  values=rownames(counts_matrix),
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id")
annotLookup <- annotLookup |> dplyr::distinct()
annotLookup |> nrow()
annotLookup <- annotLookup |> 
  dplyr::filter(gene_biotype == "protein_coding") 
annotLookup |> nrow()

# filter for protein coding genes
counts_matrix <- counts_matrix[rownames(counts_matrix) 
                        %in% annotLookup$gene_id, ]


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

# > summary(is.de) after only keeping PTC genes
# -1*H1 1*H2
# Down          371
# NotSig      16122
# Up            166

# Extract the results
tt <- topTags(qlf,n = Inf)
nrow(tt) 
head(tt)
colnames(tt$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")


tt$table <- tt$table[order(tt$table$FDR),]
tt$table$condition_1 <- names(contr[,1])[contr[,1] == -1]
tt$table$condition_2 <- names(contr[,1])[contr[,1] == 1]


# Merge the annotation with the results
merged_results <- merge(tt$table, annotLookup,
                        by = "gene_id", all.x = TRUE)

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

merged_results <- merged_results |>
  arrange(desc(abs(logFC)), FDR)

write_tsv(counts_matrix,
          file = file.path(gene_counts_dir, "H9_DNA_Seq_data_PTC_gene_counts.tsv"))

write_tsv(merged_results,
          file = file.path(gene_counts_dir, "H9_DNA_Seq_data_PTC_DE_results.tsv"))


head(merged_results, n = 10)



# Add -log10 FDR for plotting
merged_results <- merged_results |> 
  mutate(neg_log10_FDR = -log10(FDR))

# Separate out the first 20 genes for labeling
merged_results$label <- NA
merged_results$label[1:20] <- merged_results$gene_name[1:20]

merged_results$significant <- merged_results$FDR < 0.05 & abs(merged_results$logFC) > 1

# Volcano plot
pdf(file = file.path(gene_counts_dir, "H9_DNA_Seq_data_PTC_DE_volcano_plot.pdf"),
    width = 8, height = 6)

ggplot(merged_results, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()

dev.off()




