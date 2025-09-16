library(biomaRt)
library(readr)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(edgeR)

ora_plot <- function(genelist, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = names,
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  readable      = TRUE) 
  if(nrow(as.data.frame(ego)) != 0){
    ego <- enrichplot::pairwise_termsim(ego)
    # ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    write_tsv(as.data.frame(ego), file.path(output_plot_dir,
                                            paste0("simplified_ORA_ASE_DGE_genes_",analysis_type, "BP", ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("simplified_ORA_all_ASE_DGE_genes_",analysis_type, "BP", ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  } else {
    message("No significant GO terms found for ", analysis_type, " in ", "BP")}
}

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/H9_vs_EP1/filter_counts"
dir.create(dge_output_dir, recursive = TRUE, showWarnings = FALSE)

go_plot_dir <- file.path(dge_output_dir, "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2",
             "EP1-WT_ROs_D45", "EP1-BRN3B-RO", "EP1-WT_hRO_2") 

gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/H9_EP1_gene_counts_all_samples"

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

# only keep PTC genes in the OG counts matrix
nrow(counts_matrix)

# remove rows that are not in the OG counts matrix
non_ase_counts_matrix <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype"
non_ase_counts_matrix <- readRDS(file.path(non_ase_counts_matrix, "filtered_gene_counts.RDS"))
nrow(non_ase_counts_matrix)

counts_matrix <- counts_matrix[rownames(counts_matrix) %in% rownames(non_ase_counts_matrix), ]
nrow(counts_matrix)
# [1] 15799


culture_type <- c("EP1", "EP1", "EP1", "EP1", "EP1", "EP1",
                  "H9", "H9", "H9", "H9", "H9", "H9",
                  "H9", "H9", "H9", "H9", "H9", "H9",
                  "H9", "H9", "H9", "H9")

alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")

groups <- c("Stage3", "Stage3", 
            "Stage2", "Stage2",
            "Stage1", "Stage1",
            "Stage2", "Stage2",
            "Stage3", "Stage3",
            "Stage2", "Stage2",
            "Stage1", "Stage1",
            "FT", "FT",  
            "FT", "FT",
            "RGC", "RGC",
            "RGC", "RGC")

column_data <- data.frame(
  row.names = colnames(counts_matrix),
  culture_type = culture_type,
  allele = alleles,
  group = groups
)

design <- data.frame(
  sample = colnames(counts_matrix),
  group = factor(groups, levels = c("Stage1", "Stage2", "Stage3", "FT", "RGC")),
  allele = factor(alleles, levels = c("H1", "H2")),
  culture_type  = factor(culture_type , levels = c("H9", "EP1"))
)

design <- model.matrix(~ 0 + group*allele + allele*culture_type + group*allele,
                       data = design)
# > colnames(design)
# [1] "groupStage1"              "groupStage2"             
# [3] "groupStage3"              "groupFT"                 
# [5] "groupRGC"                 "alleleH2"                
# [7] "culture_typeEP1"          "groupStage2:alleleH2"    
# [9] "groupStage3:alleleH2"     "groupFT:alleleH2"        
# [11] "groupRGC:alleleH2"        "alleleH2:culture_typeEP1"

### H1 vs. H2 is coef 6
### H9 vs EP1 is coef 7

y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             group = groups,
             allele = alleles,
             genes = rownames(counts_matrix))

# mincount = 3 given the seq depth is about 30% of the non-ASE specific samples. 
keep <- filterByExpr(y, min.count = 3)  # or other threshold
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)

# > nrow(y)
# [1] 12503


y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)

qlf_allele <- glmQLFTest(fit, coef=6)
qlf_cell_line <- glmQLFTest(fit, coef = 7)

is.de_allele <- decideTests(qlf_allele, p.value=0.05)
is.de_cell_line <- decideTests(qlf_cell_line, p.value=0.05)

# summary(is.de_allele)
# summary(is.de_cell_line)
summary(is.de_allele)
# alleleH2
# Down        202
# NotSig    12169
# Up          132
summary(is.de_cell_line)
# culture_typeEP1
# Down              2900
# NotSig            7417
# Up                2186


# Extract the results
tt_allele <- topTags(qlf_allele,n = Inf)
nrow(tt_allele) 
head(tt_allele)
colnames(tt_allele$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")


tt_allele$condition_1 <- "H1"
tt_allele$condition_2 <- "H2"

tt_cell_line <- topTags(qlf_cell_line, n = Inf)
nrow(tt_cell_line)
colnames(tt_cell_line$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")
tt_cell_line$condition_1 <- "H9"
tt_cell_line$condition_2 <- "EP1"


merge_table <- function(tt, annotLookup) {
  # Merge the annotation with the results
  merged_results <- merge(tt$table, annotLookup,
                          by = "gene_id", all.x = TRUE)
  
  # Filter significant results
  merged_results_sig <- merged_results |> 
    dplyr::filter(FDR < 0.05 & abs(logFC) >= 1)
  
  # Add -log10 FDR for plotting
  merged_results <- merged_results |> 
    mutate(neg_log10_FDR = -log10(FDR))
  
  # Separate out the first 20 genes for labeling
  merged_results$label <- NA
  merged_results <- merged_results |>
    arrange(desc(abs(logFC)), FDR)
  merged_results$label[1:20] <- merged_results$gene_name[1:20]
  
  merged_results$significant <- merged_results$FDR < 0.05 & 
    abs(merged_results$logFC) > 1
  
  return(merged_results)
}
  
merged_results_allele <- merge_table(tt_allele, annotLookup)
merged_results_cell_line <- merge_table(tt_cell_line, annotLookup)

write_tsv(merged_results_allele,
          file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_DE_results.tsv"))
write_tsv(merged_results_cell_line,
          file = file.path(dge_output_dir, "ASE_H9_vs_EP1_PTC_DE_results.tsv"))

write_tsv(counts_matrix,
          file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_gene_counts.tsv"))

# Volcano plot
pdf(file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_DE_volcano_plot.pdf"),
    width = 8, height = 6)

ggplot(merged_results_allele, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "H1 vs H2, all samples",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()
dev.off()

pdf(file = file.path(dge_output_dir, "ASE_H9_vs_EP1_PTC_DE_volcano_plot.pdf"),
    width = 8, height = 6)
ggplot(merged_results_cell_line, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "H9 vs EP1, all samples",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()
dev.off()

# GO plot 
tt <- merged_results_allele
tt <- tt |> filter(significant == TRUE)
names <- tt |> pull(gene_id)
values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
names(values) <- names
head(values)

ora_plot(genelist = names, 
         output_plot_dir = go_plot_dir,
         analysis_type = "H1_vs_H2")

tt <- merged_results_cell_line
tt <- tt |> filter(significant == TRUE)
names <- tt |> pull(gene_id)
values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
names(values) <- names
head(values)

ora_plot(genelist = names, 
         output_plot_dir = go_plot_dir,
         analysis_type = "H9_vs_EP1")
