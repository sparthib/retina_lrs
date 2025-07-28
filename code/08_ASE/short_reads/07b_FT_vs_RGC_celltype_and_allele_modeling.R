library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

samples <- c("H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2") 
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

# only keep PTC genes

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( 
    "external_gene_name",
    "gene_biotype", "ensembl_gene_id","chromosome_name"),
  filter="ensembl_gene_id",
  values=rownames(counts_matrix),
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id", "chromosome_name")
annotLookup <- annotLookup |> dplyr::distinct()
annotLookup |> nrow()
annotLookup <- annotLookup |> 
  dplyr::filter(gene_biotype == "protein_coding") 
annotLookup |> nrow()

# filter for protein coding genes
counts_matrix <- counts_matrix[rownames(counts_matrix) 
                               %in% annotLookup$gene_id, ]

#FT vs RGC samples 
counts_matrix <- counts_matrix[, 15:22]

sample_names <- colnames(counts_matrix)

groups <- c("H1_FT", "H2_FT", "H1_FT", "H2_FT",
            "H1_RGC", "H2_RGC", "H1_RGC", "H2_RGC")
cell_type <- c("FT", "FT", "FT", "FT",
               "RGC", "RGC", "RGC", "RGC")
allele_type <- c("H1", "H2", "H1", "H2",
                 "H1", "H2", "H1", "H2")

library(edgeR)

y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             genes = rownames(counts_matrix),
             group = groups)

y <- normLibSizes(y)

design <- model.matrix(~ 0 + group, 
                       data = y$samples)

colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
# 0.02222802

contr <- makeContrasts(H1_FT_vs_H2_FT = H2_FT - H1_FT, 
                       H1_RGC_vs_H2_RGC = H2_RGC - H1_RGC, 
                       H1_FT_vs_H1_RGC = H1_RGC - H1_FT,
                       H2_FT_vs_H2_RGC = H2_RGC - H2_FT,
                       levels=design)
fit <- glmQLFit(y, design, robust=TRUE)

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/FT_vs_RGC"
if (!dir.exists(dge_output_dir)) {
  dir.create(dge_output_dir, recursive = TRUE)
}


for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$gene_id <- gsub("\\..*", "", tt$genes)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"), 
                       filter="ensembl_gene_id", 
                       values=tt$gene_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype", "chromosome_name")
  
  tt <- merge(tt, annotLookup, by="gene_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == -1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == 1]
  tt$neg_log10_FDR <- -log10(tt$FDR)
  tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
  tt <- tt |>
    arrange(desc(abs(logFC)), FDR)
  tt$label <- NA
  tt$label[1:20] <- tt$gene_name[1:20]
  file <- file.path(dge_output_dir, paste0( colnames(contr)[i], "_DGEs.tsv"))
  write_tsv(tt, file)
}


for (i in seq_len(ncol(contr))) {
  read_file <- file.path(dge_output_dir, paste0( colnames(contr)[i], "_DGEs.tsv"))
  tt <- read_tsv(read_file)
  pdf(file = file.path(dge_output_dir, paste0(colnames(contr)[i], "_volcano_plot.pdf")),
      width = 8, height = 6)
  ggplot(tt, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = colnames(contr)[i],
         x = "log2 Fold Change",
         y = "-log10 FDR") +
    theme_minimal()
  dev.off()
}

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


go_plot_dir <- file.path(dge_output_dir, "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


for (i in seq_len(ncol(contr))) {
  read_file <- file.path(dge_output_dir, paste0(colnames(contr)[i], "_DGEs.tsv"))
  tt <- read_tsv(read_file)
  tt <- tt |> filter(significant == TRUE)
  
  # GO plot 
  names <- tt |> pull(gene_id)
  values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
  names(values) <- names

  ora_plot(genelist = names, 
           output_plot_dir = go_plot_dir,
           analysis_type = colnames(contr)[i])
}

