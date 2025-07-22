library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Rsubread)
library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "EP1-WT_ROs_D45", "EP1-BRN3B-RO", "EP1-WT_hRO_2") 
gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/gene_counts_RO_samples"

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

# counts_matrix |> colnames()
# "EP1.BRN3B.RO_h1"   "EP1.BRN3B.RO_h2"  
#"EP1.WT_hRO_2_h1"   [4] "EP1.WT_hRO_2_h2" 
#"EP1.WT_ROs_D45_h1" "EP1.WT_ROs_D45_h2"
# "H9.BRN3B_hRO_2_h1" "H9.BRN3B_hRO_2_h2" 
#"H9.BRN3B.RO_h1"    "H9.BRN3B.RO_h2"    
# "H9.CRX_hRO_2_h1"   "H9.CRX_hRO_2_h2"  
#  "H9.CRX_ROs_D45_h1" "H9.CRX_ROs_D45_h2"

groups <- c("H1_Stage3", "H2_Stage_3", 
            "H1_Stage2", "H2_Stage2",
            "H1_Stage1", "H2_Stage1",
            "H1_Stage2", "H2_Stage2",
            "H1_Stage3", "H2_Stage3",
            "H1_Stage2", "H2_Stage2",
            "H1_Stage1", "H2_Stage1")

alleles <- c("H1", "H2", "H1", "H2", 
            "H1", "H2", "H1", "H2",
            "H1", "H2", "H1", "H2",
            "H1", "H2")
stages <- c("Stage3", "Stage3", 
            "Stage2", "Stage2",
            "Stage1", "Stage1", 
            "Stage2", "Stage2",
            "Stage3", "Stage3",
            "Stage2", "Stage2",
            "Stage1", "Stage1")

design <- data.frame(
  sample = colnames(counts_matrix),
  group = groups,
  allele = factor(alleles, levels = c("H1", "H2")),
  stage = factor(stages, levels = c("Stage1", "Stage2", "Stage3"))
)

design <- model.matrix(~ allele*stage, data = design)

# > colnames(design)
# [1] "(Intercept)"          "alleleH2"             "stageStage2"         
# [4] "stageStage3"          "alleleH2:stageStage2" "alleleH2:stageStage3"

#coefficient 2: H2 vs. H1 at stage 1 
#coefficient 3: Stage 2 vs. Stage 1
#coefficient 4: Stage 3 vs. Stage 1
#coefficient 5: (H2 - H1 at Stage2) - (H2 - H1 at Stage1)
#coefficient 6: (H2 - H1 at Stage3) - (H2 - H1 at Stage2)


y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

qlf_H1_vs_H2_stage_1 <- glmQLFTest(fit, coef=2) # H2 vs H1 across all samples 
qlf_Stage_2_vs_Stage_1_H1 <- glmQLFTest(fit, coef=3) # Stage 2 vs Stage 1 across all samples
qlf_Stage_3_vs_Stage_1_H1 <- glmQLFTest(fit, coef=4) # Stage 3 vs Stage 1 across all samples
qlf_H2_vs_H1_stage2 <- glmQLFTest(fit, coef=5) # (H2 - H1 at Stage2) - (H2 - H1 at Stage1)
qlf_H2_vs_H1_stage3 <- glmQLFTest(fit, coef=6) # (H2 - H1 at Stage3) - (H2 - H1 at Stage2)
# Get top tags for each comparison
top_H1_vs_H2_stage_1 <- topTags(qlf_H1_vs_H2_stage_1, n = Inf)$table
top_Stage_2_vs_Stage_1_H1 <- topTags(qlf_Stage_2_vs_Stage_1_H1, n = Inf)$table
top_Stage_3_vs_Stage_1_H1 <- topTags(qlf_Stage_3_vs_Stage_1_H1, n = Inf)$table
top_H2_vs_H1_stage2 <- topTags(qlf_H2_vs_H1_stage2, n = Inf)$table
top_H2_vs_H1_stage3 <- topTags(qlf_H2_vs_H1_stage3, n = Inf)$table

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/ROs/"
dir.create(dge_output_dir, showWarnings = FALSE)


add_gene_names <- function(tt, comparison_name) {

  tt$gene_id <- gsub("\\..*", "", tt$genes)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"), 
                       filter="ensembl_gene_id", 
                       values=tt$gene_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype", "chromosome_name")
  
  tt <- merge(tt, annotLookup, by="gene_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  tt$neg_log10_FDR <- -log10(tt$FDR)
  tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
  tt <- tt |>
    arrange(desc(abs(logFC)), FDR)
  tt$label <- NA
  tt$label[1:20] <- tt$gene_name[1:20]
  file <- paste0(dge_output_dir, comparison_name, "_DGEs.tsv")
  write_tsv(tt, file)
}

# Save results to files
add_gene_names(top_H1_vs_H2_stage_1,  "H1_vs_H2_stage_1")
add_gene_names(top_Stage_2_vs_Stage_1_H1, "Stage_2_vs_Stage_1_H1")
add_gene_names(top_Stage_3_vs_Stage_1_H1, "Stage_3_vs_Stage_1_H1")
add_gene_names(top_H2_vs_H1_stage2, "H2_vs_H1_stage2")
add_gene_names(top_H2_vs_H1_stage3,  "H2_vs_H1_stage3")

###### manual #########
y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             group = groups,
             genes = rownames(counts_matrix))

y <- normLibSizes(y)


colnames(design) <- gsub("group", "", colnames(design))

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

# Fit the model and create contrasts
fit <- glmQLFit(y, design, robust=TRUE)

                       #within stage, between allele
contr <- makeContrasts(H1_Stage1_vs_H2_Stage1 = H2_Stage1 - H1_Stage1, 
                      H1_Stage2_vs_H2_Stage2 = H2_Stage2 - H1_Stage2,
                      H1_Stage3_vs_H2_Stage3 = H2_Stage3 - H1_Stage3,
                      #between stage, within allele
                      H1_Stage1_vs_H1_Stage2 = H1_Stage2 - H1_Stage1,
                      H1_Stage2_vs_H1_Stage3 = H1_Stage3 - H1_Stage2,
                      H1_Stage1_vs_H1_Stage3 = H1_Stage3 - H1_Stage1,
                      H2_Stage1_vs_H2_Stage2 = H2_Stage2 - H2_Stage1,
                      H2_Stage2_vs_H2_Stage3 = H2_Stage3 - H2_Stage2,
                      H2_Stage1_vs_H2_Stage3 = H2_Stage3 - H2_Stage1,
                       levels=design)



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
  file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  write_tsv(tt, file)
}


for (i in seq_len(ncol(contr))) {
  read_file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  tt <- read_tsv(read_file)
  pdf(file = paste0(dge_output_dir, colnames(contr)[i], "_volcano_plot.pdf"),
      width = 8, height = 6)
  
  print(
    ggplot(tt, aes(x = logFC, y = neg_log10_FDR)) +
      geom_point(aes(color = significant), size = 1.5) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      scale_color_manual(values = c("gray", "red")) +
      labs(title = colnames(contr)[i],
           x = "log2 Fold Change",
           y = "-log10 FDR") +
      theme_minimal()
  )
  
  dev.off()
}

ora_plot <- function(genelist, ont, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = names(genelist),
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = ont,
                  pAdjustMethod = "fdr",
                  minGSSize     = 100,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01,
                  readable      = TRUE) 
  ego <- enrichplot::pairwise_termsim(ego)
  ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  if(nrow(as.data.frame(ego2)) != 0){
    write_tsv(as.data.frame(ego2), file.path(output_plot_dir,
                                             paste0("simplified_ORA_ASE_DGE_genes_",analysis_type, ont, ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("simplified_ORA_all_ASE_DGE_genes_",analysis_type, ont, ".pdf")))
    print(dotplot(ego2, showCategory = 15))
    dev.off()
    
  }
}

go_plot_dir <- file.path(dge_output_dir, "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


for (i in seq_len(ncol(contr))) {
  read_file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  tt <- read_tsv(read_file)
  tt <- tt |> filter(significant == TRUE)
  
  # GO plot 
  names <- tt |> pull(gene_id)
  values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
  names(values) <- names
  head(values)
  
  ora_plot(genelist = values, 
           ont = "BP", 
           output_plot_dir = go_plot_dir,
           analysis_type = colnames(contr)[i])
}
