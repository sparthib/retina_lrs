library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)

samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
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
# nrow(counts_matrix)
# 15678

colnames(counts_matrix) <- gsub(".bam", "", colnames(counts_matrix))

#remove version number from row names
rownames(counts_matrix) <- gsub("\\..*", "", rownames(counts_matrix))

# only keep PTC genes in non ASE matrix

non_ase_matrix <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/bambu/ROs/filtered_by_counts_and_biotype"
non_ase_matrix <- readRDS(file.path(non_ase_matrix, "filtered_gene_counts.RDS"))
rownames(non_ase_matrix) <- gsub("\\..*", "", rownames(non_ase_matrix))

counts_matrix <- counts_matrix[rownames(counts_matrix) %in% rownames(non_ase_matrix), ]
# nrow(counts_matrix)
# 15678
counts_matrix <- counts_matrix[, c(1:14)]


cell_lines <- c("EP1", "EP1", "EP1",
                "EP1", "EP1", "EP1", 
               "H9", "H9", "H9", 
               "H9", "H9", "H9",
               "H9", "H9")

groups <- c("H1_Stage3", "H2_Stage3", 
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

# Create a proper samples data frame
samples_df <- data.frame(
  sample = colnames(counts_matrix),
  cell_line = factor(cell_lines, levels = c("H9", "EP1")),
  group = factor(groups),
  allele = factor(alleles, levels = c("H1", "H2")),
  stage = factor(stages, levels = c("Stage1", "Stage2", "Stage3"))
)

# Create DGEList with proper structure
y <- DGEList(
  counts = counts_matrix,
  samples = samples_df
)

# Now create the design matrix
design <- model.matrix(~ allele*stage + cell_line*stage + allele*cell_line, 
                                      data = y$samples)
keep <- filterByExpr(y, min.count = 3)
y <- y[keep,, keep.lib.sizes = FALSE]
# nrow(y)
# 13100
y <- normLibSizes(y)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)


dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/ROs/filtered_counts"
dir.create(dge_output_dir, showWarnings = FALSE)


# Fix the column names to make them syntactically valid
colnames(design) <- make.names(colnames(design))

# Check the new names
colnames(design)

# Now create the contrasts using the new names
contrasts <- makeContrasts(
  # Between allele, within stage (H2 - H1)
  H1_Stage1_vs_H2_Stage1 = alleleH2,
  H1_Stage2_vs_H2_Stage2 = alleleH2 + alleleH2.stageStage2,
  H1_Stage3_vs_H2_Stage3 = alleleH2 + alleleH2.stageStage3,
  
  # Between stage, within allele H1 (later stage - earlier stage)
  H1_Stage1_vs_H1_Stage2 = stageStage2,
  H1_Stage2_vs_H1_Stage3 = stageStage3 - stageStage2,
  H1_Stage1_vs_H1_Stage3 = stageStage3,
  
  # Between stage, within allele H2 (later stage - earlier stage)
  H2_Stage1_vs_H2_Stage2 = stageStage2 + alleleH2.stageStage2,
  H2_Stage2_vs_H2_Stage3 = stageStage3 + alleleH2.stageStage3 - stageStage2 - alleleH2.stageStage2,
  H2_Stage1_vs_H2_Stage3 = stageStage3 + alleleH2.stageStage3,
  
  levels = design
)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

for (i in seq_len(ncol(contrasts))) {
  qlf <- glmQLFTest(fit, contrast = contrasts[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  colnames(tt) <- c("logFC", "logCPM", "F" , "PValue", "FDR")
  tt$gene_id <- rownames(tt)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"), 
                       filter="ensembl_gene_id", 
                       values=tt$gene_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype", "chromosome_name")
  
  tt <- merge(tt, annotLookup, by="gene_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  
  tt$condition_1 <- gsub("_vs_.*", "", colnames(contrasts)[i])
  tt$condition_2 <- gsub(".*_vs_", "", colnames(contrasts)[i])
  tt$neg_log10_FDR <- -log10(tt$FDR)
  tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
  tt <- tt |>
    arrange(desc(abs(logFC)), FDR)
  tt$label <- NA
  tt$label[1:20] <- tt$gene_name[1:20]
  file <- file.path(dge_output_dir, paste0(colnames(contrasts)[i], "_DGEs.tsv"))
  write_tsv(tt, file)
}


source("/users/sparthib/retina_lrs/code/08_ASE/short_reads/02_feature_counts_analysis/volcano_helper.R")

contrast_names <- c(
  "H1_Stage1_vs_H2_Stage1",
  "H1_Stage2_vs_H2_Stage2", 
  "H1_Stage3_vs_H2_Stage3",
  "H1_Stage1_vs_H1_Stage2",
  "H1_Stage2_vs_H1_Stage3",
  "H1_Stage1_vs_H1_Stage3",
  "H2_Stage1_vs_H2_Stage2",
  "H2_Stage2_vs_H2_Stage3", 
  "H2_Stage1_vs_H2_Stage3"
)

# Generate volcano plots
generate_allele_volcano_plots(
  input_dir = dge_output_dir,
  contrast_names = contrast_names,
  table_type = "DGE"
)





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


go_plot_dir <- file.path(dge_output_dir, "plots", "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


for (i in seq_len(ncol(contrasts))) {
  read_file <- file.path(dge_output_dir, paste0(colnames(contrasts)[i], "_DGEs.tsv"))
  tt <- read_tsv(read_file)
  tt <- tt |> filter(significant == TRUE)
  
  # GO plot 
  names <- tt |> pull(gene_id)
  values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
  names(values) <- names
  
  ora_plot(genelist = names, 
           output_plot_dir = go_plot_dir,
           analysis_type = colnames(contrasts)[i])
}
