library(here)
library(VennDiagram)


#significance means FDR < 0.05 and abs(log2FC) > 1
non_ASE_RO_dir <-file.path(code_dir, "processed_data", "dtu", "bambu", "ROs",
                       "protein_coding", "DGE")
non_ASE_FT_RGC_dir <-file.path(code_dir, "processed_data", "dtu", "bambu", "FT_vs_RGC",
                           "protein_coding") 


## NON-ASE

# Stage 1 vs. Stage 2 
non_ase_stage_1_vs_stage_2 <- read.table(file.path(non_ASE_RO_dir,
                                                  "D45_vs_D100_DGEs.tsv"),
                                            header = TRUE,
                                            sep = "\t",
                                            stringsAsFactors = FALSE)

# Stage 2 vs. Stage 3
non_ase_stage_2_vs_stage_3 <- read.table(file.path(non_ASE_RO_dir,
                                                  "D100_vs_D200_DGEs.tsv"),
                                            header = TRUE,
                                            sep = "\t",
                                            stringsAsFactors = FALSE)

# Stage 1 vs. Stage 3 
non_ase_stage_1_vs_stage_3 <- read.table(file.path(non_ASE_RO_dir,
                                                  "D45_vs_D200_DGEs.tsv"),
                                            header = TRUE,
                                            sep = "\t",
                                            stringsAsFactors = FALSE)


# FT vs. RGC
non_ase_FT_vs_RGC <- read.table(file.path(non_ASE_FT_RGC_dir,
                                           "DGE_table.tsv"),
                                     header = TRUE,
                                     sep = "\t",
                                     stringsAsFactors = FALSE)



H1_dir <-file.path(code_dir, "processed_data", "ASE", 
               "bambu_counts_matrices", "DGE", "between_stage_in_H1_allele")
H2_dir <-file.path(code_dir, "processed_data", "ASE", 
               "bambu_counts_matrices", "DGE", "between_stage_in_H2_allele")


H1_stage_1_vs_stage_2 <- read.table(file.path(H1_dir,
                                              "H1_Stage1_vs_H1_Stage2_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H1_stage_2_vs_stage_3 <- read.table(file.path(H1_dir,
                                              "H1_Stage2_vs_H1_Stage3_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H1_stage_1_vs_stage_3 <- read.table(file.path(H1_dir,
                                              "H1_Stage1_vs_H1_Stage3_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H1_FT_vs_RGC <- read.table(file.path(H1_dir,
                                       "H1_FT_vs_H1_RGC_DGEs.tsv"),
                                 header = TRUE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)




H2_stage_1_vs_stage_2 <- read.table(file.path(H2_dir,
                                              "H2_Stage1_vs_H2_Stage2_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H2_stage_2_vs_stage_3 <- read.table(file.path(H2_dir,
                                              "H2_Stage2_vs_H2_Stage3_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H2_stage_1_vs_stage_3 <- read.table(file.path(H2_dir,
                                              "H2_Stage1_vs_H2_Stage3_DGEs.tsv"),
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)
H2_FT_vs_RGC <- read.table(file.path(H2_dir,
                                       "H2_FT_vs_H2_RGC_DGEs.tsv"),
                                 header = TRUE,
                                 sep = "\t",
                                 stringsAsFactors = FALSE)

get_list_of_sig_genes <- function(df){
  sig_genes <- df[df$FDR < 0.05 & abs(df$logFC) > 1, "gene_id"]
  return(sig_genes)
}


# Get lists of significant genes for each comparison
non_ase_sig_stage_1_vs_stage_2_genes <- get_list_of_sig_genes(non_ase_stage_1_vs_stage_2)
non_ase_sig_stage_2_vs_stage_3_genes <- get_list_of_sig_genes(non_ase_stage_2_vs_stage_3)
non_ase_sig_stage_1_vs_stage_3_genes <- get_list_of_sig_genes(non_ase_stage_1_vs_stage_3)
non_ase_sig_FT_vs_RGC_genes <- get_list_of_sig_genes(non_ase_FT_vs_RGC)
H1_sig_stage_1_vs_stage_2_genes <- get_list_of_sig_genes(H1_stage_1_vs_stage_2)
H1_sig_stage_2_vs_stage_3_genes <- get_list_of_sig_genes(H1_stage_2_vs_stage_3)
H1_sig_stage_1_vs_stage_3_genes <- get_list_of_sig_genes(H1_stage_1_vs_stage_3)
H1_sig_FT_vs_RGC_genes <- get_list_of_sig_genes(H1_FT_vs_RGC)


H2_sig_stage_1_vs_stage_2_genes <- get_list_of_sig_genes(H2_stage_1_vs_stage_2)
H2_sig_stage_2_vs_stage_3_genes <- get_list_of_sig_genes(H2_stage_2_vs_stage_3)
H2_sig_stage_1_vs_stage_3_genes <- get_list_of_sig_genes(H2_stage_1_vs_stage_3)
H2_sig_FT_vs_RGC_genes <- get_list_of_sig_genes(H2_FT_vs_RGC)


# Create the Venn diagram
venn.plot.dir <-file.path(code_dir, "processed_data", "ASE", "bambu_counts_matrices","venn_plots")
dir.create(venn.plot.dir, recursive = TRUE, showWarnings = FALSE)

venn.diagram(
  x = list(
    non_ase = non_ase_sig_stage_1_vs_stage_2_genes,
    H1 = H1_sig_stage_1_vs_stage_2_genes,
    H2 = H2_sig_stage_1_vs_stage_2_genes
  ),
  filename = file.path(venn.plot.dir,
                               "sig_genes_overlap_stage_1_vs_2.png"),
  category.names = c("non_ase", "H1", "H2"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("red", "purple", "#fde725ff"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(
    non_ase = non_ase_sig_stage_2_vs_stage_3_genes,
    H1 = H1_sig_stage_2_vs_stage_3_genes,
    H2 = H2_sig_stage_2_vs_stage_3_genes
  ),
  filename = file.path(venn.plot.dir,
                               "sig_genes_overlap_stage_2_vs_3.png"),
  category.names = c("non_ase", "H1", "H2"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("red", "purple", "#fde725ff"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(
    non_ase = non_ase_sig_stage_1_vs_stage_3_genes,
    H1 = H1_sig_stage_1_vs_stage_3_genes,
    H2 = H2_sig_stage_1_vs_stage_3_genes
  ),
  filename = file.path(venn.plot.dir,
                               "sig_genes_overlap_stage_1_vs_3.png"),
  category.names = c("non_ase", "H1", "H2"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("red", "purple", "#fde725ff"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(
    non_ase = non_ase_sig_FT_vs_RGC_genes,
    H1 = H1_sig_FT_vs_RGC_genes,
    H2 = H2_sig_FT_vs_RGC_genes
  ),
  filename = file.path(venn.plot.dir,
                               "sig_genes_overlap_FT_vs_RGC.png"),
  category.names = c("non_ase", "H1", "H2"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("red", "purple", "#fde725ff"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)



# Stage 1 vs Stage 2
stage_1_vs_2_genes <- unique(c(
  non_ase_sig_stage_1_vs_stage_2_genes,
  H1_sig_stage_1_vs_stage_2_genes,
  H2_sig_stage_1_vs_stage_2_genes
))

stage_1_vs_2_df <- data.frame(gene = stage_1_vs_2_genes)
stage_1_vs_2_df$non_ase <- ifelse(stage_1_vs_2_df$gene %in% non_ase_sig_stage_1_vs_stage_2_genes, "yes", "no")
stage_1_vs_2_df$H1 <- ifelse(stage_1_vs_2_df$gene %in% H1_sig_stage_1_vs_stage_2_genes, "yes", "no")
stage_1_vs_2_df$H2 <- ifelse(stage_1_vs_2_df$gene %in% H2_sig_stage_1_vs_stage_2_genes, "yes", "no")

# Stage 2 vs Stage 3
stage_2_vs_3_genes <- unique(c(
  non_ase_sig_stage_2_vs_stage_3_genes,
  H1_sig_stage_2_vs_stage_3_genes,
  H2_sig_stage_2_vs_stage_3_genes
))

stage_2_vs_3_df <- data.frame(gene = stage_2_vs_3_genes)
stage_2_vs_3_df$non_ase <- ifelse(stage_2_vs_3_df$gene %in% non_ase_sig_stage_2_vs_stage_3_genes, "yes", "no")
stage_2_vs_3_df$H1 <- ifelse(stage_2_vs_3_df$gene %in% H1_sig_stage_2_vs_stage_3_genes, "yes", "no")
stage_2_vs_3_df$H2 <- ifelse(stage_2_vs_3_df$gene %in% H2_sig_stage_2_vs_stage_3_genes, "yes", "no")

# Stage 1 vs Stage 3
stage_1_vs_3_genes <- unique(c(
  non_ase_sig_stage_1_vs_stage_3_genes,
  H1_sig_stage_1_vs_stage_3_genes,
  H2_sig_stage_1_vs_stage_3_genes
))

stage_1_vs_3_df <- data.frame(gene = stage_1_vs_3_genes)
stage_1_vs_3_df$non_ase <- ifelse(stage_1_vs_3_df$gene %in% non_ase_sig_stage_1_vs_stage_3_genes, "yes", "no")
stage_1_vs_3_df$H1 <- ifelse(stage_1_vs_3_df$gene %in% H1_sig_stage_1_vs_stage_3_genes, "yes", "no")
stage_1_vs_3_df$H2 <- ifelse(stage_1_vs_3_df$gene %in% H2_sig_stage_1_vs_stage_3_genes, "yes", "no")

# FT vs RGC
FT_vs_RGC_genes <- unique(c(
  non_ase_sig_FT_vs_RGC_genes,
  H1_sig_FT_vs_RGC_genes,
  H2_sig_FT_vs_RGC_genes
))

FT_vs_RGC_df <- data.frame(gene = FT_vs_RGC_genes)
FT_vs_RGC_df$non_ase <- ifelse(FT_vs_RGC_df$gene %in% non_ase_sig_FT_vs_RGC_genes, "yes", "no")
FT_vs_RGC_df$H1 <- ifelse(FT_vs_RGC_df$gene %in% H1_sig_FT_vs_RGC_genes, "yes", "no")
FT_vs_RGC_df$H2 <- ifelse(FT_vs_RGC_df$gene %in% H2_sig_FT_vs_RGC_genes, "yes", "no")

readr::write_tsv(stage_1_vs_2_df,
                   file.path(venn.plot.dir,
                             "sig_genes_overlap_stage_1_vs_2.tsv"))
readr::write_tsv(stage_2_vs_3_df,
                   file.path(venn.plot.dir,
                             "sig_genes_overlap_stage_2_vs_3.tsv"))
readr::write_tsv(stage_1_vs_3_df,
                   file.path(venn.plot.dir,
                             "sig_genes_overlap_stage_1_vs_3.tsv"))
readr::write_tsv(FT_vs_RGC_df,
                   file.path(venn.plot.dir,
                             "sig_genes_overlap_FT_vs_RGC.tsv"))






