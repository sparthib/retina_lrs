library(here)
library(VennDiagram)

# > here::here()
# [1] "/users/sparthib/retina_lrs"

#significance means FDR < 0.05 and abs(log2FC) > 1
EP1_allele_DGE <- here("processed_data/ASE/bambu_counts_matrices/DGE/cell_line_vs_allele/H1_EP1_vs_H2_EP1_DGEs.tsv")
H9_allele_DGE <- here("processed_data/ASE/bambu_counts_matrices/DGE/cell_line_vs_allele/H1_H9_vs_H2_H9_DGEs.tsv")

EP1_allele_DGE <- read.table(EP1_allele_DGE,
                                         header = TRUE,
                                         sep = "\t",
                                         stringsAsFactors = FALSE)

H9_allele_DGE  <- read.table(H9_allele_DGE,
                             header = TRUE,
                                         sep = "\t",
                                         stringsAsFactors = FALSE)

get_list_of_sig_genes <- function(df){
  sig_genes <- df[df$FDR < 0.05 & abs(df$logFC) >= 1, "gene_id"]
  return(sig_genes)
}


# Get lists of significant genes for each comparison
EP1_sig_genes <- get_list_of_sig_genes(EP1_allele_DGE)
H9_sig_genes <- get_list_of_sig_genes(H9_allele_DGE)


# Create the Venn diagram
venn.plot.dir <- here("processed_data", "ASE", "bambu_counts_matrices","venn_plots")
dir.create(venn.plot.dir, recursive = TRUE, showWarnings = FALSE)

venn.diagram(
  x = list(
    EP1 = EP1_sig_genes,
    H9 = H9_sig_genes),
  filename = file.path(venn.plot.dir,
                       "sig_genes_overlap_H1_vs_H2_in_H9_and_EP1.png"),
  category.names = c( "EP1", "H9"),
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  col = "black",
  fill = c("red", "purple"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)


H1_vs_H2_in_H9_and_EP1_genes <- unique(c(
  EP1_sig_genes,
  H9_sig_genes
))

df <- data.frame(gene = H1_vs_H2_in_H9_and_EP1_genes)
df$EP1 <- ifelse(df$gene %in%  EP1_sig_genes, "yes", "no")
df$H9 <- ifelse(df$gene %in% H9_sig_genes, "yes", "no")


readr::write_tsv(df,
                 file.path(venn.plot.dir,
                           "sig_genes_overlap_H1_vs_H2_in_H9_and_EP1.tsv"))


