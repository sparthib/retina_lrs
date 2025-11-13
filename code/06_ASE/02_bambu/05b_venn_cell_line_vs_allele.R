library(here)
library(VennDiagram)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db)
library(readr)
library(clusterProfiler)

code_dir <- Sys.getenv("retina_lrs_code")
data_dir <- Sys.getenv("retina_lrs_dir")

EP1_allele_DGE <- file.path(code_dir, "processed_data/ASE/bambu_counts_matrices/DGE/cell_line_vs_allele/H1_EP1_vs_H2_EP1_DGEs.tsv")
H9_allele_DGE <- file.path(code_dir, "processed_data/ASE/bambu_counts_matrices/DGE/cell_line_vs_allele/H1_H9_vs_H2_H9_DGEs.tsv")

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
venn.plot.dir <- file.path(code_dir, "processed_data", "ASE", "bambu_counts_matrices","venn_plots")
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

##how many EP1_ASE genes are IDR? 
sig_genes <- readr::read_tsv(file.path(venn.plot.dir,
                           "sig_genes_overlap_H1_vs_H2_in_H9_and_EP1.tsv")) 

ird_genes <- readr::read_tsv(file = file.path(code_dir, "processed_data/dtu/retnet_disease_genes.tsv"), 
                             col_names = TRUE)


# get ensembl names for ird_genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ird_genes_annot <- getBM(attributes = c( 'ensembl_gene_id', 'hgnc_symbol', "chromosome_name"),
                           filters = 'hgnc_symbol',
                           values = ird_genes$gene_name,
                           mart = ensembl)

ird_genes$ensembl_gene_id <- ird_genes_annot$ensembl_gene_id[match(ird_genes$gene_name,
                                                             ird_genes_annot$hgnc_symbol)]
ird_genes$chromosome <- ird_genes_annot$chromosome_name[match(ird_genes$gene_name,
                                                             ird_genes_annot$hgnc_symbol)]

readr::write_tsv(ird_genes,
                 file.path(code_dir, "processed_data/dtu/retnet_disease_genes.tsv"))


# overlaps between EP1 ASE ge
EP1_sig_genes <- sig_genes |> dplyr::filter(EP1 == "yes") |> dplyr::pull(gene)

# how many EP1_sig_genes are in ird_genes$ensembl_gene_id?
overlap_EP1_ird <- EP1_sig_genes[EP1_sig_genes %in% ird_genes$ensembl_gene_id]

EP1_IRD_df <- ird_genes[ird_genes$ensembl_gene_id %in% overlap_EP1_ird, ]

write.table(EP1_IRD_df,
            file = file.path(venn.plot.dir,
                             "EP1_ASE_IDR_genes.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


## GO plot for the shared 168 genes 

shared_genes <- sig_genes |> 
  dplyr::filter(EP1 == "yes" & H9 == "yes") |>
  dplyr::pull(gene)


go_plot_dir <- file.path(code_dir, "processed_data", "ASE", "bambu_counts_matrices","GO_plots")


ora_plot <- function(genelist, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = genelist,
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  readable      = TRUE) 
  if(nrow(as.data.frame(ego)) != 0){
    #ego <- enrichplot::pairwise_termsim(ego)
    # ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    write_tsv(as.data.frame(ego), file.path(output_plot_dir,
                                            paste0("simplified_ORA_ASE_DGE_genes_",analysis_type, "BP", ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("simplified_ORA_all_ASE_DGE_genes_",analysis_type, "BP", ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  } else {
    message("No significant GO terms found for ", analysis_type, " in ", "BP")}
}


ora_plot(shared_genes,
         go_plot_dir,
         "sig_genes_in_H9_and_EP1")

EP1_only_sig_genes <- sig_genes |>
  dplyr::filter(EP1 == "yes" & H9 == "no") |> 
  dplyr::pull(gene)

ora_plot(EP1_sig_genes,
         go_plot_dir,
         "sig_genes_in_EP1_only")



