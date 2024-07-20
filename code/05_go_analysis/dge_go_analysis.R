library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

# Load bambu DGE data

FT_vs_RGC_bambu_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","bambu", "FT_vs_RGC",
                                    "DGEs.tsv"))
D100_vs_D45_bambu_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","bambu", "ROs",
                                    "D100_vs_D45_DGEs.tsv"))
D200_vs_D45_bambu_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","bambu", "ROs",
                                    "D200_vs_D45_DGEs.tsv"))
D100_vs_D200_bambu_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","bambu", "ROs",
                                    "D200_vs_D100_DGEs.tsv"))

# Load isoquant DGE data
FT_vs_RGC_isoquant_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","isoquant", "FT_vs_RGC",
                                    "DGEs.tsv"))
D100_vs_D45_isoquant_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","isoquant", "ROs",
                                    "D100_vs_D45_DGEs.tsv"))
D200_vs_D45_isoquant_dge  <- read_tsv(here("processed_data","dge",
                                    "edgeR","isoquant", "ROs",
                                    "D200_vs_D45_DGEs.tsv"))
D100_vs_D200_isoquant_dge <- read_tsv(here("processed_data","dge",
                                    "edgeR","isoquant", "ROs",
                                    "D200_vs_D100_DGEs.tsv"))


# get geneList function 

getGeneList <- function(dge, fdr = 0.05, log2fc_cutoff = 1) {
  values <- dge |>
    filter(FDR < fdr, abs(logFC) > log2fc_cutoff) |>
    pull(logFC) |> as.vector()
  names <- dge |> 
    filter(FDR < fdr, abs(logFC) > log2fc_cutoff) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE)
}

# get geneList for bambu DGE data
FT_vs_RGC_bambu_geneList <- getGeneList(FT_vs_RGC_bambu_dge)
length(FT_vs_RGC_bambu_geneList)
#1645
D100_vs_D45_bambu_geneList <- getGeneList(D100_vs_D45_bambu_dge)
length(D100_vs_D45_bambu_geneList)
#847
D200_vs_D45_bambu_geneList <- getGeneList(D200_vs_D45_bambu_dge)
length(D200_vs_D45_bambu_geneList)
#4629
D100_vs_D200_bambu_geneList <- getGeneList(D100_vs_D200_bambu_dge)
length(D100_vs_D200_bambu_geneList)
#1289

# get geneList for isoquant DGE data
FT_vs_RGC_isoquant_geneList <- getGeneList(FT_vs_RGC_isoquant_dge)
length(FT_vs_RGC_isoquant_geneList)
# 2449
D100_vs_D45_isoquant_geneList <- getGeneList(D100_vs_D45_isoquant_dge)
length(D100_vs_D45_isoquant_geneList)
# 1150
D200_vs_D45_isoquant_geneList <- getGeneList(D200_vs_D45_isoquant_dge)
length(D200_vs_D45_isoquant_geneList)
# 2326
D100_vs_D200_isoquant_geneList <- getGeneList(D100_vs_D200_isoquant_dge)
length(D100_vs_D200_isoquant_geneList)
# 1175


########make plots for over representation analysis ########
ora_data <- function(genelist, ont, output_data_dir){
  ego <- enrichGO(gene          = names(genelist),
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = ont,
                  pAdjustMethod = "fdr",
                  minGSSize     = 100,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01,
                  readable      = TRUE) 
  write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("DGE_ora_", ont, ".tsv")))
}

ora_data(FT_vs_RGC_bambu_geneList, "MF", here("processed_data", "go",
                                              "bambu", "FT_vs_RGC"))
ora_data(D100_vs_D45_bambu_geneList, "MF", here("processed_data", "go",
                                                "bambu", "RO_D100_vs_D45"))
ora_data(D200_vs_D45_bambu_geneList, "MF", here("processed_data", "go",
                                                "bambu", "RO_D200_vs_D45"))
ora_data(D100_vs_D200_bambu_geneList, "MF", here("processed_data", "go",
                                                 "bambu", "RO_D100_vs_D200"))

ora_data(FT_vs_RGC_bambu_geneList, "BP", here("processed_data", "go",
                                              "bambu", "FT_vs_RGC"))
ora_data(D100_vs_D45_bambu_geneList, "BP", here("processed_data", "go",
                                                "bambu", "RO_D100_vs_D45"))
ora_data(D200_vs_D45_bambu_geneList, "BP", here("processed_data", "go",
                                                "bambu", "RO_D200_vs_D45"))
ora_data(D100_vs_D200_bambu_geneList, "BP", here("processed_data", "go",
                                                 "bambu", "RO_D100_vs_D200"))


ora_data(FT_vs_RGC_bambu_geneList, "CC", here("processed_data", "go",
                                              "bambu", "FT_vs_RGC"))
ora_data(D100_vs_D45_bambu_geneList, "CC", here("processed_data", "go",
                                                "bambu", "RO_D100_vs_D45"))
ora_data(D200_vs_D45_bambu_geneList, "CC", here("processed_data", "go",
                                                "bambu", "RO_D200_vs_D45"))
ora_data(D100_vs_D200_bambu_geneList, "CC", here("processed_data", "go",
                                                 "bambu", "RO_D100_vs_D200"))
# 
# ora_data(FT_vs_RGC_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "FT_vs_RGC"))
# ora_data(D100_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D45"))
# ora_data(D200_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D45"))
# ora_data(D100_vs_D200_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D200"))
#### ####
#### ####
ora_plots <- function(genelist, ont, output_plot_dir){
  ego <- enrichGO(gene          = names(genelist),
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = ont,
                  pAdjustMethod = "fdr",
                  minGSSize     = 100,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01,
                  readable      = TRUE) 
  # write_tsv(ego, file.path(output_data_dir, paste0("ora_", ont, ".tsv")))
  #dotplot
  pdf(file.path(output_plot_dir, paste0("DGE_ora_dotplot_", ont, ".pdf")))
  print(dotplot(ego, showCategory = 15))
  dev.off()
  
  #cnetplot
  # pdf(file.path(output_plot_dir, paste0("ora_cnetplot_", ont, ".pdf")))
  # print(cnetplot(ego, foldChange=genelist, colorEdge = TRUE))
  # dev.off()
  # 
  # #emapplot
  # ego <- enrichplot::pairwise_termsim(ego)
  # pdf(file.path(output_plot_dir, paste0("ora_emapplot_", ont, ".pdf")))
  # print(emapplot(ego, layout="kk", showCategory=15))
  # dev.off()
  # 
  #goplot
  pdf(file.path(output_plot_dir, paste0("DGE_ora_goplot_", ont, ".pdf")))
  print(goplot(ego, showCategory = 8))
  dev.off()

}

ora_plots(FT_vs_RGC_bambu_geneList, "MF", here("plots", "go_analysis",
                                               "bambu", "FT_vs_RGC"))
ora_plots(D100_vs_D45_bambu_geneList, "MF", here("plots", "go_analysis",
                                                 "bambu", "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_bambu_geneList, "MF", here("plots", "go_analysis",
                                                 "bambu", "RO_D200_vs_D45"))
ora_plots(D100_vs_D200_bambu_geneList, "MF", here("plots", "go_analysis",
                                                  "bambu", "RO_D100_vs_D200"))

ora_plots(FT_vs_RGC_bambu_geneList, "BP", here("plots", "go_analysis",
                                               "bambu", "FT_vs_RGC"))
ora_plots(D100_vs_D45_bambu_geneList, "BP", here("plots", "go_analysis",
                                                 "bambu", "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_bambu_geneList, "BP", here("plots", "go_analysis",
                                                 "bambu", "RO_D200_vs_D45"))
ora_plots(D100_vs_D200_bambu_geneList, "BP", here("plots", "go_analysis",
                                                  "bambu", "RO_D100_vs_D200"))

ora_plots(FT_vs_RGC_bambu_geneList, "CC", here("plots", "go_analysis",
                                               "bambu", "FT_vs_RGC"))
ora_plots(D100_vs_D45_bambu_geneList, "CC", here("plots", "go_analysis",
                                                 "bambu", "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_bambu_geneList, "CC", here("plots", "go_analysis",
                                                 "bambu", "RO_D200_vs_D45"))
ora_plots(D100_vs_D200_bambu_geneList, "CC", here("plots", "go_analysis",
                                                  "bambu", "RO_D100_vs_D200"))

# 
# 
# ora_plots(FT_vs_RGC_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "FT_vs_RGC"))
# ora_plots(D100_vs_D45_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D100_vs_D45"))
# ora_plots(D200_vs_D45_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D200_vs_D45"))
# ora_plots(D100_vs_D200_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D100_vs_D200"))
# #### save gseGO data ####
# 
# gse_data <- function(genelist, ont, output_data_dir){
#   ego <- gseGO(geneList     = genelist,
#                 OrgDb        = org.Hs.eg.db,
#                 keyType  = "ENSEMBL",
#                 ont          = ont,
#                pAdjustMethod = "fdr",
#                minGSSize     = 100,
#                pvalueCutoff  = 0.001,
#                 eps = 0,
#                 verbose = FALSE)
#   ego <- setReadable(ego, org.Hs.eg.db, keyType = "ENSEMBL")
#   write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("gse_", ont, ".tsv")))
# }
# 
# gse_data(FT_vs_RGC_bambu_geneList, "MF", here("processed_data", "go",
#                                                         "bambu", "FT_vs_RGC"))
# gse_data(D100_vs_D45_bambu_geneList, "MF", here("processed_data", "go",
#                                                         "bambu", "RO_D100_vs_D45"))
# gse_data(D200_vs_D45_bambu_geneList, "MF", here("processed_data", "go",
#                                                         "bambu", "RO_D200_vs_D45"))
# gse_data(D100_vs_D200_bambu_geneList, "MF", here("processed_data", "go",
#                                                         "bambu", "RO_D100_vs_D200"))
# 
# gse_data(FT_vs_RGC_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "FT_vs_RGC"))
# gse_data(D100_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D45"))
# gse_data(D200_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D45"))
# gse_data(D100_vs_D200_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D200"))
# 
# #make plots for enrichment analysis 
# 
# gsea_plots <- function(genelist, ont, output_plot_dir){
#   ego <- gseGO(geneList     = genelist,
#                 OrgDb        = org.Hs.eg.db,
#                 keyType  = "ENSEMBL",
#                 ont          = ont,
#                pAdjustMethod = "fdr",
#                minGSSize     = 100,
#                pvalueCutoff  = 0.001,
#                 eps = 0,
#                 verbose = FALSE)
#   ego <- setReadable(ego, org.Hs.eg.db, keyType = "ENSEMBL")
#   #dotplot
#   pdf(file.path(output_plot_dir, paste0("gsea_dotplot_", ont, ".pdf")))
#   print(dotplot(ego, showCategory = 15))
#   dev.off()
#   
#   #cnetplot
#   pdf(file.path(output_plot_dir, paste0("gsea_cnetplot_", ont, ".pdf")))
#   print(cnetplot(ego, foldChange=genelist, colorEdge = TRUE))
#   dev.off()
#   
#   #emapplot
#   ego <- enrichplot::pairwise_termsim(ego)
#   pdf(file.path(output_plot_dir, paste0("gsea_emapplot_", ont, ".pdf")))
#   print(emapplot(ego, layout="kk", showCategory=15))
#   dev.off()
#   
#   #ridgeplot
#   pdf(file.path(output_plot_dir, paste0("gsea_ridgeplot_", ont, ".pdf")))
#   print(enrichplot::ridgeplot(ego))
#   dev.off()
#   
# }
# 
# gsea_plots(FT_vs_RGC_bambu_geneList, "MF", here("plots", "go_analysis",
#                                                         "bambu", "FT_vs_RGC"))
# gsea_plots(D100_vs_D45_bambu_geneList, "MF", here("plots", "go_analysis",
#                                                         "bambu", "RO_D100_vs_D45"))
# gsea_plots(D200_vs_D45_bambu_geneList, "MF", here("plots", "go_analysis",
#                                                         "bambu", "RO_D200_vs_D45"))
# gsea_plots(D100_vs_D200_bambu_geneList, "MF", here("plots", "go_analysis",
#                                                         "bambu", "RO_D100_vs_D200"))
# 
# gsea_plots(FT_vs_RGC_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "FT_vs_RGC"))
# gsea_plots(D100_vs_D45_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D100_vs_D45"))
# gsea_plots(D200_vs_D45_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D200_vs_D45"))
# gsea_plots(D100_vs_D200_isoquant_geneList, "MF", here("plots", "go_analysis",
#                                                         "isoquant", "RO_D100_vs_D200"))
# 



