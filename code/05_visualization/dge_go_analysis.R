library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

# Load DGE data

D100_vs_D45_dge <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                      "Isoquant", "ROs",
                                      "DGE_DTU_DTE.tsv"))
D100_vs_D45_dge <- D100_vs_D45_dtu |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")

D200_vs_D45_dge <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                      "Isoquant", "ROs",
                                      "DGE_DTU_DTE.tsv"))
D200_vs_D45_dge <- D200_vs_D45_dtu |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")

D200_vs_D100_dge <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                       "Isoquant", "ROs",
                                       "DGE_DTU_DTE.tsv"))
D200_vs_D100_dge <- D200_vs_D100_dtu |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")


# get geneList function 

dge=D200_vs_D100_dtu
fdr = 0.05
log2fc_cutoff = 1

getGeneList <- function(dge, fdr = 0.05, log2fc_cutoff = 1) {
  values <- dge |> select(DGE_log2FC, DGE_qval, gene_id) |> distinct() |>
    filter(DGE_qval < fdr, abs(DGE_log2FC) >= log2fc_cutoff) |>
    pull(DGE_log2FC) |> as.vector()
  names <- dge |> select(DGE_log2FC, DGE_qval, gene_id) |> distinct() |>
    filter(DGE_qval < fdr, abs(DGE_log2FC) >= log2fc_cutoff) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE) 
}

# get geneList for bambu DGE data
# FT_vs_RGC_geneList <- getGeneList(FT_vs_RGC_dge)
# length(FT_vs_RGC_geneList)
#1645

D100_vs_D45_geneList <- getGeneList(D100_vs_D45_dge)
length(D100_vs_D45_geneList)

D200_vs_D45_geneList <- getGeneList(D200_vs_D45_dge)
length(D200_vs_D45_geneList)

D200_vs_D100_geneList <- getGeneList(D200_vs_D100_dge)
length(D200_vs_D100_geneList)



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
  if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }
  write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("DGE_ora_", ont, ".tsv")))
}

# ora_data(FT_vs_RGC_geneList, "MF", here("processed_data", "go",
#                                                "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "MF",file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                  "RO_D200_vs_D100"))

# ora_data(FT_vs_RGC_geneList, "BP", here("processed_data", "go",
#                                                "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "BP",file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                  "RO_D200_vs_D100"))

# 
# ora_data(FT_vs_RGC_geneList, "CC", here("processed_data", "go",
#                                                "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                  "RO_D200_vs_D100"))
# 
# ora_data(FT_vs_RGC_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "FT_vs_RGC"))
# ora_data(D100_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D45"))
# ora_data(D200_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D45"))
# ora_data(D200_vs_D100_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D100"))
#### ####
#### ####

genelist = D100_vs_D45_geneList
ont = "MF"
output_data_dir = file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45")
output_plot_dir = file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D100_vs_D45")

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
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  pdf(file.path(output_plot_dir, paste0("DGE_ora_dotplot_", ont, ".pdf")))
  print(dotplot(ego, showCategory = 15))
  dev.off()
  
 

}


ora_plots(FT_vs_RGC_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                              "FT_vs_RGC"))
ora_plots(D100_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D200_vs_D45"))
ora_plots(D200_vs_D100_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                 "RO_D200_vs_D100"))

ora_plots(FT_vs_RGC_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                              "FT_vs_RGC"))
ora_plots(D100_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D200_vs_D45"))
ora_plots(D200_vs_D100_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                 "RO_D200_vs_D100"))

ora_plots(FT_vs_RGC_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                              "FT_vs_RGC"))
ora_plots(D100_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "RO_D200_vs_D45"))
ora_plots(D200_vs_D100_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                 "RO_D200_vs_D100"))



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
# pdf(file.path(output_plot_dir, paste0("DGE_ora_goplot_", ont, ".pdf")))
# print(goplot(ego, showCategory = 8))
# dev.off()




# 
# 
# ora_plots(FT_vs_RGC_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "FT_vs_RGC"))
# ora_plots(D100_vs_D45_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D100_vs_D45"))
# ora_plots(D200_vs_D45_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D200_vs_D45"))
# ora_plots(D200_vs_D100_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D200_vs_D100"))
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
# gse_data(FT_vs_RGC_geneList, "MF", here("processed_data", "go",
#                                                          "FT_vs_RGC"))
# gse_data(D100_vs_D45_geneList, "MF", here("processed_data", "go",
#                                                          "RO_D100_vs_D45"))
# gse_data(D200_vs_D45_geneList, "MF", here("processed_data", "go",
#                                                          "RO_D200_vs_D45"))
# gse_data(D200_vs_D100_geneList, "MF", here("processed_data", "go",
#                                                          "RO_D200_vs_D100"))
# 
# gse_data(FT_vs_RGC_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "FT_vs_RGC"))
# gse_data(D100_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D100_vs_D45"))
# gse_data(D200_vs_D45_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D45"))
# gse_data(D200_vs_D100_isoquant_geneList, "MF", here("processed_data", "go",
#                                                         "isoquant", "RO_D200_vs_D100"))
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
# gsea_plots(FT_vs_RGC_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                          "FT_vs_RGC"))
# gsea_plots(D100_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                          "RO_D100_vs_D45"))
# gsea_plots(D200_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                          "RO_D200_vs_D45"))
# gsea_plots(D200_vs_D100_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                          "RO_D200_vs_D100"))
# 
# gsea_plots(FT_vs_RGC_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "FT_vs_RGC"))
# gsea_plots(D100_vs_D45_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D100_vs_D45"))
# gsea_plots(D200_vs_D45_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D200_vs_D45"))
# gsea_plots(D200_vs_D100_isoquant_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
#                                                         "isoquant", "RO_D200_vs_D100"))
# 



