library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

# Load Isoquant DGE data


D100_vs_D45_dtu <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                       "Isoquant", "ROs",
                                       "DGE_DTU_DTE.tsv"))
D100_vs_D45_dtu <- D100_vs_D45_dtu |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")

D200_vs_D45_dtu <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                       "Isoquant", "ROs",
                                       "DGE_DTU_DTE.tsv"))
D200_vs_D45_dtu <- D200_vs_D45_dtu |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")

D200_vs_D100_dtu <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                        "Isoquant", "ROs",
                                        "DGE_DTU_DTE.tsv"))
D200_vs_D100_dtu <- D200_vs_D100_dtu |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")


# get geneList function 

getGeneList <- function(dtu, fdr = 0.05, dIF = 0.1) {
  values <- dtu |>
    filter(DTU_qval < 0.05, abs(DTU_dIF) >= 0.1) |>
    pull(DTU_dIF) |> as.vector()
  names <- dtu |> 
    filter(DTU_qval < 0.05, abs(DTU_dIF) >= 0.1) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*",  names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE)
  
}

# get geneList for  DTU data
# FT_vs_RGC_geneList <- getGeneList(FT_vs_RGC_dtu)
# length(unique(names(FT_vs_RGC_geneList)))

D100_vs_D45_geneList <- getGeneList(D100_vs_D45_dtu)
length(unique(names(D100_vs_D45_geneList)))

D200_vs_D45_geneList <- getGeneList(D200_vs_D45_dtu)
length(unique(names(D200_vs_D45_geneList)))

D200_vs_D100_geneList <- getGeneList(D200_vs_D100_dtu)
length(unique(names(D200_vs_D100_geneList)))



########make plots for over representation analysis ########
ora_data <- function(genelist, ont, output_data_dir){
  ego <- enrichGO(gene          = unique(names(genelist)),
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = ont,
                  pAdjustMethod = "fdr",
                  minGSSize     = 100,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01,
                  readable      = TRUE) 
  if (!dir.exists(output_data_dir)){
    dir.create(output_data_dir, recursive = TRUE)
  }
  write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("DTU_ora_", ont, ".tsv")))
}

# dir.create(file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go"),
#            recursive = TRUE, showWarnings = FALSE)

ora_data(FT_vs_RGC_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                               "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "MF", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D100"))

ora_data(FT_vs_RGC_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                               "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "BP", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                  "RO_D200_vs_D100"))


ora_data(FT_vs_RGC_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                               "FT_vs_RGC"))
ora_data(D100_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D100_vs_D45"))
ora_data(D200_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                 "RO_D200_vs_D45"))
ora_data(D200_vs_D100_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "processed_data", "isoquant","go",
                                                  "RO_D200_vs_D100"))

#### ####
ora_plots <- function(genelist, ont, output_plot_dir){
  ego <- enrichGO(gene          = unique(names(genelist)),
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
  
  if (!dir.exists(output_plot_dir)){
    dir.create(output_plot_dir, recursive = TRUE)
  }
  pdf(file.path(output_plot_dir, paste0("DTU_ora_dotplot_", ont, ".pdf")))
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
ora_plots(D200_vs_D100_geneList, "BP",file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                   "RO_D200_vs_D100"))

ora_plots(FT_vs_RGC_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                "FT_vs_RGC"))

ora_plots(D100_vs_D45_geneList, "CC",file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                  "RO_D100_vs_D45"))
ora_plots(D200_vs_D45_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                  "RO_D200_vs_D45"))
ora_plots(D200_vs_D100_geneList, "CC", file.path("/users", "sparthib", "retina_lrs", "plots", "go_analysis", "isoquant",
                                                   "RO_D200_vs_D100"))

