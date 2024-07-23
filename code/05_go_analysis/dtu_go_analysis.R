library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

# Load bambu DGE data

FT_vs_RGC_bambu_dtu <- read_tsv(here("processed_data","dtu",
                                     "DTU_gandall","bambu", "FT_vs_RGC",
                                     "DGE_DTU_DTE.tsv"))

D100_vs_D45_bambu_dtu <- read_tsv(here("processed_data","dtu",
                                       "DTU_gandall","bambu", "ROs",
                                       "DGE_DTU_DTE.tsv"))
D100_vs_D45_bambu_dtu <- D100_vs_D45_bambu_dtu |> filter(condition_1 == "RO_D100" & condition_2 == "RO_D45")

D200_vs_D45_bambu_dtu <- read_tsv(here("processed_data","dtu",
                                       "DTU_gandall","bambu", "ROs",
                                       "DGE_DTU_DTE.tsv"))
D200_vs_D45_bambu_dtu <- D200_vs_D45_bambu_dtu |> filter(condition_1 == "RO_D200" & condition_2 == "RO_D45")
D100_vs_D200_bambu_dtu <- read_tsv(here("processed_data","dtu",
                                        "DTU_gandall","bambu", "ROs",
                                        "DGE_DTU_DTE.tsv"))
D100_vs_D200_bambu_dtu <- D100_vs_D200_bambu_dtu |> filter(condition_1 == "RO_D100" & condition_2 == "RO_D200")


# get geneList function 

getGeneList <- function(dtu, fdr = 0.05, dIF = 0.1) {
  values <- dtu |>
    filter(DTU_qval < 0.05, abs(DTU_dIF) > 0.1) |>
    pull(DTU_dIF) |> as.vector()
  names <- dtu |> 
    filter(DTU_qval < 0.05, abs(DTU_dIF) > 0.1) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE)
  
}

# get geneList for bambu DGE data
FT_vs_RGC_bambu_geneList <- getGeneList(FT_vs_RGC_bambu_dtu)
length(unique(names(FT_vs_RGC_bambu_geneList)))

D100_vs_D45_bambu_geneList <- getGeneList(D100_vs_D45_bambu_dtu)
length(unique(names(D100_vs_D45_bambu_geneList)))

D200_vs_D45_bambu_geneList <- getGeneList(D200_vs_D45_bambu_dtu)
length(unique(names(D200_vs_D45_bambu_geneList)))


D100_vs_D200_bambu_geneList <- getGeneList(D100_vs_D200_bambu_dtu)
length(unique(names(D100_vs_D200_bambu_geneList)))



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
  write_tsv(as.data.frame(ego), file.path(output_data_dir, paste0("DTU_ora_", ont, ".tsv")))
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

#### ####
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
  pdf(file.path(output_plot_dir, paste0("DTU_ora_dotplot_", ont, ".pdf")))
  print(dotplot(ego, showCategory = 15))
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

