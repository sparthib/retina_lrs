library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)


### GO Analysis of upregulated genes or genes associated to upregulated 
### DTE.

## or genes undergoing significant switching. 

get_dge_genelist <- function(df) {
  values <- df |> select(DGE_log2FC, DGE, gene_id) |> distinct() |>
    filter(DGE == TRUE) |> filter(DGE_log2FC > 0) |>
    pull(DGE_log2FC) |> as.vector()
  
  names <- df |> select(DGE_log2FC, DGE, gene_id) |> distinct() |>
    filter(DGE == TRUE) |> filter(DGE_log2FC > 0) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE) 
}


get_dtu_genelist <- function(df) {
  values <- df|> select(dIF, DTU, gene_id) |> distinct() |>
    filter(DTU) |> 
    pull(dIF) |> as.vector()
  names <- df |> select(dIF, DTU, gene_id) |> distinct() |>
    filter(DTU) |> 
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE) 
}


get_dte_genelist <- function(df) {
  values <- df |> select(DTE_log2FC, DTE, gene_id) |> distinct() |>
    filter(DTE) |> filter(DTE_log2FC > 0) |>
    pull(DTE_log2FC) |> as.vector()
  names <- df |> select(DTE_log2FC, DTE, gene_id) |> distinct() |>
    filter(DTE) |> filter(DTE_log2FC > 0) |>
    pull(gene_id) |> as.vector()
  names(values) <- names
  #remove version number in names
  names <- gsub("\\..*", "", names)
  #sort values in decreasing order using sort function
  values <- sort(values, decreasing = TRUE) 
}


ora_plot <- function(genelist, ont, output_plot_dir, analysis_type, conditions){

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
  if(nrow(as.data.frame(ego)) != 0){
    write_tsv(as.data.frame(ego), file.path(output_plot_dir, paste0(conditions, analysis_type,"_ora_", ont, ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0(conditions, analysis_type,"_ora_", ont, ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  }
}

  
run_all_go <- function(method, comparison, ont = "BP"){ 
  DGE_DTU_DTE <- read_tsv(file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                                    method, comparison,
                                    "DGE_DTE_DTU.tsv"))
  
  output_plot_dir <- file.path("/users", "sparthib", "retina_lrs", "processed_data","dtu",
                               method, comparison, "plots", "go_analysis")
  if (!dir.exists(output_plot_dir)) {
    dir.create(output_plot_dir, recursive = TRUE)
  }
  
  if (comparison == "ROs") { 
    D100_vs_D45 <- DGE_DTU_DTE |> filter(condition_1 == "B_RO_D100" & condition_2 == "C_RO_D45")
    D200_vs_D45 <- DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" & condition_2 == "C_RO_D45")
    D200_vs_D100 <-DGE_DTU_DTE |> filter(condition_1 == "A_RO_D200" & condition_2 == "B_RO_D100")
    
    # D100_vs_D45_dge <- get_dge_genelist(D100_vs_D45)
    # D200_vs_D45_dge <- get_dge_genelist(D200_vs_D45)
    # D200_vs_D100_dge <- get_dge_genelist(D200_vs_D100)
    # 
    # 
    # D100_vs_D45_dte <- get_dte_genelist(D100_vs_D45)
    # D200_vs_D45_dte <- get_dte_genelist(D200_vs_D45)
    # D200_vs_D100_dte <- get_dte_genelist(D200_vs_D100)
    
    D100_vs_D45_dtu <- get_dtu_genelist(D100_vs_D45)
    D200_vs_D45_dtu <- get_dtu_genelist(D200_vs_D45)
    D200_vs_D100_dtu <- get_dtu_genelist(D200_vs_D100)
    
    
    # ora_plot(D100_vs_D45_dge, ont, output_plot_dir, "DGE", "D100_vs_D45")
    # ora_plot(D200_vs_D45_dge, ont, output_plot_dir, "DGE", "D200_vs_D45")
    # ora_plot(D200_vs_D100_dge, ont, output_plot_dir, "DGE", "D200_vs_D100")
    # 
    # 
    # 
    # ora_plot(D100_vs_D45_dte, ont, output_plot_dir, "DTE", "D100_vs_D45")
    # ora_plot(D200_vs_D45_dte, ont, output_plot_dir, "DTE", "D200_vs_D45")
    # ora_plot(D200_vs_D100_dte, ont, output_plot_dir, "DTE", "D200_vs_D100")
    
    ora_plot(D100_vs_D45_dtu, ont, output_plot_dir, "DTU", "D100_vs_D45")
    ora_plot(D200_vs_D45_dtu, ont, output_plot_dir, "DTU", "D200_vs_D45")
    ora_plot(D200_vs_D100_dtu, ont, output_plot_dir, "DTU", "D200_vs_D100")
    
  } else if (comparison == "FT_vs_RGC") {
    
    # dge <- get_dge_genelist(DGE_DTU_DTE)
    # 
    # dte <- get_dte_genelist(DGE_DTU_DTE)
    
    dtu <- get_dtu_genelist(DGE_DTU_DTE)
    
    ora_plot(dge, ont, output_plot_dir, "DGE", "FT_vs_RGC")
    ora_plot(dte, ont, output_plot_dir, "DTE", "FT_vs_RGC")
    
    # ora_plot(dtu, ont, output_plot_dir, "DTU", "FT_vs_RGC")
  
  }
  
  
  }


methods <- c("bambu", "Isoquant")
comparisons <- c("ROs", "FT_vs_RGC")
ontologies <- c("BP", "CC", "MF")

for (method in methods) {
  for (comparison in comparisons) {
    for (ont in ontologies) {
      run_all_go(method, comparison, ont)
    }
  }
}




