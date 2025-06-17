# 166 genes with DTU in alll 3
#comparisons (from the upset plot) --
#GO analysis and see if these 166 genes
#are enriched in any particular biological process(es)?

library(clusterProfiler)
library(org.Hs.eg.db)
library(here)
library(readr)
library(dplyr)
library(tidyr)

method <- "bambu"
comparison <- "ROs"

input_data_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison, "protein_coding")
plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison,"protein_coding", "plots", "upset")

gene_overlaps <- readr::read_tsv( file.path(plots_dir,
                                               "gene_overlaps_DTU_genes.tsv"))


all_stage_genes <- gene_overlaps |> dplyr::filter(DTU_Stage_1_vs_Stage_2 == TRUE & 
                                             DTU_Stage_1_vs_Stage_3 == TRUE & 
                                             DTU_Stage_2_vs_Stage_3 == TRUE)

## GO analysis 

DGE_DTU_DTE <- read_tsv(file.path(input_data_dir, "DGE_DTE_DTU.tsv"))



# get gene_id, DGE_log2FC

get_dtu_genelist <- function(data) {
  data |> 
    dplyr::select(gene_id, DGE_log2FC) |> 
    dplyr::mutate(gene_id = gsub("\\..*", "", gene_id)) |> 
    dplyr::filter(gene_id %in% all_stage_genes$gene_id) |> distinct()
  
  # for each gene ID keep the max log2FC value
  dtu_genelist <- dtu_genelist |> 
    dplyr::group_by(gene_id) |> 
    dplyr::summarise(DGE_log2FC = max(DGE_log2FC, na.rm = TRUE)) |> 
    dplyr::ungroup()
  
  names <- dtu_genelist |> pull(gene_id)
  values <- dtu_genelist |> pull(DGE_log2FC) |> sort(decreasing = TRUE)
  names(values) <- names
  
  return(values)
}

dtu_genelist <- get_dtu_genelist(DGE_DTU_DTE)

# Convert gene IDs to 

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
  if(nrow(as.data.frame(ego)) != 0){
    write_tsv(as.data.frame(ego), file.path(output_plot_dir,
                                            paste0("ORA_all_stage_DTU_genes_", ont, ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("ORA_all_stage_DTU_genes_", ont, ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  }
}

ora_plot(dtu_genelist, "BP", plots_dir, "DTU")

