library(IsoformSwitchAnalyzeR)
library(biomaRt)


switchlist_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/"
comparison <- "ROs"
switchlist_dir <- file.path(switchlist_dir, comparison, "protein_coding", "rds", "SwitchList_part2.rds")

switchAnalysisObject <- readRDS(switchlist_dir)
DGE_DTU_DTE_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/protein_coding/"
DGE_DTU_DTE <- read_tsv(file.path(DGE_DTU_DTE_dir, "DGE_DTE_DTU.tsv"))


gene_names <- c("RP1", "CRX", "PROM1", "CRB1")

condition_pairs <- list(c("Stage_1", "Stage_2"), 
                        c("Stage_1", "Stage_3"),
                        c("Stage_2", "Stage_3"))

gene_ids <- DGE_DTU_DTE |>
  filter(gene_name %in% gene_names) |>
  pull(gene_id) |>
  unique()

most_switched_genes <- extractTopSwitches(
  switchAnalysisObject,
  filterForConsequences = TRUE,
  n = 10000,
  extractGenes=TRUE
)

retnet_switched_genes  <- most_switched_genes |> filter(gene_name %in% gene_names)

most_switched_isoforms <- extractTopSwitches(
  switchAnalysisObject,
  filterForConsequences = TRUE,
  n = 40000,
  extractGenes=FALSE
)
retnet_switched_isoforms <- most_switched_isoforms |> filter(gene_name %in% gene_names)
retnet_switched_isoforms$isoform_id
# [1] "ENST00000636932" "ENST00000636932" "ENST00000220676" "ENST00000540805"
# [5] "ENST00000447510" "ENST00000539194" "ENST00000508167" "ENST00000540805"
# [9] "ENST00000220676" "ENST00000447510" "ENST00000613299" "ENST00000508167"
# [13] "ENST00000221996" "ENST00000505450" "ENST00000505450" "ENST00000367397"
# [17] "ENST00000513448" "ENST00000538660"

# print as character vector 

retnet_plot_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/protein_coding/plots/retnet"
pdf(file = file.path(retnet_plot_dir, "switch_plots_with_consequences.pdf"), 
    width = 8, height = 6)

# Loop through each gene
for (i in 1:nrow(retnet_switched_genes)) {
  switchPlot(
    switchAnalyzeRlist = switchAnalysisObject,
    gene = retnet_switched_genes$gene_id[i],
    condition1 = retnet_switched_genes$condition_1[i],
    condition2 = retnet_switched_genes$condition_2[i]
  )
}

# Close PDF device
dev.off()


