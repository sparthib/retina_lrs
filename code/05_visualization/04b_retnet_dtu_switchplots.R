library(IsoformSwitchAnalyzeR)

#list of DTU genes to plot switch plots for 

retnet_dtu_transcripts <- c(
  "ENST00000220676", "ENST00000233596", "ENST00000245157",
  "ENST00000307340", "ENST00000366998", "ENST00000375735",
  "ENST00000380656", "ENST00000395479", "ENST00000436393",
  "ENST00000447510", "ENST00000461727", "ENST00000503581",
  "ENST00000540805", "ENST00000553885", "ENST00000562520",
  "ENST00000563137", "ENST00000619831", "ENST00000620466",
  "ENST00000622513", "ENST00000636932", "ENST00000642395",
  "ENST00000676363", "ENST00000676883"
)

switchlist_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/"
comparison <- "ROs"
switchlist_dir <- file.path(switchlist_dir, comparison,  "rds", "SwitchList_part2.rds")

switchAnalysisObject <- readRDS(switchlist_dir)


library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Get gene names for the transcripts
gene_ids <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
  filters = "ensembl_transcript_id",
  values = retnet_dtu_transcripts,
  mart = mart
)
gene_ids <- gene_ids |> 
  dplyr::pull(ensembl_gene_id) |> unique()
condition_pairs <- list(c("Stage_1", "Stage_2"), 
                        c("Stage_1", "Stage_3"),
                        c("Stage_2", "Stage_3"))

## SwitchAnalysis object 
retnet_plot_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/bambu/ROs/protein_coding/plots/retnet"
pdf(file = file.path(retnet_plot_dir, "switch_plots.pdf"), width = 8, height = 6)

for(condition_vec in condition_pairs){ 
  for (gene in gene_ids) {
    switchPlot(
      ### Core arguments
      switchAnalysisObject,
      gene = gene,
      condition1 = condition_vec[1],
      condition2 = condition_vec[2],
    )
  }
  }
  
dev.off()



