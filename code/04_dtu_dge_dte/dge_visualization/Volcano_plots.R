library(ggplot2)
library(here)


volcano_plot <- function(SwitchList_part1, output_dir, cond1, cond2){ 
  df <- SwitchList_part1$isoformFeatures |> 
    dplyr::filter(condition_1 == cond1 & condition_2 == cond2) |>
    dplyr::select(isoform_id, gene_id, gene_name, condition_1, condition_2, dIF, isoform_switch_q_value) |> 
    na.omit()
  
  # Order and filter by dIF
  SwitchList_part1_top_20 <- df[order(abs(df$isoform_switch_q_value)),]
  SwitchList_part1_top_20 <- SwitchList_part1_top_20[abs(SwitchList_part1_top_20$dIF) > 0.5,]
  if (nrow(SwitchList_part1_top_20) > 20) {
    SwitchList_part1_top_20 <- SwitchList_part1_top_20[1:20,]
  }
  
  pdf(paste0(output_dir, "isoform_volcano_", cond1, "_vs_", cond2, ".pdf"))
  
  p <- ggplot(data=df, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
    geom_point(aes(color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), size=1) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') +
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') +
    scale_color_manual('Significant\nIsoform Switch', values = c('grey', 'skyblue')) +
    labs(x='dIF', y='-Log10 (Isoform Switch Q Value)') +
    theme_bw() +
    geom_text(aes(label=ifelse(gene_name %in% SwitchList_part1_top_20$gene_name & abs(dIF) > 0.5, 
                               gene_name, '')), hjust=1, vjust=1, size=2)
  
  print(p)
  dev.off()
}

dtu_rdata_path  = here("processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/rds/DexSeqDTUDGESwitchList.rds")
SwitchList_part1 <- readRDS(dtu_rdata_path)

volcano_plot(SwitchList_part1, here("processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/"), "FT" , "RGC")

dtu_rdata_path  = here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/DexSeqDTUDGESwitchList.rds")
SwitchList_part1 <- readRDS(dtu_rdata_path)

volcano_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D100" , "RO_D45")
volcano_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D200", "RO_D45")
volcano_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D100", "RO_D200")


####.SWITCH VS GENE CHANGES ####


switch_vs_gene_plot <- function(SwitchList_part1, output_dir, cond1, cond2){ 
  df <- SwitchList_part1$isoformFeatures |> 
    dplyr::filter(condition_1 == cond1 & condition_2 == cond2) |>
    dplyr::select(isoform_id, gene_id, gene_name, condition_1, condition_2, dIF, isoform_switch_q_value, gene_switch_q_value,
                  gene_log2_fold_change) |> na.omit()
  
  SwitchList_part1_top_20 <- df[order(abs(df$isoform_switch_q_value)),]
  SwitchList_part1_top_20 <- SwitchList_part1_top_20[abs(SwitchList_part1_top_20$dIF) > 0.5,]
  if (nrow(SwitchList_part1_top_20) > 20) {
    SwitchList_part1_top_20 <- SwitchList_part1_top_20[1:20,]
  }

  pdf(paste0(output_dir,"switch_vs_degs_", cond1, "_vs_", cond2, ".pdf"))
  s <- ggplot(data=df, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
      size=1
    )  + 
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw() + 
    geom_text(aes(label=ifelse(gene_name %in% SwitchList_part1_top_20$gene_name & abs(dIF) > 0.5, 
                               gene_name,'')),
              hjust=0,vjust=0, size = 2)
  print(s)
  dev.off()
  
}

dtu_rdata_path  = here("processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/rds/DexSeqDTUDGESwitchList.rds")
SwitchList_part1 <- readRDS(dtu_rdata_path)

switch_vs_gene_plot(SwitchList_part1, here("processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/"), "FT" , "RGC")

dtu_rdata_path  = here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/DexSeqDTUDGESwitchList.rds")
SwitchList_part1 <- readRDS(dtu_rdata_path)


switch_vs_gene_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D100" , "RO_D45")
switch_vs_gene_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D200", "RO_D45")
switch_vs_gene_plot(SwitchList_part1,here("processed_data/dtu/DTU_gandall/bambu/ROs/plots/"), "RO_D100", "RO_D200")


microexon_data <- function(SwitchList_part1, output_dir){ 
  exons <- as.data.frame(SwitchList_part1$exons)
  microexons <- exons |> dplyr::filter(width < 27)
  
  #remove version number from gene_id
  microexons$gene_id <- gsub("\\..*", "", microexons$gene_id)
  microexons <- merge(microexons, annotLookup, by="gene_id", all.x=TRUE)
  microexons$gene_name <- microexons$ensembl_gene_name
  #drop ensemble_gene_name
  microexons <- microexons |> dplyr::select(-ensembl_gene_name)
  
  #all microexons found in all genes that are expressed in the samples. 
  write_tsv(microexons, paste0(output_dir, "/microexons/all_micro_exons.tsv"))
  
  DEXSeq <- df |>  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |> 
    dplyr::select(isoform_id, dIF, isoform_switch_q_value,gene_id, gene_name, condition_1, condition_2)
  
  ## all microexons found in genes that have some isoform that showed differential usage across any pairwise comparison
  microexons_DTU_genes <- microexons |> filter(gene_id %in% DEXSeq$gene_id)
  nrow(microexons_DTU_genes)
  
  write_tsv(microexons_DTU_genes, paste0(output_dir, "/microexons/microexons_DTU_genes.tsv"))
  
  ## all microexons found in isoforms  that showed differential usage across any pairwise comparison 
  microexons_DTU_isoforms <- microexons |> filter(isoform_id %in% DEXSeq$isoform_id)
  nrow(microexons_DTU_isoforms)
  
  
  write_tsv(microexons_DTU_isoforms, paste0(output_dir, "/microexons/microexons_DTU_isoforms.tsv"))
  
}

microexon_data(SwitchList_part1_D100_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/")
microexon_data(SwitchList_part1_D200_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/")
microexon_data(SwitchList_part1_D100_D200, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/")


###### MAKE SWITCH PLOTS #######

## D100 vs D45
top_genes_D100_D45 <- extractTopSwitches(SwitchList_part1_D100_D45, n=500, 
                                         alpha = 0.05,
                                         dIFcutoff = 0.1)
output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/"
write_tsv(top_genes_D100_D45,
          file = paste0(output_data_dir,  "top_genes.tsv"))

#count number of nas in gene_name
sum(is.na(top_genes_D100_D45$gene_name))


output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D100_vs_RO_D45/"
pdf(paste0(output_plots_dir,"switch_RO_D100_vs_RO_D45_MAP4.pdf"))
plot <- switchPlot(
  SwitchList_part1_D100_D45,
  gene="MAP4",
  plotTopology=FALSE
)
print(plot)
dev.off()

pdf(paste0(output_plots_dir,"switch_RO_D100_vs_RO_D45_top_500_dtu_events_genes.pdf"))
for(gene_id in top_genes_D100_D45$gene_id){
  plot <- switchPlot(
    SwitchList_part1_D100_D45,
    gene=gene_id,
    plotTopology=FALSE
  )
  
}
print(plot)
dev.off()
