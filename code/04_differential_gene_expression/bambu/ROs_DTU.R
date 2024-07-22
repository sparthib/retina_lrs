library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(here)

bambu_dir <- here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation")

#read counts_transcript.txt table from bambu_dir

counts <- read.table(file.path(bambu_dir, "counts_transcript.txt"),
                     header = TRUE)
head(counts)

#remove "_primary_over_30_chr_only_sorted" in column names 
colnames(counts) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(counts))
colnames(counts)[1] <- "isoform_id"
#remove GENE_ID (2nd column)
counts <- counts[, -2]
head(counts)


#read CPM_transcript.txt
cpm <- read.table(file.path(bambu_dir, "CPM_transcript.txt"),
                  header = TRUE)
colnames(cpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(cpm))
colnames(cpm)[1] <- "isoform_id"
cpm <- cpm[, -2]
head(cpm)

myDesign  <- data.frame(sampleID = c("EP1.BRN3B.RO" , "EP1.WT_hRO_2", "EP1.WT_ROs_D45", 
                                     "H9.BRN3B_hRO_2",  "H9.BRN3B.RO", "H9.CRX_hRO_2", "H9.CRX_ROs_D45") ,
                        condition = c("RO_D200", "RO_D100", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45"),
                        stringsAsFactors = FALSE)

dir.create(here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/"), showWarnings = T,
           recursive = T)
rdata_path = here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/SwitchList.rds")

if(!file.exists(rdata_path)){
SwitchList <- importRdata(isoformCountMatrix   = counts,
                          isoformRepExpression = cpm,
                          designMatrix         = myDesign,
                          isoformExonAnnoation = here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/extended_annotations.gtf"),
                          isoformNtFasta       = here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/sqanti3_qc/ROs_corrected.fasta"),
                          removeNonConvensionalChr = TRUE,
                          ignoreAfterBar = TRUE,
                          ignoreAfterPeriod = FALSE,
                          showProgress = TRUE)
#if cds gtf doesn't exist, convert gff to gtf
if(!file.exists(here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/sqanti3_qc/ROs_corrected.gtf.cds.gtf"))) {
  gff <- rtracklayer::import(here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/sqanti3_qc/ROs_corrected.gtf.cds.gff"))
  rtracklayer::export(gff, here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/sqanti3_qc/ROs_corrected.gtf.cds.gtf"),
                      "gtf")
  rm(gff)
}


SwitchList <- addORFfromGTF(
  switchAnalyzeRlist     = SwitchList,
  pathToGTF              = here("processed_data/DTU_Gandall/bambu/ROs_extended_annotation/sqanti3_qc/ROs_corrected.gtf.cds.gtf")
)

SwitchList$isoformFeatures$gene_id <- gsub("\\..*", "", SwitchList$isoformFeatures$gene_id)
require("biomaRt")
us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
mart <- useDataset("hsapiens_gene_ensembl", us_mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_gene_id",
                "hgnc_symbol"),
  filter="ensembl_gene_id",
  values=SwitchList$isoformFeatures$gene_id,
  uniqueRows=TRUE)

colnames(annotLookup) <- c("gene_id", "ensembl_gene_name")
SwitchList$isoformFeatures <- merge(SwitchList$isoformFeatures, annotLookup, by="gene_id", all.x=TRUE)


head(SwitchList$isoformFeatures)
SwitchList$isoformFeatures$gene_name <- SwitchList$isoformFeatures$ensembl_gene_name

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist         = SwitchList,
  geneExpressionCutoff       = 1,     # default
  isoformExpressionCutoff    = 0,     # default
  removeSingleIsoformGenes   = TRUE  # default
)

saveRDS(SwitchListFiltered, file = rdata_path)
#save DEXSeq switchlist

}else { 
  SwitchListFiltered <- readRDS(rdata_path)
}

# if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/DTE_table.tsv")){
#   
dtu_rdata_path  = here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/DexSeqDTUDGESwitchList.rds")
dir.create(here("processed_data/dtu/DTU_gandall/bambu/ROs/fastas/"), showWarnings = T, recursive = T)
if(!file.exists(dtu_rdata_path)){
  SwitchList_part1 <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist         = SwitchListFiltered,
    reduceToSwitchingGenes     = FALSE
  )
  
  idx = match(SwitchList_part1$isoformFeatures$isoform_id, DTE_table$isoform_id)
  SwitchList_part1$isoformFeatures$iso_q_value = DTE_table$FDR[idx]
  
  idx = match(SwitchList_part1$isoformFeatures$gene_id, DGE_table$gene_id)
  SwitchList_part1$isoformFeatures$gene_q_value = DGE_table$FDR[idx]
  
  SwitchList_part1$isoformFeatures <- SwitchList_part1$isoformFeatures |> distinct() 
  
  # Switching features:
  #   Comparison Isoforms Switches Genes
  # 1  RO_D100 vs RO_D45     2076     1588  1495
  # 2 RO_D100 vs RO_D200     1684     1389  1290
  # 3  RO_D200 vs RO_D45     2490     2055  1775
  # 4           Combined     4804     4417  3188
  
  SwitchList_part1 <- analyzeORF(SwitchList_part1, genomeObject = Hsapiens)
  
  SwitchList_part1$aaSequence = NULL
  
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = here("processed_data/dtu/DTU_gandall/bambu/ROs/fastas/"),
    extractNTseq       = TRUE, #for CPC2
    extractAAseq       = TRUE,
    removeShortAAseq   = TRUE,
    removeLongAAseq    = TRUE, #FOR PFAM, SignalP
    onlySwitchingGenes = TRUE,
    alsoSplitFastaFile=FALSE
  )
  
  saveRDS(SwitchList_part1, file = dtu_rdata_path)
}else{
  SwitchList_part1 <- readRDS(dtu_rdata_path)
}


if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE.tsv")){
  
  
  ### load DTEs ###
  D100_vs_D200_DTE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DTE/D100_vs_D200_DTEs.tsv")
  D200_vs_D45_DTE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DTE/D200_vs_D45_DTEs.tsv")
  D100_vs_D45_DTE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DTE/D100_vs_D45_DTEs.tsv")
  
  #rbind 
  DTE_table <- rbind(D100_vs_D200_DTE_table, D200_vs_D45_DTE_table, D100_vs_D45_DTE_table)
  
  ### load DGEs ###
  D100_vs_D200_DGE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE/D100_vs_D200_DGEs.tsv")
  D200_vs_D45_DGE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE/D200_vs_D45_DGEs.tsv")
  D100_vs_D45_DGE_table <- read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE/D100_vs_D45_DGEs.tsv")
  
  #rbind
  DGE_table <- rbind(D100_vs_D200_DGE_table, D200_vs_D45_DGE_table, D100_vs_D45_DGE_table)
  
  #check for duplicates 
  SwitchList_part1$isoformFeatures <- SwitchList_part1$isoformFeatures |> distinct() 
  
  DGE_DTU_DTE = SwitchList_part1$isoformFeatures |>
    as_tibble() |>
    dplyr::select(isoform_id, gene_id, gene_name, condition_1, condition_2) |>
    left_join(
      SwitchList_part1$isoformSwitchAnalysis |> dplyr::select(isoform_id, dIF, pvalue, padj, condition_1, condition_2),
      by = c("isoform_id", "condition_1", "condition_2")
    ) |>
    dplyr::rename(
      DTU_dIF    = "dIF",
      DTU_pval   = "pvalue",
      DTU_qval   = "padj"
    ) |>
    mutate(
      DTU = DTU_qval < 0.05 # & abs(DTU_dIF) > 0.1
    ) 
  
  DGE_DTU_DTE$isoform_id <- gsub("\\..*", "", DGE_DTU_DTE$isoform_id)
  
  DGE_DTU_DTE  <- DGE_DTU_DTE |>
    left_join(
      DTE_table , by = c("isoform_id", "condition_1", "condition_2")
    ) |>
    dplyr::rename(
      DTE_log2FC = "logFC",
      DTE_pval   = "PValue",
      DTE_qval   = "FDR"
    ) |>
    mutate(
      DTE = DTE_qval < 0.05
    ) 
  DGE_DTU_DTE <- DGE_DTU_DTE |> dplyr::select(-c(logCPM, F))
  
  DGE_DTU_DTE$gene_id <- DGE_DTU_DTE$gene_id.x
  DGE_DTU_DTE$gene_name <- DGE_DTU_DTE$gene_name.x
  DGE_DTU_DTE <- DGE_DTU_DTE |> dplyr::select(-c(gene_id.x,gene_name.x, gene_id.y, gene_name.y))
  colnames(DGE_DTU_DTE) 
  
  DGE_DTU_DTE  <- DGE_DTU_DTE |>
    left_join(
      DGE_table, by = c("gene_id", "condition_1", "condition_2")
    ) |>
    dplyr::rename(
      DGE_log2FC = "logFC",
      DGE_pval   = "PValue",
      DGE_qval   = "FDR"
    ) |>
    mutate(
      DGE = DGE_qval < 0.05
    )
  colnames(DGE_DTU_DTE)
  

  DGE_DTU_DTE$gene_name <- DGE_DTU_DTE$gene_name.x
  DGE_DTU_DTE <- DGE_DTU_DTE |> dplyr::select(-c(gene_name.x,  gene_name.y, gene_name, gene_biotype, transcript_biotype))
  
  annotLookup <- getBM(
    mart=mart,
    attributes=c( "ensembl_gene_id",
                  "external_gene_name",
                  "gene_biotype"),
    filter="ensembl_gene_id",
    values=DGE_DTU_DTE$gene_id,
    uniqueRows=TRUE)
  
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype")
  
  DGE_DTU_DTE <- merge(DGE_DTU_DTE, annotLookup,
                    by="gene_id", all.x=TRUE)
  
  # Set the Ensembl mirror to "useast"
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  
  # Ensure DGE_DTU_DTE$isoform_id is a vector
  isoform_ids <- as.vector(DGE_DTU_DTE$isoform_id)
  
  # Query the database
  # annotLookup <- getBM(
  #   mart = mart,
  #   attributes = c("ensembl_transcript_id", "transcript_biotype"),
  #   filters = "ensembl_transcript_id",
  #   values = isoform_ids,
  #   uniqueRows = TRUE
  # )
  # 
  # colnames(annotLookup) <- c("isoform_id", "transcript_biotype")
  # tt$table <- merge(tt$table, annotLookup,
  #                   by="isoform_id", all.x=TRUE)
  # 
  
  
  
  
  write_tsv( DGE_DTU_DTE, file = "./processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE.tsv")
} else {
  DGE_DTU_DTE <- readr::read_tsv("./processed_data/dtu/DTU_gandall/bambu/ROs/DGE_DTU_DTE.tsv")
}



#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)

saveRDS(SwitchList_part2, file = here("processed_data/dtu/DTU_gandall/bambu/ROs/rds/SwitchList_part2.rds"))

# dir.create("./processed_data/dtu/DTU_gandall/bambu/ROs/plots", showWarnings = T, recursive = T)
pdf("./processed_data/dtu/DTU_gandall/bambu/ROs/plots/Splicing_Summary.pdf")
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()
pdf("./processed_data/dtu/DTU_gandall/bambu/ROs/plots/Splicing_Enrichment.pdf",
    width = 10, height = 7)
splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = F ,
  onlySigIsoforms = T,
  countGenes = F
)
print(splicing_enrichment)
dev.off()



################################################################################



volcano_plot(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D100_vs_RO_D45/")
volcano_plot(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D200_vs_RO_D45/")
volcano_plot(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D100_vs_RO_D200/")


volcano_df <- function(dexseq_switchlist, output_dir){
  
  Volcano_df <- dexseq_switchlist$isoformFeatures |> 
    dplyr::filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)
  
  write_tsv(Volcano_df,
            file =  paste0(output_dir, "switches.tsv"))
  
}

volcano_df(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/")
volcano_df(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/")
volcano_df(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/")



####.SWITCH VS GENE CHANGES ####

switch_vs_gene_plot <- function(dexseq_switchlist, output_dir){ 
  dexseq_switchlist_top_20 <- dexseq_switchlist$isoformFeatures[order(abs(dexseq_switchlist$isoformFeatures$isoform_switch_q_value)),]
  #filter by dIF
  dexseq_switchlist_top_20 <- dexseq_switchlist_top_20[abs(dexseq_switchlist_top_20$dIF) > 0.5,]
  #top 20
  dexseq_switchlist_top_20 <- dexseq_switchlist_top_20[1:20,]
  
  pdf(paste0(output_dir,"switch_vs_degs.pdf"))
  s <- ggplot(data=dexseq_switchlist$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
      size=1
    )  + 
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw() + 
    geom_text(aes(label=ifelse(gene_name %in% dexseq_switchlist_top_20$gene_name & abs(dIF) > 0.5, 
                               gene_name,'')),
              hjust=0,vjust=0, size = 1)
    print(s)
    dev.off()
    
}

switch_vs_gene_plot(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D100_vs_RO_D45/")
switch_vs_gene_plot(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D200_vs_RO_D45/")
switch_vs_gene_plot(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/plots/de/switch_analyzer/bambu/RO_D100_vs_RO_D200/")


switch_vs_deg_data <- function(dexseq_switchlist, output_dir){
  
  switch_vs_degs <- dexseq_switchlist$isoformFeatures |> 
    dplyr::filter(abs(dIF) > 0.1 & abs(gene_log2_fold_change) < 2 & isoform_switch_q_value < 0.05)
  
  write_tsv(switch_vs_degs,
            file =  paste0(output_dir, "switch_vs_degs.tsv"))
  
}

switch_vs_deg_data(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/")
switch_vs_deg_data(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/")
switch_vs_deg_data(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/")


microexon_data <- function(dexseq_switchlist, output_dir){ 
  exons <- as.data.frame(dexseq_switchlist$exons)
  microexons <- exons |> dplyr::filter(width < 27)
  
  #remove version number from gene_id
  microexons$gene_id <- gsub("\\..*", "", microexons$gene_id)
  microexons <- merge(microexons, annotLookup, by="gene_id", all.x=TRUE)
  microexons$gene_name <- microexons$ensembl_gene_name
  #drop ensemble_gene_name
  microexons <- microexons |> dplyr::select(-ensembl_gene_name)
  
  #all microexons found in all genes that are expressed in the samples. 
  write_tsv(microexons, paste0(output_dir, "/microexons/all_micro_exons.tsv"))
  
  DEXSeq <- dexseq_switchlist$isoformFeatures |>  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |> 
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

microexon_data(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/")
microexon_data(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/")
microexon_data(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/")


###### MAKE SWITCH PLOTS #######

## D100 vs D45
top_genes_D100_D45 <- extractTopSwitches(DEXSeq_SwitchList_D100_D45, n=500, 
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
  DEXSeq_SwitchList_D100_D45,
  gene="MAP4",
  plotTopology=FALSE
)
print(plot)
dev.off()

pdf(paste0(output_plots_dir,"switch_RO_D100_vs_RO_D45_top_500_dtu_events_genes.pdf"))
for(gene_id in top_genes_D100_D45$gene_id){
  plot <- switchPlot(
    DEXSeq_SwitchList_D100_D45,
    gene=gene_id,
    plotTopology=FALSE
  )
  
}
print(plot)
dev.off()
