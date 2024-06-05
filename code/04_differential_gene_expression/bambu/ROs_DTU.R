library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)



bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation"

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

SwitchList <- importRdata(isoformCountMatrix   = counts,
                          isoformRepExpression = cpm,
                          designMatrix         = myDesign,
                          isoformExonAnnoation = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/ROs_extended_annotation/extended_annotations.gtf",
                          isoformNtFasta       = "/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa",
                          removeNonConvensionalChr = TRUE,
                          ignoreAfterBar = TRUE,
                          ignoreAfterPeriod = FALSE,
                          showProgress = TRUE)


#remove version number from gene_id
SwitchList$isoformFeatures$gene_id <- gsub("\\..*", "", SwitchList$isoformFeatures$gene_id)

###add gene names 

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)


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

# > unique(SwitchList$isoformFeatures$condition_1)
# [1] "RO_D100" "RO_D200"
# > unique(SwitchList$isoformFeatures$condition_2)
# [1] "RO_D200" "RO_D45" 

D100 vs D45 
D200 vs D45
D100 vs D200


#check if "NF1" is in SwitchList$isoformFeatures$gene_name
"NF1" %in% SwitchList$isoformFeatures$gene_name

SwitchList$isoformFeatures |> dplyr::filter(is.na(IF1) & IF2 > 0) |> 
  dplyr::select(isoform_id, gene_name, dIF, isoform_switch_q_value,
                iso_value_2, IF1, IF2,
                condition_1, condition_2) |> nrow() 
#38515

SwitchList$isoformFeatures |> dplyr::filter(is.na(IF1) & IF2 > 0) |> 
  dplyr::select(gene_name) |> unique() |> nrow()
#5071

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)


SwitchListFiltered$isoformFeatures |> dplyr::filter(gene_name == "NF1") |> 
  nrow()

#split by conditions 

SwitchList_D100_D45 <- SwitchList
SwitchList_D200_D45 <- SwitchList
SwitchList_D100_D200 <- SwitchList


write_tsv(SwitchList_D100_D45$isoformFeatures, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/all_isoform_Features.tsv")
write_tsv(SwitchList_D200_D45$isoformFeatures, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/all_isoform_Features.tsv")
write_tsv(SwitchList_D100_D200$isoformFeatures, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/all_isoform_Features.tsv")


SwitchList_D100_D45$isoformFeatures <- SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RO_D100" & condition_2 == "RO_D45")

SwitchList_D200_D45$isoformFeatures <- SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RO_D200" & condition_2 == "RO_D45")

SwitchList_D100_D200$isoformFeatures <- SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RO_D100" & condition_2 == "RO_D200")


non_DTUs_D100_D45 <- SwitchList_D100_D45$isoformFeatures |> dplyr::filter((iso_value_1 > 0 & gene_value_2 == 0) | (gene_value_1 == 0 & iso_value_2 > 0)) |>
  dplyr::select(isoform_id, gene_id, gene_name, iso_value_1, iso_value_2, IF1, IF2)
write_tsv(non_DTUs_D100_D45, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/non_DTUs.tsv")

non_DTUs_D200_D45 <- SwitchList_D200_D45$isoformFeatures |> dplyr::filter((iso_value_1 > 0 & gene_value_2 == 0) | (gene_value_1 == 0 & iso_value_2 > 0)) |>
  dplyr::select(isoform_id, gene_id, gene_name, iso_value_1, iso_value_2, IF1, IF2)
write_tsv(non_DTUs_D200_D45, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/non_DTUs.tsv")

non_DTUs_D100_D200 <- SwitchList_D100_D200$isoformFeatures |> dplyr::filter((iso_value_1 > 0 & gene_value_2 == 0) | (gene_value_1 == 0 & iso_value_2 > 0)) |>
  dplyr::select(isoform_id, gene_id, gene_name, iso_value_1, iso_value_2, IF1, IF2)
write_tsv(non_DTUs_D100_D200, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/non_DTUs.tsv")


SwitchList_D100_D45_Filtered <- preFilter(
  switchAnalyzeRlist = SwitchList_D100_D45,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

SwitchList_D200_D45_Filtered <- preFilter(
  switchAnalyzeRlist = SwitchList_D200_D45,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

SwitchList_D100_D200_Filtered <- preFilter(
  switchAnalyzeRlist = SwitchList_D100_D200,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

DEXSeq_SwitchList_D100_D45 <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchList_D100_D45_Filtered,  
                                             reduceToSwitchingGenes=TRUE)

DEXSeq_SwitchList_D200_D45 <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchList_D200_D45_Filtered,  
                                             reduceToSwitchingGenes=TRUE)

DEXSeq_SwitchList_D100_D200 <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchList_D100_D200_Filtered,  
                                             reduceToSwitchingGenes=TRUE)


DEX_Seq_output_writing <- function(dexseq_switchlist, output_dir) { 
  dexseq_switchlist$isoformFeatures$neg_log_10_q <- -log10(dexseq_switchlist$isoformFeatures$isoform_switch_q_value)
  write_tsv(dexseq_switchlist$isoformFeatures,
            file = paste0(output_dir, "DEXSeqSwitchList.tsv"))
  dexseq_switchlist$isoformCountMatrix |> write_tsv(file =  paste0(output_dir,  "isoform_counts.tsv"))
  dexseq_switchlist$isoformRepExpression |> write_tsv(file =  paste0(output_dir, "isoform_abundance.tsv"))

  
  }

DEX_Seq_output_writing(DEXSeq_SwitchList_D100_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D45/")
DEX_Seq_output_writing(DEXSeq_SwitchList_D200_D45, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D200_vs_RO_D45/")
DEX_Seq_output_writing(DEXSeq_SwitchList_D100_D200, "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/RO_D100_vs_RO_D200/")

## Volcano Plots ##

volcano_plot <- function(dexseq_switchlist, output_dir){ 
  dexseq_switchlist_top_20 <- dexseq_switchlist$isoformFeatures[order(abs(dexseq_switchlist$isoformFeatures$isoform_switch_q_value)),]
  #filter by dIF
  dexseq_switchlist_top_20 <- dexseq_switchlist_top_20[abs(dexseq_switchlist_top_20$dIF) > 0.5,]
  #top 20
  dexseq_switchlist_top_20 <- dexseq_switchlist_top_20[1:20,]
  
  pdf(paste0(output_dir,"isoform_volcano.pdf"))
  p <- ggplot(data=dexseq_switchlist$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
    geom_point(
      aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
      size=1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw() + 
    geom_text(aes(label=ifelse(gene_name %in% dexseq_switchlist_top_20$gene_name & abs(dIF) > 0.5, 
                               gene_name,'')),
              hjust=0,vjust=0, size = 1)
  
  print(p)
  dev.off()
  
  }

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

