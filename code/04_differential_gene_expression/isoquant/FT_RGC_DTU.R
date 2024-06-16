library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)



input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT"
output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/isoquant/FT_vs_RGC"

#read counts 
counts <- read.table(file.path(input_dir, "OUT.transcript_model_grouped_counts.tsv"),
                     header = FALSE, sep = "\t")
colnames(counts) <- c( "isoform_id" , "H9.FT_1_primary_over_30_chr_only_sorted", "H9.FT_2_primary_over_30_chr_only_sorted", "H9.hRGC_1_primary_over_30_chr_only_sorted","H9.hRGC_2_primary_over_30_chr_only_sorted")
#remove "_primary_over_30_chr_only_sorted" in column names
colnames(counts) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(counts))
colnames(counts)[1] <- "isoform_id"


#read transcript_model_grouped_tpm.tsv
tpm <- read.table(file.path(input_dir, "OUT.transcript_model_grouped_tpm.tsv"),
                     header = FALSE, sep = "\t")
head(tpm)
colnames(tpm) <- c( "isoform_id" , "H9.FT_1_primary_over_30_chr_only_sorted", "H9.FT_2_primary_over_30_chr_only_sorted", "H9.hRGC_1_primary_over_30_chr_only_sorted","H9.hRGC_2_primary_over_30_chr_only_sorted")
#remove "_primary_over_30_chr_only_sorted" in column names 
colnames(tpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(tpm))
colnames(tpm)[1] <- "isoform_id"

head(tpm)


myDesign  <- data.frame(sampleID = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2") ,
                        condition = c( "FT", "FT", "RGC", "RGC"),
                        stringsAsFactors = FALSE)


SwitchList <- importRdata(isoformCountMatrix   = counts,
                          isoformRepExpression = tpm,
                          designMatrix         = myDesign,
                          isoformExonAnnoation = paste0(input_dir, "/OUT.transcript_models.gtf"),
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


#only where one of the conditions has iso_value > 0 and the other = 0

non_DTUs <- SwitchList$isoformFeatures |> dplyr::filter((iso_value_1 > 0 & gene_value_2 == 0) | (gene_value_1 == 0 & iso_value_2 > 0)) |>
  dplyr::select(isoform_id, gene_id, gene_name, iso_value_1, iso_value_2, IF1, IF2)
output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/isoquant/FT_vs_RGC/"
write_tsv(non_DTUs, file = paste0(output_data_dir,"non_DTUs.tsv"))

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)



DEXSeq_SwitchList <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchListFiltered,  
                                             reduceToSwitchingGenes=TRUE)

#FDR cutoff = 0.05 by default

output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/isoquant/FT_vs_RGC"

DEXSeq_SwitchList$isoformFeatures$neg_log_10_q <- -log10(DEXSeq_SwitchList$isoformFeatures$isoform_switch_q_value)
write_tsv(DEXSeq_SwitchList$isoformFeatures,
          file = paste0(output_data_dir, "DEXSeqSwitchList.tsv"))


write_tsv(DEXSeq_SwitchList$isoformCountMatrix,
          file =  paste0(output_data_dir,  "isoform_counts.tsv"))

write_tsv(DEXSeq_SwitchList$isoformRepExpression,
          file =  paste0(output_data_dir, "isoform_abundance.tsv"))


top_genes <- extractTopSwitches(DEXSeq_SwitchList, n=500, 
                                alpha = 0.05,
                                dIFcutoff = 0.1)
write_tsv(top_genes,
          file = paste0(output_data_dir,  "top_genes.tsv"))

#count number of nas in gene_name
sum(is.na(top_genes$gene_name))

###### MAKE SWITCH PLOTS #######
## FT vs RGC 

output_plots_dir <- "/users/sparthib/retina_lrs/plots/de/switch_analyzer/isoquant/FT_vs_RGC/"
# pdf(paste0(output_plots_dir,"switch_FT_vs_RGC_top_dtu_events_genes.pdf"))
# for(gene_id in top_genes$gene_id){
#   plot <- switchPlot(
#     DEXSeq_SwitchList,
#     gene=gene_id,
#     condition1 = 'FT',
#     condition2 = 'RGC',
#     plotTopology=FALSE
#   )
#   
# }
# print(plot)
# dev.off()

## Volcano Plots ##
DEXSeq_SwitchList_top_20 <- DEXSeq_SwitchList$isoformFeatures[order(abs(DEXSeq_SwitchList$isoformFeatures$isoform_switch_q_value)),]
#filter by dIF 
DEXSeq_SwitchList_top_20 <- DEXSeq_SwitchList_top_20[abs(DEXSeq_SwitchList_top_20$dIF) > 0.5,]
#top 20
DEXSeq_SwitchList_top_20 <- DEXSeq_SwitchList_top_20[1:20,]

pdf(paste0(output_plots_dir,"isoform_volcano.pdf"))
p <- ggplot(data=DEXSeq_SwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw() + 
  geom_text(aes(label=ifelse(gene_name %in% DEXSeq_SwitchList_top_20$gene_name, 
                             gene_name,'')),
            hjust=0,vjust=0, size = 1)

print(p)
dev.off()


nrow(DEXSeq_SwitchList$isoformFeatures |> filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05))

Volcano_df <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "FT" & condition_2 == "RGC") |>
  dplyr::filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)

write_tsv(Volcano_df,
          file =  paste0(output_data_dir, "switches.tsv"))



####.SWITCH VS GENE CHANGES ####

pdf(paste0(output_plots_dir,"switch_vs_degs.pdf"))
s <- ggplot(data=DEXSeq_SwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  )  + 
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw() + 
  geom_text(aes(label=ifelse(gene_name %in% DEXSeq_SwitchList_top_20$gene_name, 
                             gene_name,'')),
            hjust=0,vjust=0, size = 1)
print(s)
dev.off()


switch_vs_degs_FT_vs_RGC <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "FT" & condition_2 == "RGC") |>
  dplyr::filter(abs(dIF) > 0.1 & abs(gene_log2_fold_change) < 1 & isoform_switch_q_value < 0.05)

write_tsv(switch_vs_degs_FT_vs_RGC ,
          file = paste0(output_data_dir,"switch_vs_degs.tsv"))



########## MICROEXONS ###############


exons <- as.data.frame(DEXSeq_SwitchList$exons)

microexons <- exons |> dplyr::filter(width < 27)

#remove version number from gene_id
microexons$gene_id <- gsub("\\..*", "", microexons$gene_id)
microexons <- merge(microexons, annotLookup, by="gene_id", all.x=TRUE)
microexons$gene_name <- microexons$ensembl_gene_name

#merge microexons with DEXSeq gene names

#all microexons found in all genes that are expressed in the samples. 
write_tsv(microexons, paste0(output_data_dir, "/microexons/all_micro_exons.tsv"))

DEXSeq <- DEXSeq_SwitchList$isoformFeatures |>  filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1) |> 
  dplyr::select(isoform_id, dIF, isoform_switch_q_value,gene_id, gene_name, condition_1, condition_2)

## all microexons found in genes that have some isoform that showed differential usage across any pairwise comparison
microexons_DTU_genes <- microexons |> filter(gene_id %in% DEXSeq$gene_id)
nrow(microexons_DTU_genes)

write_tsv(microexons_DTU_genes, paste0(output_data_dir, "/microexons/microexons_DTU_genes.tsv"))

## all microexons found in isoforms  that showed differential usage across any pairwise comparison 
microexons_DTU_isoforms <- microexons |> filter(isoform_id %in% DEXSeq$isoform_id)
nrow(microexons_DTU_isoforms)


write_tsv(microexons_DTU_isoforms, paste0(output_data_dir, "/microexons/microexons_DTU_isoforms.tsv"))



