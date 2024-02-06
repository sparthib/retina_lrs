library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)

# 
# salmonDf <- prepareSalmonFileDataFrame("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/mapping_mode")
# 
# salmonDf$condition <- c("RGC", "RO_D209", "RO_D45", "RO_D209" ,"RO_D45",
#                         "RGC", "RGC", "RGC" )
# 
# # salmonDf
# 
# # names condition
# # 1     DG-WT-hRGC       RGC
# # 2   EP1-BRN3B-RO   RO_D209
# # 3 EP1-WT_ROs_D45    RO_D45
# # 4    H9-BRN3B-RO   RO_D209
# # 5 H9-CRX_ROs_D45    RO_D45
# # 6      H9-hRGC_1       RGC
# # 7      H9-hRGC_2       RGC
# # 8           hRGC       RGC


salmonQuant <- importIsoformExpression(parentDir = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/mapping_mode",
                                       addIsofomIdAsColumn = TRUE)

# Step 3 of 3: Normalizing abundance values (not counts) via edgeR...
# Done

head(salmonQuant$abundance, 2)
head(salmonQuant$counts, 2)

myDesign  <- data.frame(sampleID = c("DG-WT-hRGC", "EP1-BRN3B-RO", "EP1-WT_ROs_D45", "H9-BRN3B-RO"  ,
                                     "H9-CRX_ROs_D45","H9-hRGC_1" , "H9-hRGC_2", "hRGC"),
                        condition = c("RGC", "RO_D209", "RO_D45", "RO_D209" ,"RO_D45",  "RGC", "RGC", "RGC" ),
                        stringsAsFactors = FALSE)

#### IMPORT DATA ####


SwitchList <- importRdata(isoformCountMatrix   = salmonQuant$counts,
                               isoformRepExpression = salmonQuant$abundance,
                               designMatrix         = myDesign,
                               isoformExonAnnoation = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz",
                               isoformNtFasta       = "/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa",
                               removeNonConvensionalChr = TRUE,
                               ignoreAfterBar = TRUE,
                               ignoreAfterPeriod = TRUE,
                               showProgress = TRUE)


# 
# Warning messages:
#   1: In importRdata(isoformCountMatrix = salmonQuant$counts, isoformRepExpression = salmonQuant$abundance,  :
#                       The annotation and quantification (count/abundance matrix and isoform annotation) Seem to be slightly different.
#                     Specifically:
#                       305 isoforms were only found in the annotation
#                     
#                     Please make sure this is on purpouse since differences will cause inaccurate quantification and thereby skew all analysis.
#                     If you have quantified with Salmon this could be normal since it as default only keep one copy of identical sequnces (can be prevented using the --keepDuplicates option)
#                     We strongly encurage you to go back and figure out why this is the case.
#                     
                    
                    
#### FILTER BASED ON SINGLE ISOFORM GENES AND GENE COUNTS ####
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = SwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

# The filtering removed 44601 ( 46.5% of ) transcripts. There is now 51310 isoforms left

#### DEXSeq SwitchList ####

DEXSeq_SwitchList <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchListFiltered,  
                                             reduceToSwitchingGenes=TRUE)

#FDR cutoff = 0.05 by default


DEXSeq_SwitchList$isoformFeatures$neg_log_10_q <- -log10(DEXSeq_SwitchList$isoformFeatures$isoform_switch_q_value)
write_tsv(DEXSeq_SwitchList$isoformFeatures,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/DEXSeqSwitchList_Feb4.tsv")


write_tsv(DEXSeq_SwitchList$isoformCountMatrix,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/isoform_counts.tsv")

write_tsv(DEXSeq_SwitchList$isoformRepExpression,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/isoform_abundance.tsv")


top_genes <- extractTopSwitches(DEXSeq_SwitchList, n=500, 
                                alpha = 0.05,
                                dIFcutoff = 0.1)
write_tsv(top_genes,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/top_genes.tsv")



###### MAKE SWITCH PLOTS #######
## RGC vs RO_D209

top_genes_RGC_vs_ROD209 <- top_genes |> 
  dplyr::filter(condition_1 =="RGC" & condition_2 == "RO_D209") 



pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/mapping_mode/switch_RGC_vs_ROD209.pdf")

for(gene_name in top_genes_RGC_vs_ROD209$gene_name){
  plot <- switchPlot(
    DEXSeq_SwitchList,
    gene=gene_name,
    condition1 = 'RGC',
    condition2 = 'RO_D209',
    plotTopology=FALSE
  ) 
  
}
print(plot)
dev.off()


## RGC vs RO_D45 ##
top_genes_RGC_vs_RO_D45 <- top_genes |> 
  dplyr::filter(condition_1 =="RGC" & condition_2 == "RO_D45")


pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/mapping_mode/switch_RGC_vs_ROD45.pdf")
for(gene_name in top_genes_RGC_vs_RO_D45$gene_name){
  plot <- switchPlot(
    DEXSeq_SwitchList,
    gene=gene_name,
    condition1 = 'RGC',
    condition2 = 'RO_D45',
    plotTopology=FALSE
  ) 
  
}
print(plot)
dev.off()

## RO_D209 vs RO_D45

top_genes_RO_D209_vs_RO_D45 <- top_genes |> 
  dplyr::filter(condition_1 =="RO_D209" & condition_2 == "RO_D45")


pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/mapping_mode/switch_RO_D209_vs_D45.pdf")
for(gene_name in top_genes_RO_D209_vs_RO_D45$gene_name){
  plot <- switchPlot(
    DEXSeq_SwitchList,
    gene=gene_name,
    condition1 = 'RO_D209',
    condition2 = 'RO_D45',
    plotTopology=FALSE
  ) 
  
}
print(plot)
dev.off()


## Volcano Plots ##

pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/mapping_mode/isoform_volcano.pdf")

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
  geom_text(aes(label=ifelse(abs(dIF)> 0.5 & abs(isoform_switch_q_value) < 0.05 ,
                             gene_name,'')),
            hjust=0,vjust=0, size = 1)

print(p)
dev.off()

nrow(DEXSeq_SwitchList$isoformFeatures |> filter(abs(dIF) > 0.1 & isoform_switch_q_value < 0.05))

Volcano_df_RGC_vs_ROD209 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RGC" & condition_2 == "RO_D209") |>
  dplyr::filter(isoform_switch_q_value < 0.05 & dIF > 0.1)

write_tsv(Volcano_df_RGC_vs_ROD209,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/RGC_vs_ROD209_switches.tsv")



Volcano_df_RGC_vs_ROD45 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RGC" & condition_2 == "RO_D45") |>
  dplyr::filter(isoform_switch_q_value < 0.05 & dIF > 0.1)
write_tsv(Volcano_df_RGC_vs_ROD45,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/RGC_vs_ROD45_switches.tsv")



Volcano_df_ROD209_vs_ROD45 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RO_D209" & condition_2 == "RO_D45") |>
  dplyr::filter(isoform_switch_q_value < 0.05 & dIF > 0.1)

write_tsv(Volcano_df_ROD209_vs_ROD45,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/ROD209_vs_ROD45_switches.tsv")




####.SWITCH VS GENE CHANGES ####
pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/mapping_mode/switch_vs_degs.pdf")
s <- ggplot(data=DEXSeq_SwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('grey','skyblue')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw() + 
  geom_text(aes(label=ifelse(abs(dIF)> 0.5 & abs(gene_log2_fold_change) < 2 ,
                             gene_name,'')),
            hjust=0,vjust=0, size = 1)
print(s)
dev.off()



switch_vs_degs_RGC_vs_ROD209 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RGC" & condition_2 == "RO_D209") |>
  dplyr::filter(abs(dIF)> 0.5 & abs(gene_log2_fold_change) < 2)

write_tsv(switch_vs_degs_RGC_vs_ROD209 ,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/switch_vs_degs_RGC_vs_ROD209.tsv")



switch_vs_degs_RGC_vs_ROD45 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RGC" & condition_2 == "RO_D45") |>
  dplyr::filter(abs(dIF)> 0.5 & abs(gene_log2_fold_change) < 2)
write_tsv(switch_vs_degs_RGC_vs_ROD45,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/switch_vs_degs_RGC_vs_ROD209.tsv")



switch_vs_degs_ROD209_vs_ROD45 <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(condition_1 == "RO_D209" & condition_2 == "RO_D45") |>
  dplyr::filter(abs(dIF)> 0.5 & abs(gene_log2_fold_change) < 2)
write_tsv(switch_vs_degs_ROD209_vs_ROD45,
          file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/salmon_mapping_mode/switch_vs_degs_RGC_vs_ROD209.tsv")



########## MICROEXONS ###############


exons <- as.data.frame(DEXSeq_SwitchList$exons)

micro_exons <- exons |> dplyr::filter(width <= 27)

nrow(micro_exons)
# 1122

Switched_genes <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(isoform_switch_q_value < 0.05 & color=abs(dIF) > 0.1 )




##### PVALUES ######

zero_qvals <- DEXSeq_SwitchList$isoformFeatures |> 
  dplyr::filter(isoform_switch_q_value == 0 )

pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/qval_hist.pdf")
p <- hist(DEXSeq_SwitchList$isoformFeatures$isoform_switch_q_value,
          col = "skyblue", main = "Histogram of Q Values",
          breaks = 100, xlab = "Q Vals")
print(p)
dev.off()





