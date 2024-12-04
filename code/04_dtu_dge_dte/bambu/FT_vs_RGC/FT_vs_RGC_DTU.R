library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(rtracklayer)
library(edgeR)
library(tidyr)
library(dplyr)
library(ggVennDiagram)
library('BSgenome.Hsapiens.UCSC.hg38')


method <- "bambu"
comparison <- "FT_vs_RGC"
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison)
counts <- file.path(matrix_dir, "isoform_counts.RDS") 
counts <- readRDS(counts)
cpm <- file.path(matrix_dir, "isoform_cpm.RDS")
cpm <- readRDS(cpm)


#read CPM_transcript.txt


myDesign  <- data.frame(sampleID = colnames(counts) ,
                        condition = c( "FT", "FT", "RGC", "RGC"),
                        stringsAsFactors = FALSE)

bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"

rdata_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
method, comparison, "rds", "SwitchList.rds")
if(!file.exists(rdata_path)){ 
  SwitchList <- importRdata(isoformCountMatrix   = counts,
                            isoformRepExpression = cpm,
                            designMatrix         = myDesign,
                            isoformExonAnnoation = paste0(bambu_dir, "/extended_annotations.gtf"),
                            isoformNtFasta       = paste0(bambu_dir, "/sqanti3_qc/all_samples_corrected.fasta"),
                            removeNonConvensionalChr = TRUE,
                            ignoreAfterBar = TRUE,
                            ignoreAfterPeriod = FALSE,
                            showProgress = TRUE)
  
  
  #if cds gtf doesn't exist, convert gff to gtf
  if(!file.exists( paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf"))) {
    gff <- rtracklayer::import(paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gff"))
    rtracklayer::export(gff, paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf"),
                        "gtf")
    rm(gff)
  }
  
  
  SwitchList <- addORFfromGTF(
    switchAnalyzeRlist     = SwitchList,
    pathToGTF              =  paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf")
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

summary(SwitchList)

#### DTE ####
DTE_table <- readr::read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                             method, comparison, "DTE_table.tsv"))

#### DGE ####
DGE_table <- readr::read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                                         method, comparison, "DGE_table.tsv"))
#### DEXSeq ####
dtu_rdata_path  = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison, "rds/DexSeqDTUDGESwitchList.rds")
if(!file.exists(dtu_rdata_path)){
  
  SwitchList_part1 <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist         = SwitchListFiltered,
    reduceToSwitchingGenes     = FALSE
  )

  ### Add DTE/DGE to switchList
  # SwitchList_part1$isoformFeatures$isoform_id <- gsub("\\..*", "", SwitchList_part1$isoformFeatures$isoform_id)
  idx = match(SwitchList_part1$isoformFeatures$isoform_id, DTE_table$isoform_id)
  SwitchList_part1$isoformFeatures$iso_q_value = DTE_table$FDR[idx]
  
  idx = match(SwitchList_part1$isoformFeatures$gene_id, DGE_table$gene_id)
  SwitchList_part1$isoformFeatures$gene_q_value = DGE_table$FDR[idx]
  
  
  SwitchList_part1 <- analyzeORF(SwitchList_part1, genomeObject = Hsapiens)

  SwitchList_part1$aaSequence = NULL
 
  
  if(!file.exists(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "fastas"))){
    dir.create(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "fastas"))
  }
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "fastas"),
    extractNTseq       = TRUE,
    extractAAseq       = TRUE,
    removeShortAAseq   = TRUE,
    removeLongAAseq    = TRUE, #FOR PFAM 
    onlySwitchingGenes = TRUE,
    alsoSplitFastaFile=TRUE #FOR PFAM
  )
  saveRDS(SwitchList_part1, file = dtu_rdata_path)

}else{
  SwitchList_part1 <- readRDS(dtu_rdata_path)
}

summary(SwitchList_part1)


###### extra #######
  
# #similar to table s3 in the paper
# 
# if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")){
#   # SwitchList_part1$isoformSwitchAnalysis$isoform_id <- gsub("\\..*", "", SwitchList_part1$isoformSwitchAnalysis$isoform_id)
#   DGE_DTU_DTE = SwitchList_part1$isoformFeatures |>
#     as_tibble() |>
#     dplyr::select(isoform_id, gene_id, gene_name, condition_1, condition_2) |>
#     left_join(
#       SwitchList_part1$isoformSwitchAnalysis |> dplyr::select(isoform_id, dIF, pvalue, padj)
#     ) |>
#     dplyr::rename(
#       DTU_dIF    = "dIF",
#       DTU_pval   = "pvalue",
#       DTU_qval   = "padj"
#     ) |>
#     mutate(
#       DTU = DTU_qval < 0.05 # & abs(DTU_dIF) > 0.1
#     ) |>
#     left_join(
#       DTE_table |> dplyr::select(isoform_id, logFC, PValue, FDR)
#     ) |>
#     dplyr::rename(
#       DTE_log2FC = "logFC",
#       DTE_pval   = "PValue",
#       DTE_qval   = "FDR"
#     ) |>
#     mutate(
#       DTE = DTE_qval < 0.05
#     ) |>
#     left_join(
#       DGE_table |> dplyr::select(gene_id, logFC, PValue, FDR)
#     ) |>
#     dplyr::rename(
#       DGE_log2FC = "logFC",
#       DGE_pval   = "PValue",
#       DGE_qval   = "FDR"
#     ) |>
#     mutate(
#       DGE = DGE_qval < 0.05
#     )
#   
#   write_tsv( DGE_DTU_DTE, file = "/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")
# } else {
#   library(here)
#   DGE_DTU_DTE <- readr::read_tsv(here("processed_data/dtu/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv"))
# }



#### VENN DIAGRAM ####
# 
# gene_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> 
#   summarise(DTE = any(DTE), DGE=any(DGE), DTU=any(DTU))  |> dplyr::select(-gene_id)
# 
# dge_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> dplyr::select( -DTE) |> 
#   summarise(DGE=any(DGE), DTU=any(DTU))  
# 
# dte_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> dplyr::select(-DGE) |> 
#   summarise(DTE = any(DTE), DTU=any(DTU)) 
# #total 10673 genes after filtering 
# dge_contingency_table <- as_tibble(xtabs(~ DGE + DTU, data = dge_overlaps))
# dge_contingency_table <- plyr::ddply(dge_contingency_table, ~DGE, transform, Prop = n / sum(n))
# 
# 
# dge_bar <- ggplot(dge_contingency_table, aes(x = DGE, y = n, fill = DTU)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(title = "DGE and DTU genes", x = "DGE", y = "Count") + 
#   scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
#   theme_minimal() +
#   geom_text(aes(label = n), 
#             position = position_stack(vjust = 0.5), 
#             size = 2)
# 
# dte_contingency_table <- xtabs(~ DTE + DTU, data = dte_overlaps)
# dte_contingency_table <- plyr::ddply(as_tibble(dte_contingency_table), ~DTE, transform, Prop = n / sum(n))
# 
# dte_bar <- ggplot(dte_contingency_table, aes(x = DTE, y = n, fill = DTU)) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(title = "DTE and DTU genes", x = "DTE", y = "Count") + 
#   scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
#   theme_minimal() +
#   geom_text(aes(label = n), 
#             position = position_stack(vjust = 0.5), 
#             size = 2)
# library(patchwork)
# p <- dge_bar + dte_bar + plot_annotation(title = 'FT vs RGC')
# ggsave(path = "/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/plots/",
#        device = "pdf", plot = p, filename = "DGE_DTU_DTE_barplot.pdf")
# 
# # save venn diagram as pdf 
# pdf("/users/sparthib/retina_lrs/processed_data/dtu/bambu/FT_vs_RGC/plots/DGE_DTU_DTE_venn.pdf")
# fig <- ggVennDiagram(list(DTU = which(gene_overlaps$DTU), 
#                           DGE = which(gene_overlaps$DGE),
#                           DTE = which(gene_overlaps$DTE))) + 
#   scale_fill_gradient(low="grey",high = "red")
# print(fig)
# dev.off()



#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)


switchlist_part2_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                  method, comparison, "rds", "SwitchList_part2.rds")


# 
# SwitchList_part2 <- readRDS(switchlist_part2_path)

plots_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "plots")
if(!file.exists(plots_dir)){
  dir.create(plots_dir)
}
pdf(file.path(plots_dir, "Splicing_Summary.pdf"))
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = T,
                                           plot = F)
write_tsv(splicing_summary, file = file.path(plots_dir, "Splicing_Summary.tsv"))

pdf(file.path(plots_dir, "Splicing_Enrichment.pdf"))
splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = F ,
  onlySigIsoforms = T,
  countGenes = F
)
print(splicing_enrichment)
dev.off()




external_protein_analyses_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                           method, comparison, "external_protein_analyses")
read_csv(file.path(external_protein_analyses_dir, "pfam_results.csv")) |> 
  write_tsv(file.path(external_protein_analyses_dir, "pfam_results.txt"))

#### Consequence Switch Plots ####

SwitchList_part2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_part2, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = file.path(external_protein_analyses_dir, "CPC2_output.txt"),
  pathToPFAMresultFile      = file.path(external_protein_analyses_dir,"pfam_results.txt"),
  pathToSignalPresultFile   = file.path(external_protein_analyses_dir,"prediction_results.txt"),
  outputPlots               = TRUE,
  pathToOutput              = file.path(external_protein_analyses_dir,"switchplots_with_consequences"),
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  ))

saveRDS(SwitchList_part2, file = switchlist_part2_path)



pdf(file.path(plots_dir, "Splicing_Summary.pdf"))
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()

splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = T,
                                           plot = F)  
write_tsv(splicing_summary, file = file.path(plots_dir, "Splicing_Summary.tsv"))


pdf(file.path(plots_dir, "Splicing_Enrichment.pdf"))
splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = F ,
  onlySigIsoforms = T,
  countGenes = F
)
print(splicing_enrichment)
dev.off()

splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = T,
  onlySigIsoforms = T,
  countGenes = F
)
write_tsv(splicing_enrichment, file = file.path(plots_dir, "Splicing_Enrichment.tsv"))

pdf(file.path(plots_dir, "Consequence_Enrichment.pdf"),
    width = 10, height = 7)
p <- extractConsequenceEnrichment(
  SwitchList_part2,
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  ),
  analysisOppositeConsequence = TRUE,
  localTheme = theme_bw(base_size = 14), # Increase font size in vignette
  returnResult = F, # if TRUE returns a data.frame with the summary statistics
  countGenes = F
)
print(p)
dev.off()

consequences <- extractConsequenceEnrichment(
  SwitchList_part2,
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  ),
  analysisOppositeConsequence = TRUE,
  localTheme = theme_bw(base_size = 14), # Increase font size in vignette
  returnResult = T, # if TRUE returns a data.frame with the summary statistics
  countGenes = F
)
write_tsv(consequences, file = file.path(plots_dir, "Consequence_Enrichment.tsv"))
