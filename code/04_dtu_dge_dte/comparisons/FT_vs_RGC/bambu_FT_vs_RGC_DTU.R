library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(rtracklayer)
library(edgeR)
library(tidyr)
library(dplyr)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")

method <- "bambu"
comparison <- "FT_vs_RGC"
matrix_dir <- file.path(data_dir, "06_quantification/counts_matrices/bambu/FT_vs_RGC/filtered_by_counts_and_biotype")
counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)
cpm <- file.path(matrix_dir, "filtered_isoform_cpm.RDS")
cpm <- readRDS(cpm)

nrow(counts)
nrow(cpm)
#read CPM_transcript.txt
myDesign  <- data.frame(sampleID = colnames(counts) ,
                        condition = c( "FT", "FT", "RGC", "RGC"),
                        stringsAsFactors = FALSE)

bambu_dir <- file.path(data_dir, 
                       "06_quantification/bambu/all_samples_extended_annotation_track_reads")

rdata_path = file.path(code_dir, "processed_data/dtu/",
method, comparison, "protein_coding", "rds", "SwitchList.rds")
dir.create(file.path(code_dir, "processed_data/dtu/",
                     method, comparison, "protein_coding", "rds"), showWarnings = FALSE, recursive =  TRUE)
if(!file.exists(rdata_path)){ 
  SwitchList <- importRdata(isoformCountMatrix   = counts,
                            isoformRepExpression = cpm,
                            designMatrix         = myDesign,
                            isoformExonAnnoation = paste0(bambu_dir, "/FT_vs_RGC_protein_coding_annotations.gtf"),
                            isoformNtFasta       = paste0(bambu_dir, "/sqanti3_qc/all_samples_corrected.fasta"),
                            removeNonConvensionalChr = TRUE,
                            ignoreAfterBar = TRUE,
                            ignoreAfterPeriod = TRUE,
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
    geneExpressionCutoff       = NULL,     
    isoformExpressionCutoff    = NULL,     
    removeSingleIsoformGenes   = FALSE,
    IFcutoff=0
  )
  
  saveRDS(SwitchListFiltered, file = rdata_path)
  SwitchListFiltered$isoformFeatures <- SwitchListFiltered$isoformFeatures |> 
    distinct(across(-gene_name), .keep_all = TRUE)
  
  write_tsv(SwitchListFiltered$isoformFeatures,
            file = file.path(code_dir, "processed_data/dtu",
                             method, comparison, "protein_coding","isoformFeatures.tsv"))
  
  #save DEXSeq switchlist
  
}else { 
  SwitchListFiltered <- readRDS(rdata_path)
    }

summary(SwitchListFiltered)

#### DTE ####
DTE_table <- readr::read_tsv(file.path(code_dir, "processed_data/dtu/", 
                             method, comparison, "protein_coding","DTE_table.tsv"))

#### DGE ####
DGE_table <- readr::read_tsv(file.path(code_dir, "processed_data/dtu/", 
                                         method, comparison, "protein_coding","DGE_table.tsv"))

#### DEXSeq ####
dtu_rdata_path  = file.path(code_dir, "processed_data/dtu/",
                            method, comparison,"protein_coding", "rds/DexSeqDTUDGESwitchList.rds")
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
 
  
  if(!file.exists(file.path(code_dir, "processed_data/dtu/", method, comparison, "protein_coding","fastas"))){
    dir.create(file.path(code_dir,"processed_data/dtu/", method, comparison, "protein_coding","fastas"))
  }
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = file.path(code_dir,"processed_data/dtu/", method, 
                                   comparison, "protein_coding","fastas"),
    extractNTseq       = TRUE,
    extractAAseq       = TRUE,
    removeShortAAseq   = TRUE,
    removeLongAAseq    = TRUE, #FOR PFAM 
    onlySwitchingGenes = TRUE,
    alsoSplitFastaFile=FALSE#FOR PFAM
  )
  saveRDS(SwitchList_part1, file = dtu_rdata_path)

}else{
  SwitchList_part1 <- readRDS(dtu_rdata_path)
}

summary(SwitchList_part1)

# Switching features:
#   Comparison Isoforms Switches Genes
# 1  FT vs RGC     1624     1314  1081

#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)

switchlist_part2_path = file.path(code_dir,"processed_data/dtu/",
                                  method, comparison,"protein_coding", "rds", "SwitchList_part2.rds")


# 
# SwitchList_part2 <- readRDS(switchlist_part2_path)

plots_dir <- file.path(code_dir,"processed_data/dtu/",
                       method, comparison,"protein_coding", "plots")
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

external_protein_analyses_dir <- file.path(code_dir,"processed_data/dtu/",
                                           method, comparison, "external_protein_analyses")
external_protein_ptc_analyses_dir <- file.path(code_dir,"processed_data/dtu/",
                                           method, comparison, "protein_coding","external_protein_analyses")
dir.create(external_protein_ptc_analyses_dir, showWarnings = FALSE, recursive =  TRUE)


external_protein_analyses_dir <- file.path(code_dir,"processed_data/dtu/",
                                           method, comparison, "external_protein_analyses")
protein_analysis <- read_csv(file.path(external_protein_analyses_dir, "pfam_results.csv")) 
protein_analysis$`seq_id` <- gsub("\\..*", "", protein_analysis$`seq_id`)
write_tsv(protein_analysis, 
          file = file.path(external_protein_analyses_dir, "pfam_results.txt"), 
          col_names = TRUE)
#### CPC2 output ####
cpc2_output <- read_tsv(file.path(external_protein_analyses_dir, "CPC2_output.txt"))
cpc2_output$`#ID` <- gsub("\\..*", "", cpc2_output$`#ID`)
write_tsv(cpc2_output, 
          file = file.path(external_protein_analyses_dir, "CPC2_output.txt"), 
          col_names = TRUE)

#### SignalP output ####
signalp_output <- read_tsv(file.path(external_protein_analyses_dir, "prediction_results.txt"),
                           skip = 1)
signalp_output$`# ID` <- gsub("\\..*", "", signalp_output$`# ID`)
write_tsv(signalp_output, 
          file = file.path(external_protein_analyses_dir, "prediction_results.txt"), 
          col_names = TRUE)





#### Consequence Switch Plots ####
SwitchList_part2 <- readRDS(switchlist_part2_path)

SwitchList_part2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_part2, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = file.path(external_protein_analyses_dir, "CPC2_output.txt"),
  pathToPFAMresultFile      = file.path(external_protein_analyses_dir,"pfam_results.txt"),
  pathToSignalPresultFile   = file.path(external_protein_analyses_dir,"prediction_results.txt"),
  outputPlots               = TRUE,
  pathToOutput              = file.path(external_protein_ptc_analyses_dir,"switchplots_with_consequences"),
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

SwitchList_part2$isoformFeatures <- SwitchList_part2$isoformFeatures |>
  distinct(across(-gene_name), .keep_all = TRUE)

write_tsv(SwitchList_part2$isoformFeatures, file = file.path(code_dir,"processed_data/dtu/",
                                                             method, comparison,"protein_coding", "isoformFeatures_part2.tsv"))
SwitchList_part2 <- readRDS(switchlist_part2_path)

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

splicing_enrichment <- splicing_enrichment +
  # Change y-axis text size
  theme(
    axis.text.y = element_text(size = 8, angle = 45, hjust = 0.7, vjust = 1),   
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 8),  # Adjust size as needed
    axis.title.x = element_text(size = 16)  # Also adjust y-axis title if desired
  ) +
  # Change the color scale to use light blue instead of red
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "lightgray"),  # Light blue for colorblind-friendly
    name = "FDR < 0.05",
    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")
  ) +
  # Wrap y-axis labels to multiple lines (alternative to angle)
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20))
# Display the plot
splicing_enrichment
dev.off()

splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = T,
  onlySigIsoforms = T,
  countGenes = F,
  plot = F
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
p <- p +
  # Change y-axis text size
  theme(
    axis.text.y = element_text(size = 8, angle = 30, hjust = 0.5, vjust = 1),   
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 8),  # Adjust size as needed
    axis.title.x = element_text(size = 16)  # Also adjust y-axis title if desired
  ) +
  # Change the color scale to use light blue instead of red
  scale_color_manual(
    values = c("TRUE" = "black", "FALSE" = "lightgray"),  # Light blue for colorblind-friendly
    name = "FDR < 0.05",
    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")
  ) +
  # Wrap y-axis labels to multiple lines (alternative to angle)
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25))
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
  countGenes = F,
  plot = F
)
write_tsv(consequences, file = file.path(plots_dir,
                                         "Consequence_Enrichment.tsv"))


write_tsv(SwitchList_part2$switchConsequence,
           file = file.path(code_dir,"processed_data/dtu/",
                                      method, comparison, "protein_coding", "switchConsequence.tsv")
)


