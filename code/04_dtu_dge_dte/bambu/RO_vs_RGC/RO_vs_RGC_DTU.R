library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(here)

method <- "bambu"
comparison <- "RO_vs_RGC"

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")

matrix_dir <- file.path(data_dir,"06_quantification/counts_matrices/bambu/RO_vs_RGC/filtered_by_counts_and_biotype")

counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)
cpm <- file.path(matrix_dir, "filtered_isoform_cpm.RDS")
cpm <- readRDS(cpm)

nrow(counts)
nrow(cpm)

myDesign  <- data.frame(sampleID = colnames(counts) ,
                        condition = c("Stage_1", "Stage_1", "Stage_2","Stage_2",
                                      "Stage_2", "Stage_3","Stage_3", "RGC","RGC"),
                        stringsAsFactors = FALSE)
comparisions <- data.frame(condition_1 = c("RGC", "RGC", "RGC"),
                           condition_2 = c("Stage_1", "Stage_2", "Stage_3"))
                       
    
if(!dir.exists(file.path(code_dir,"processed_data/dtu/",
                         method, comparison, "protein_coding","rds"))){
  dir.create(file.path(code_dir,"processed_data/dtu/",
                       method, comparison,"protein_coding", "rds"), 
             showWarnings = T, recursive = T)
}

rdata_path = file.path(code_dir,"processed_data/dtu/",
                       method, comparison, "protein_coding", "rds", "SwitchList.rds")


bambu_dir <- file.path(data_dir,"06_quantification/bambu/all_samples_extended_annotation_track_reads")



if(!file.exists(rdata_path)){
  
  SwitchList <- importRdata(isoformCountMatrix   = counts,
                            isoformRepExpression = cpm ,
                            designMatrix         = myDesign,
                            isoformExonAnnoation = paste0(bambu_dir, "/RO_vs_RGC_protein_coding_annotations.gtf"),
                            isoformNtFasta       = paste0(bambu_dir, "/sqanti3_qc/all_samples_corrected.fasta"),
                            removeNonConvensionalChr = TRUE,
                            ignoreAfterBar = TRUE,
                            ignoreAfterPeriod = TRUE,
                            showProgress = TRUE,
                            comparisonsToMake = comparisions)
  #if cds gtf doesn't exist, convert gff to gtf
  if(!file.exists( paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf"))) {
    gff <- rtracklayer::import(paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gff"))
    rtracklayer::export(gff, paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf"),
                        "gtf")
    rm(gff)
  }
  
  SwitchList <- addORFfromGTF(
    switchAnalyzeRlist     = SwitchList,
    pathToGTF              =  paste0( bambu_dir, "/sqanti3_qc/all_samples_corrected.gtf.cds.gtf"),
    ignoreAfterPeriod = TRUE,
    overwriteExistingORF = TRUE
  )
  
  # comparison estimated_genes_with_dtu
  # 1 RGC vs Stage_1                416 - 693
  # 2 RGC vs Stage_2                195 - 324
  # 3 RGC vs Stage_3                255 - 425
  
  SwitchList$isoformFeatures$gene_id <- gsub("\\..*", "",
                                             SwitchList$isoformFeatures$gene_id)
  
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
  SwitchList$isoformFeatures <- merge(SwitchList$isoformFeatures, 
                                      annotLookup, by="gene_id", all.x=TRUE)
  
  
  head(SwitchList$isoformFeatures)
  SwitchList$isoformFeatures$gene_name <- SwitchList$isoformFeatures$ensembl_gene_name
  
  SwitchListFiltered <- preFilter(
    switchAnalyzeRlist         = SwitchList,
    geneExpressionCutoff       = NULL,    
    isoformExpressionCutoff    = NULL,     
    IFcutoff=0,
    removeSingleIsoformGenes   = FALSE)

## 53688 isoforms remaining
#### SwitchListFiltered$isoformFeatures ####
  
  saveRDS(SwitchListFiltered, file = rdata_path)
  SwitchListFiltered$isoformFeatures <- SwitchListFiltered$isoformFeatures |> 
    distinct(across(-gene_name), .keep_all = TRUE)
  
  write_tsv(SwitchListFiltered$isoformFeatures,
            file = file.path(code_dir, "processed_data/dtu/",
                             method, comparison, "protein_coding",  "isoformFeatures.tsv"))
  #save DEXSeq switchlist
  
}else { 
  SwitchListFiltered <- readRDS(rdata_path)
}


### load DTEs ###
D45_vs_RGC_DTE_table <- read_tsv(  file.path(code_dir, "processed_data/dtu/", method, comparison, "protein_coding", "DTE" , "D45_vs_RGC_DTEs.tsv"))
D100_vs_RGC_DTE_table <- read_tsv(  file.path(code_dir, "processed_data/dtu/", method, comparison, "protein_coding", "DTE" , "D100_vs_RGC_DTEs.tsv"))
D200_vs_RGC_DTE_table <- read_tsv(  file.path(code_dir, "processed_data/dtu/", method, comparison, "protein_coding", "DTE" , "D200_vs_RGC_DTEs.tsv"))
DTE_table <- rbind(D45_vs_RGC_DTE_table, D100_vs_RGC_DTE_table, D200_vs_RGC_DTE_table)

write_tsv(DTE_table, file = file.path(code_dir, "processed_data/dtu/", 
                                      method, comparison,"protein_coding", "DTE_table.tsv"))

### load DGEs ###
D45_vs_RGC_DGE_table <- read_tsv(file.path(code_dir, "processed_data/dtu/", method, comparison,"protein_coding",  "DGE" , "D45_vs_RGC_DGEs.tsv"))
D100_vs_RGC_DGE_table <- read_tsv(file.path(code_dir, "processed_data/dtu/", method, comparison,"protein_coding",  "DGE" , "D100_vs_RGC_DGEs.tsv"))
D200_vs_RGC_DGE_table <- read_tsv(file.path(code_dir, "processed_data/dtu/", method, comparison, "protein_coding", "DGE" , "D200_vs_RGC_DGEs.tsv"))

DGE_table <- rbind(D45_vs_RGC_DGE_table, D100_vs_RGC_DGE_table, D200_vs_RGC_DGE_table)
write_tsv(DGE_table, file = file.path(code_dir, "processed_data/dtu/", 
                                      method, comparison, "protein_coding","DGE_table.tsv"))


dtu_rdata_path <- file.path(code_dir, "processed_data/dtu/",
                            method, comparison, "protein_coding", "rds", "DexSeqDTUDGESwitchList.rds")
dir.create(file.path(code_dir, "processed_data/dtu/",
                     method, comparison,"protein_coding", "fastas"), 
           showWarnings = T, recursive = T)

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
  # 1 RGC vs Stage_2     1764     1382  1241
  # 2 RGC vs Stage_1     2813     2215  1850
  # 3 RGC vs Stage_3     2423     1971  1638
  
  SwitchList_part1 <- analyzeORF(SwitchList_part1, genomeObject = Hsapiens)
  
  SwitchList_part1$aaSequence = NULL
  
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = file.path(code_dir, "processed_data/dtu/",
                                   method, comparison, "protein_coding", "fastas"),
    extractNTseq       = TRUE, #for CPC2
    extractAAseq       = TRUE,
    removeShortAAseq   = TRUE,
    removeLongAAseq    = TRUE, #FOR PFAM, SignalP
    onlySwitchingGenes = TRUE,
    alsoSplitFastaFile=FALSE
  )

  #check for duplicates 
  SwitchList_part1$isoformFeatures <- SwitchList_part1$isoformFeatures |> 
    distinct(across(-gene_name), .keep_all = TRUE)
  
  saveRDS(SwitchList_part1, file = dtu_rdata_path)
  saveRDS(SwitchList_part1$isoformFeatures, file.path(code_dir, "processed_data/dtu/",
                                                      method, comparison, "protein_coding" ,"rds", "isoformFeatures.rds"))
  
  write_tsv(SwitchList_part1$isoformFeatures,
            file = file.path(code_dir, "processed_data/dtu/",
                             method, comparison, "protein_coding",  "isoformFeatures.tsv"))
  
}else{
  SwitchList_part1 <- readRDS(dtu_rdata_path)
}


#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)


plots_dir <- file.path(code_dir, "processed_data/dtu/", 
                       method, comparison, "protein_coding", "plots")
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


#### Consequence Switch Plots ####

external_protein_analyses_dir <- file.path(code_dir, "processed_data/dtu/",
                                           method, comparison, "protein_coding","external_protein_analyses")
read_csv(file.path(external_protein_analyses_dir, "pfam_results.csv")) |> 
  write_tsv(file.path(external_protein_analyses_dir, "pfam_results.txt"))

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



switchlist_part2_path = file.path(code_dir, "processed_data/dtu/",
                                  method, comparison, "protein_coding", "rds", "SwitchList_part2.rds")

saveRDS(SwitchList_part2, file = switchlist_part2_path)


SwitchList_part2$isoformFeatures <- SwitchList_part2$isoformFeatures |> 
  distinct(across(-gene_name), .keep_all = TRUE)

write_tsv(SwitchList_part2$isoformFeatures, file = file.path(code_dir, "processed_data/dtu/",
                                                             method, comparison,"protein_coding", "isoformFeatures_part2.tsv"))



pdf(file.path(plots_dir, "Splicing_Summary.pdf"))
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()

splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',
                                           dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = T,
                                           plot = F)  
write_tsv(splicing_summary, file = file.path(plots_dir, 
                                             "Splicing_Summary.tsv"))


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
  countGenes = F,
  plot = F
)
write_tsv(splicing_enrichment, 
          file = file.path(plots_dir, "Splicing_Enrichment.tsv"))

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
  countGenes = F,
  plot = F
)
write_tsv(consequences, file = file.path(plots_dir, "Consequence_Enrichment.tsv"))
