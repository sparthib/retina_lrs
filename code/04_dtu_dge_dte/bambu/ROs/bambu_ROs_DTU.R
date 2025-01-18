library(IsoformSwitchAnalyzeR)
library(tximeta)
library(readr)
library(sessioninfo)
library(here)

method <- "bambu"
comparison <- "ROs"
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered")

counts <- file.path(matrix_dir, "isoform_counts.RDS") 
counts <- readRDS(counts)
cpm <- file.path(matrix_dir, "isoform_cpm.RDS")
cpm <- readRDS(cpm)


myDesign  <- data.frame(sampleID = colnames(counts) ,
                        condition = c("Stage_1", "Stage_1", "Stage_2","Stage_2","Stage_2", "Stage_3","Stage_3"),
                        stringsAsFactors = FALSE)

if(!dir.exists(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison, "rds"))){
  dir.create(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "rds"), 
             showWarnings = T, recursive = T)
}

rdata_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "rds", "SwitchList.rds")

  
bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/all_samples_extended_annotation_track_reads"

if(!file.exists(rdata_path)){
  
SwitchList <- importRdata(isoformCountMatrix   = counts,
                          isoformRepExpression = cpm ,
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
SwitchList$isoformFeatures <- merge(SwitchList$isoformFeatures, 
                                    annotLookup, by="gene_id", all.x=TRUE)


head(SwitchList$isoformFeatures)
SwitchList$isoformFeatures$gene_name <- SwitchList$isoformFeatures$ensembl_gene_name

SwitchListFiltered <- preFilter(
  switchAnalyzeRlist         = SwitchList,
  geneExpressionCutoff       = 1,     # default
  isoformExpressionCutoff    = 0,     # only keep the isoforms that are expressed more than 10 CPM in atleast one of the conditions
  removeSingleIsoformGenes   = TRUE  # default
)

saveRDS(SwitchListFiltered, file = rdata_path)
SwitchListFiltered$isoformFeatures <- SwitchListFiltered$isoformFeatures |> 
  distinct(across(-gene_name), .keep_all = TRUE)

write_tsv(SwitchListFiltered$isoformFeatures,
          file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                           method, comparison, "isoformFeatures.tsv"))
#save DEXSeq switchlist

}else { 
  SwitchListFiltered <- readRDS(rdata_path)
}

### load DTEs ###
D100_vs_D200_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D100_vs_D200_DTEs.tsv"))
D45_vs_D200_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D45_vs_D200_DTEs.tsv"))
D45_vs_D100_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D45_vs_D100_DTEs.tsv"))

DTE_table <- rbind(D100_vs_D200_DTE_table, D45_vs_D200_DTE_table, D45_vs_D100_DTE_table)
write_tsv(DTE_table, file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                                      method, comparison, "DTE_table.tsv"))

### load DGEs ###
D100_vs_D200_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D100_vs_D200_DGEs.tsv"))
D45_vs_D200_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D45_vs_D200_DGEs.tsv"))
D45_vs_D100_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D45_vs_D100_DGEs.tsv"))

DGE_table <- rbind(D100_vs_D200_DGE_table, D45_vs_D200_DGE_table, D45_vs_D100_DGE_table)
write_tsv(DGE_table, file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                                      method, comparison, "DGE_table.tsv"))


dtu_rdata_path <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                            method, comparison, "rds", "DexSeqDTUDGESwitchList.rds")
dir.create(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                     method, comparison, "fastas"), 
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
  # 1  B_RO_D100 vs C_RO_D45     2151     1679  1561
  # 2 A_RO_D200 vs B_RO_D100     1785     1490  1357
  # 3  A_RO_D200 vs C_RO_D45     2508     2102  1788
  # 4               Combined     4983     4572  3306
  
  SwitchList_part1 <- analyzeORF(SwitchList_part1, genomeObject = Hsapiens)
  
  SwitchList_part1$aaSequence = NULL
  
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                   method, comparison, "fastas"),
    extractNTseq       = TRUE, #for CPC2
    extractAAseq       = TRUE,
    removeShortAAseq   = TRUE,
    removeLongAAseq    = TRUE, #FOR PFAM, SignalP
    onlySwitchingGenes = TRUE,
    alsoSplitFastaFile=FALSE
  )
  
  saveRDS(SwitchList_part1, file = dtu_rdata_path)
  saveRDS(SwitchList_part1$isoformFeatures, file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                                            method, comparison, "rds", "isoformFeatures.rds"))
}else{
  SwitchList_part1 <- readRDS(dtu_rdata_path)
}

  #check for duplicates 
  SwitchList_part1$isoformFeatures <- SwitchList_part1$isoformFeatures |> 
    distinct(across(-gene_name), .keep_all = TRUE)
  

#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)


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


#### Consequence Switch Plots ####

external_protein_analyses_dir <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                           method, comparison, "external_protein_analyses")
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

switchlist_part2_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                  method, comparison, "rds", "SwitchList_part2.rds")

saveRDS(SwitchList_part2, file = switchlist_part2_path)


SwitchList_part2$isoformFeatures <- SwitchList_part2$isoformFeatures |> 
  distinct(across(-gene_name), .keep_all = TRUE)

write_tsv(SwitchList_part2$isoformFeatures, file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                                             method, comparison, "isoformFeatures.tsv"))



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
