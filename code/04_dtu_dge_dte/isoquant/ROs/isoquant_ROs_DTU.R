library(IsoformSwitchAnalyzeR)
library(readr)
library(sessioninfo)
library(here)
library(biomaRt)

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples"
sqanti_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc"


method <- "Isoquant"
comparison <- "ROs"
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison)
counts <- file.path(matrix_dir, "isoform_counts.RDS") 
counts <- readRDS(counts)
cpm <- file.path(matrix_dir, "isoform_cpm.RDS")
cpm <- readRDS(cpm)

myDesign  <- data.frame(sampleID = colnames(counts) ,
                        condition = c("C_RO_D45", "C_RO_D45", "B_RO_D100","B_RO_D100","B_RO_D100", "A_RO_D200","A_RO_D200" ),
                        stringsAsFactors = FALSE)


if(!dir.exists(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison, "rds"))){
  dir.create(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "rds"), 
             showWarnings = T, recursive = T)
}

rdata_path = file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                       method, comparison, "rds", "SwitchList.rds")


if(!file.exists(rdata_path)){
  SwitchList <- importRdata(isoformCountMatrix   = counts,
                            isoformRepExpression = cpm,
                            designMatrix         = myDesign,
                            isoformExonAnnoation = file.path(isoquant_dir, "OUT",  "OUT.transcript_models.gtf"),
                            isoformNtFasta       = file.path(sqanti_dir, "all_samples_corrected.fasta"),
                            removeNonConvensionalChr = TRUE,
                            ignoreAfterBar = TRUE,
                            ignoreAfterPeriod = FALSE,
                            showProgress = TRUE)
  #if cds gtf doesn't exist, convert gff to gtf
  if(!file.exists(file.path(sqanti_dir, "all_samples_corrected.gtf.cds.gtf"))) {
    gff <- rtracklayer::import(file.path(sqanti_dir, "all_samples_corrected.gtf.cds.gff"))
    rtracklayer::export(gff, file.path(sqanti_dir, "all_samples_corrected.gtf.cds.gtf"), "gtf")
    rm(gff)
  }
  
  SwitchList <- addORFfromGTF(
    switchAnalyzeRlist     = SwitchList,
    pathToGTF              = file.path(sqanti_dir, "all_samples_corrected.gtf.cds.gtf"),
    overwriteExistingORF=TRUE
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
  
}else { 
  SwitchListFiltered <- readRDS(rdata_path)
}

### load DTEs ###
D200_vs_D100_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D200_vs_D100_DTEs.tsv"))
D200_vs_D45_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D200_vs_D45_DTEs.tsv"))
D100_vs_D45_DTE_table <- read_tsv(  file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DTE" , "D100_vs_D45_DTEs.tsv"))

DTE_table <- rbind(D200_vs_D100_DTE_table, D200_vs_D45_DTE_table, D100_vs_D45_DTE_table)
write_tsv(DTE_table, file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                                      method, comparison, "DTE_table.tsv"))

### load DGEs ###
D200_vs_D100_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D200_vs_D100_DGEs.tsv"))
D200_vs_D45_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D200_vs_D45_DGEs.tsv"))
D100_vs_D45_DGE_table <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/", method, comparison, "DGE" , "D100_vs_D45_DGEs.tsv"))

DGE_table <- rbind( D200_vs_D100_DGE_table, D200_vs_D45_DGE_table, D100_vs_D45_DGE_table)
write_tsv(DGE_table, file = file.path("/users/sparthib/retina_lrs/processed_data/dtu/", 
                                      method, comparison, "DGE_table.tsv"))


dtu_rdata_path  = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/DexSeqDTUDGESwitchList.rds"
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
  # 1  RO_D100 vs RO_D45     2076     1588  1495
  # 2 RO_D100 vs RO_D200     1684     1389  1290
  # 3  RO_D200 vs RO_D45     2490     2055  1775
  # 4           Combined     4804     4417  3188
  
  SwitchList_part1 <- analyzeORF(SwitchList_part1, genomeObject = Hsapiens)
  
  SwitchList_part1$aaSequence = NULL
  
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/fastas/",
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