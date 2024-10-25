library(IsoformSwitchAnalyzeR)
library(readr)
library(sessioninfo)
library(here)
library(biomaRt)

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples"
#read counts_transcript.txt table from isoquant dir
tpm <- read.table(file.path(isoquant_dir, "/OUT/OUT.transcript_model_grouped_tpm.tsv"),
                             header = FALSE)
sqanti_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/high_quality/all_samples/sqanti3_qc"
#feature_id     EP1-BRN3B-RO    EP1-WT_ROs_D45    EP1-WT_hRO_2    H9-BRN3B-RO    H9-BRN3B_hRO_2    H9-CRX_ROs_D45    H9-CRX_hRO_2    H9-FT_1    H9-FT_2    H9-hRGC_1    H9-hRGC_2

counts <- read.table(file.path(isoquant_dir, "OUT","OUT.transcript_model_grouped_counts.tsv"),
                      header = FALSE)

sample_names <- c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                  "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2", "H9-FT_1", 
                  "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2")

colnames(counts) <- c("isoform_id", sample_names)
colnames(tpm) <- c("isoform_id", sample_names)

#for ROs 
RO_counts <- counts[,1:8]
RO_tpm <- tpm[,1:8]

myDesign  <- data.frame(sampleID = c("EP1-BRN3B-RO", "EP1-WT_ROs_D45", "EP1-WT_hRO_2", "H9-BRN3B-RO", 
                                     "H9-BRN3B_hRO_2", "H9-CRX_ROs_D45", "H9-CRX_hRO_2") ,
                        condition = c("RO_D200", "RO_D45", "RO_D100", "RO_D200", 
                                      "RO_D100", "RO_D45", "RO_D100"),
                        stringsAsFactors = FALSE)

dir.create(here("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds"), showWarnings = T,
           recursive = T)
rdata_path = here("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/SwitchList.rds")

### isoformExonAnnotation - annotation of only the transcripts observed in your sample.

if(!file.exists(rdata_path)){
  SwitchList <- importRdata(isoformCountMatrix   = RO_counts,
                            isoformRepExpression = RO_tpm,
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
  biomartCacheClear()
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
    geneExpressionCutoff       = 3,     
    isoformExpressionCutoff    = 1,     
    removeSingleIsoformGenes   = TRUE  # default
  )
  
  saveRDS(SwitchListFiltered, file = rdata_path)
  #save DEXSeq switchlist
  
}else { 
  SwitchListFiltered <- readRDS(rdata_path)
}

# if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/DTE_table.tsv")){
#   
dtu_rdata_path  = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/DexSeqDTUDGESwitchList.rds"
dir.create("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/fastas/", showWarnings = T, recursive = T)

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

###### DGE DTE ######
if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE_DTU_DTE.tsv")){
  
  
  ### load DTEs ###
  D200_vs_D100_DTE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DTE/D200_vs_D100_DTEs.tsv")
  D200_vs_D45_DTE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DTE/D200_vs_D45_DTEs.tsv")
  D100_vs_D45_DTE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DTE/D100_vs_D45_DTEs.tsv")
  
  #rbind 
  DTE_table <- rbind(D200_vs_D100_DTE_table, D200_vs_D45_DTE_table, D100_vs_D45_DTE_table)
  write_tsv(DTE_table, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DTE_table.tsv")
  
  ### load DGEs ###
  D200_vs_D100_DGE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE/D200_vs_D100_DGEs.tsv")
  D200_vs_D45_DGE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE/D200_vs_D45_DGEs.tsv")
  D100_vs_D45_DGE_table <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE/D100_vs_D45_DGEs.tsv")
  
  #rbind
  DGE_table <- rbind( D200_vs_D100_DGE_table, D200_vs_D45_DGE_table, D100_vs_D45_DGE_table)
  write_tsv(DGE_table, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE_table.tsv")
  
  #check for duplicates 
  SwitchList_part1$isoformFeatures <- SwitchList_part1$isoformFeatures |> distinct() 
  write_tsv(SwitchList_part1$isoformFeatures, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/isoformFeatures.tsv")
  write_tsv(SwitchList_part1$isoformSwitchAnalysis, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/isoformSwitchAnalysis.tsv")
  
  
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
  DTE_table <- DTE_table |> dplyr::select(c(isoform_id, logFC, PValue, FDR,
                                            condition_1, condition_2))
  
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
  
  DGE_table <- DGE_table |> dplyr::select(c(gene_id, logFC, PValue, FDR,
                                            condition_1, condition_2))
  

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
  
  
  us_mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
  mart <- useDataset("hsapiens_gene_ensembl", us_mart)
  
  annotLookup <- getBM(
    mart=mart,
    attributes=c( "ensembl_transcript_id",
                  "external_gene_name",
                  "gene_biotype", "transcript_biotype"),
    filter="ensembl_transcript_id",
    values=DGE_DTU_DTE$isoform_id,
    uniqueRows=TRUE)
  
  colnames(annotLookup) <- c("isoform_id", "gene_biotype", "transcript_biotype")
  
  
  DGE_DTU_DTE <- merge(DGE_DTU_DTE, annotLookup,
                       by="isoform_id", all.x=TRUE)
  
  
  
  nrow(DGE_DTU_DTE)
  
  #COUNT NUMBER OF NAs under DTE
  sum(is.na(DGE_DTU_DTE$DTE))
  sum(is.na(DGE_DTU_DTE$DTU))
  sum(is.na(DGE_DTU_DTE$DGE))
  
  #get genes that have DGE as NA 
  DGE_DTU_DTE[is.na(DGE_DTU_DTE$DGE),] |> head(5)
  
  table(DGE_DTU_DTE$DGE, useNA = "always")
  table(DGE_DTU_DTE$DTE, useNA = "always")
  table(DGE_DTU_DTE$DTU, useNA = "always")

  
  write_tsv( DGE_DTU_DTE, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE_DTU_DTE.tsv")
} else {
  DGE_DTU_DTE <- readr::read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/DGE_DTU_DTE.tsv")
}



#### Consequences ####


SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)


saveRDS(SwitchList_part2, file = "/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/SwitchList_part2.rds")

dir.create("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/plots", showWarnings = T, recursive = T)
pdf("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/plots/Splicing_Summary.pdf")
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()

pdf("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/plots/Splicing_Enrichment.pdf",
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
SwitchList_part2 <- readRDS("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/rds/SwitchList_part2.rds")

dir.create("/users/sparthib/retina_lrs/processed_data/dtu/Isoquant/ROs/external_protein_analyses/", showWarnings = T, recursive = T)

SwitchList_part2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_part2, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/external_protein_analyses/CPC2_output.txt",
  pathToPFAMresultFile      = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/external_protein_analyses/pfam_results.txt",
  pathToSignalPresultFile   = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/external_protein_analyses/prediction_results.txt",
  outputPlots               = TRUE,
  pathToOutput              = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/external_protein_analyses/switchplots_with_consequences",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  )
  
)

saveRDS(SwitchList_part2, file = "./processed_data/dtu/DTU_gandall/bambu/ROs/rds/SwitchList_part2.rds")


# The number of isoform switches with functional consequences identified were:
#   Comparison nrIsoforms nrSwitches nrGenes
# 1  RO_D100 vs RO_D45       1299       1168     898
# 2 RO_D100 vs RO_D200       1044       1027     773
# 3  RO_D200 vs RO_D45       1545       1473    1046
# 4           Combined       3050       3232    1977

pdf("./processed_data/dtu/DTU_gandall/bambu/ROs/plots/Consequence_Enrichment.pdf",
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
  returnResult = FALSE, # if TRUE returns a data.frame with the summary statistics
  countGenes = F
)
print(p)
dev.off()

# to convert pfam output to right format
# write.table(pfam_results, file = "/Users/sparthib/Documents/retina_lrs/processed_data/dtu/DTU_gandall/bambu/ROs/external_protein_analyses/pfam_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
