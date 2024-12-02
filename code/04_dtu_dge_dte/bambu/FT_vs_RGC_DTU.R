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

bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation"

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


myDesign  <- data.frame(sampleID = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2") ,
                        condition = c( "FT", "FT", "RGC", "RGC"),
                        stringsAsFactors = FALSE)

rdata_path = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/rds/SwitchList.rds"
if(!file.exists(rdata_path)){ 
  SwitchList <- importRdata(isoformCountMatrix   = counts,
                            isoformRepExpression = cpm,
                            designMatrix         = myDesign,
                            isoformExonAnnoation = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/extended_annotations.gtf",
                            isoformNtFasta       = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/FT_vs_RGC_corrected.fasta",
                            removeNonConvensionalChr = TRUE,
                            ignoreAfterBar = TRUE,
                            ignoreAfterPeriod = FALSE,
                            showProgress = TRUE)
  
  
  #if cds gtf doesn't exist, convert gff to gtf
  if(!file.exists("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/FT_vs_RGC_corrected.gtf.cds.gtf")){
    gff <- rtracklayer::import("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/FT_vs_RGC_corrected.gtf.cds.gff")
    rtracklayer::export(gff, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/FT_vs_RGC_corrected.gtf.cds.gtf",
                        "gtf")
    rm(gff)
  }
  
  SwitchList <- addORFfromGTF(
    switchAnalyzeRlist     = SwitchList,
    pathToGTF              = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/sqanti3_qc/FT_vs_RGC_corrected.gtf.cds.gtf"
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

if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DTE_table.tsv")){
  
  # counts$isoform_id <- gsub("\\..*", "", counts$isoform_id)

  counts_preFilter <- counts |> filter(counts$isoform_id %in% SwitchListFiltered$isoformFeatures$isoform_id)
  
  print(paste0("number of isoforms before filtering: ", nrow(counts)))
  print(paste0("number of isoforms after filtering: ", nrow(counts_preFilter)))
  y <- DGEList(counts = counts_preFilter[2:5],
               samples = myDesign$sampleID,
               group = myDesign$condition,
               genes = counts_preFilter[1])
  y <- normLibSizes(y)
  
  design <- model.matrix(~ 0 + group,data = y$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  design
  y <- estimateDisp(y, design, robust=TRUE)
  y$common.dispersion
  fit <- glmQLFit(y, design, robust=TRUE)
  contr <- makeContrasts(FT - RGC, levels=design)
  qlf <- glmQLFTest(fit, contrast=contr)
  is.de <- decideTests(qlf, p.value=0.05)
  print(summary(is.de))
  
  # print(summary(is.de))
  # 1*FT -1*RGC
  # Down           211
  # NotSig       48443
  # Up             199
  
  tt <- topTags(qlf,n = Inf)
  nrow(tt) #10261
  print("head of topTags")
  head(tt)
  
  tt$table <- tt$table[order(tt$table$FDR),]
  DTE_table <- tt$table
  
  write_tsv(tt$table, file = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DTE_table.tsv")
  
} else {
  DTE_table <- readr::read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DTE_table.tsv")
}



#### DGE ####
if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_table.tsv")){
  counts <- read.table(file.path(bambu_dir, "counts_gene.txt"),
                       header = TRUE)
  #remove "_primary_over_30_chr_only_sorted" in column names 
  colnames(counts) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(counts))
  colnames(counts)[1] <- "gene_id"
  head(counts)
  
  group <- factor(c("FT", "FT", "RGC", "RGC"))
  y_dge <- DGEList(counts=counts, group=group, genes=counts$gene_id,
                   samples = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2"))
  
  #filtering for DTE is based on SwitchList_preFilter, for DGE, just used edgeR function. 
  keep <- filterByExpr(y_dge)
  y_dge <- y_dge[keep, , keep.lib.sizes=FALSE]
  
  y_dge <- normLibSizes(y_dge)
  
  y_dge$samples
  y_dge <- estimateDisp(y_dge)
  
  design <- model.matrix(~ 0 + group,data = y_dge$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  design
  
  
  y_dge <- estimateDisp(y_dge, design, robust=TRUE)
  # y_dge$common.dispersion
  
  fit_dge <- glmQLFit(y_dge, design, robust=TRUE)
  contr <- makeContrasts(FT - RGC, levels=design)
  qlf_dge <- glmQLFTest(fit_dge, contrast=contr)
  
  is.de <- decideTests(qlf_dge, p.value=0.05)
  summary(is.de)
  # >   summary(is.de)
  # 1*FT -1*RGC
  # Down           783
  # NotSig       20006
  # Up             891
  
  tt_dge <- topTags(qlf_dge,n = Inf)
  tt_dge$table$gene_id <- gsub("\\..*", "", tt_dge$table$gene_id)
  DGE_table <- tt_dge$table
  write_tsv(tt_dge$table, file = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_table.tsv")
} else {
  DGE_table <- readr::read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_table.tsv")
}
  

#### DEXSeq ####
dtu_rdata_path  = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DexSeqDTUDGESwitchList.rds"
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
 
  
  SwitchList_part1 <- extractSequence(
    switchAnalyzeRlist = SwitchList_part1,
    pathToOutput       = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/fastas/",
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



  
#similar to table s3 in the paper

if(!file.exists("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")){
  # SwitchList_part1$isoformSwitchAnalysis$isoform_id <- gsub("\\..*", "", SwitchList_part1$isoformSwitchAnalysis$isoform_id)
  DGE_DTU_DTE = SwitchList_part1$isoformFeatures |>
    as_tibble() |>
    dplyr::select(isoform_id, gene_id, gene_name, condition_1, condition_2) |>
    left_join(
      SwitchList_part1$isoformSwitchAnalysis |> dplyr::select(isoform_id, dIF, pvalue, padj)
    ) |>
    dplyr::rename(
      DTU_dIF    = "dIF",
      DTU_pval   = "pvalue",
      DTU_qval   = "padj"
    ) |>
    mutate(
      DTU = DTU_qval < 0.05 # & abs(DTU_dIF) > 0.1
    ) |>
    left_join(
      DTE_table |> dplyr::select(isoform_id, logFC, PValue, FDR)
    ) |>
    dplyr::rename(
      DTE_log2FC = "logFC",
      DTE_pval   = "PValue",
      DTE_qval   = "FDR"
    ) |>
    mutate(
      DTE = DTE_qval < 0.05
    ) |>
    left_join(
      DGE_table |> dplyr::select(gene_id, logFC, PValue, FDR)
    ) |>
    dplyr::rename(
      DGE_log2FC = "logFC",
      DGE_pval   = "PValue",
      DGE_qval   = "FDR"
    ) |>
    mutate(
      DGE = DGE_qval < 0.05
    )
  
  write_tsv( DGE_DTU_DTE, file = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")
} else {
  library(here)
  DGE_DTU_DTE <- readr::read_tsv(here("processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv"))
}



#### VENN DIAGRAM ####

gene_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> 
  summarise(DTE = any(DTE), DGE=any(DGE), DTU=any(DTU))  |> dplyr::select(-gene_id)

dge_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> dplyr::select( -DTE) |> 
  summarise(DGE=any(DGE), DTU=any(DTU))  

dte_overlaps = DGE_DTU_DTE |> group_by(gene_id) |> dplyr::select(-DGE) |> 
  summarise(DTE = any(DTE), DTU=any(DTU)) 
#total 10673 genes after filtering 
dge_contingency_table <- as_tibble(xtabs(~ DGE + DTU, data = dge_overlaps))
dge_contingency_table <- plyr::ddply(dge_contingency_table, ~DGE, transform, Prop = n / sum(n))


dge_bar <- ggplot(dge_contingency_table, aes(x = DGE, y = n, fill = DTU)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "DGE and DTU genes", x = "DGE", y = "Count") + 
  scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
  theme_minimal() +
  geom_text(aes(label = n), 
            position = position_stack(vjust = 0.5), 
            size = 2)

dte_contingency_table <- xtabs(~ DTE + DTU, data = dte_overlaps)
dte_contingency_table <- plyr::ddply(as_tibble(dte_contingency_table), ~DTE, transform, Prop = n / sum(n))

dte_bar <- ggplot(dte_contingency_table, aes(x = DTE, y = n, fill = DTU)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "DTE and DTU genes", x = "DTE", y = "Count") + 
  scale_y_continuous(labels = scales::comma) +   # Display y-axis as percentage
  theme_minimal() +
  geom_text(aes(label = n), 
            position = position_stack(vjust = 0.5), 
            size = 2)
library(patchwork)
p <- dge_bar + dte_bar + plot_annotation(title = 'FT vs RGC')
ggsave(path = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/",
       device = "pdf", plot = p, filename = "DGE_DTU_DTE_barplot.pdf")

# save venn diagram as pdf 
pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/DGE_DTU_DTE_venn.pdf")
fig <- ggVennDiagram(list(DTU = which(gene_overlaps$DTU), 
                          DGE = which(gene_overlaps$DGE),
                          DTE = which(gene_overlaps$DTE))) + 
  scale_fill_gradient(low="grey",high = "red")
print(fig)
dev.off()



#### Consequences ####

SwitchList_part2 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = SwitchList_part1
)


saveRDS(SwitchList_part2, file = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/SwitchList_part2.rds")

SwitchList_part2 <- readRDS("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/SwitchList_part2.rds")

pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/Splicing_Summary.pdf")
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()
pdf("/users/sparthib/retina_lrs//processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/Splicing_Enrichment.pdf")
splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = F ,
  onlySigIsoforms = T,
  countGenes = F
)
print(splicing_enrichment)
dev.off()




#### Consequence Switch Plots ####

SwitchList_part2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList_part2, 
  n                         = 50,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPC2resultFile      = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/external_protein_analyses/CPC2_output.txt",
  pathToPFAMresultFile      = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/external_protein_analyses/pfam_results.txt",
  pathToSignalPresultFile   = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/external_protein_analyses/prediction_results.txt",
  outputPlots               = TRUE,
  pathToOutput              = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/external_protein_analyses/switchplots_with_consequences",
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  ))

saveRDS(SwitchList_part2, file = "/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/rds/SwitchList_part2.rds")


# The number of isoform switches with functional consequences identified were:
#   Comparison nrIsoforms nrSwitches nrGenes
# 1  FT vs RGC        864        899     608

#### Switch Consequence plots  ####
pdf("./processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/Splicing_Summary.pdf")
splicing_summary <- extractSplicingSummary(SwitchList_part2,
                                           splicingToAnalyze = 'all',dIFcutoff = 0.1,
                                           onlySigIsoforms = T,
                                           returnResult = F,
                                           plot = T)  
print(splicing_summary)
dev.off()



pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/Splicing_Enrichment.pdf")
splicing_enrichment <- extractSplicingEnrichment(
  SwitchList_part2,
  returnResult = F ,
  onlySigIsoforms = T,
  countGenes = F
)
print(splicing_enrichment)
dev.off()

pdf("./processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/plots/Consequence_Enrichment.pdf",
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