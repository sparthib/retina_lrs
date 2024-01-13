library(IsoformSwitchAnalyzeR)


salmonQuant <- importIsoformExpression(parentDir = "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/salmon/alignment_mode",
                                       addIsofomIdAsColumn = TRUE)

head(salmonQuant$abundance, 2)
head(salmonQuant$counts, 2)

# nrow = 252669

# Subset the data frame to remove rows where subset columns are all 0

columns_to_check <- c("DG-WT-hRGC", "EP1-BRN3B-RO", "EP1-WT_ROs_D45" , "H9-BRN3B-RO",
                      "H9-CRX_ROs_D45","hRGC","YZ-15T_hRGC","YZ-3T_hRGC")
zero_rows_counts <- rowSums(salmonQuant$counts[columns_to_check] == 0) == length(columns_to_check)
zero_rows_abundance <- rowSums(salmonQuant$abundance[columns_to_check] == 0) == length(columns_to_check)

salmon_counts_filtered <- salmonQuant$counts[!zero_rows_counts, ]
salmon_abundance_filtered <- salmonQuant$abundance[!zero_rows_abundance, ]

# 228128

#design matrix
myDesign  <- data.frame(sampleID = c("DG-WT-hRGC", "EP1-BRN3B-RO", "EP1-WT_ROs_D45" , "H9-BRN3B-RO",
                                     "H9-CRX_ROs_D45","hRGC","YZ-15T_hRGC","YZ-3T_hRGC"),
                        condition = c("RGC", "RO","RO", "RO", 
                                      "RO", "RGC", "RGC","RGC"),
                        stringsAsFactors = FALSE)

aSwitchList <- importRdata(isoformCountMatrix   = salmonQuant$counts,
                           isoformRepExpression = salmonQuant$abundance,
                           designMatrix         = myDesign,
                           isoformExonAnnoation = "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf.gz",
                           isoformNtFasta       = "/dcs04/hicks/data/sparthib/references/transcriptome/GENCODE/gencode.v44.transcripts_short_header.fa",
                           removeNonConvensionalChr = TRUE,
                           ignoreAfterBar = TRUE,
                           ignoreAfterPeriod = TRUE,
                           showProgress = TRUE)



##### sva error #####
# 18797 ( 10.03%) isoforms were removed since they were not expressed in any samples.
# Step 3 of 10: Fixing StringTie gene annoation problems...
# There were no need to rescue any annotation
# 28453 genes_id were assigned their original gene_id instead of the StringTie gene_id.
# This was only done when it could be done unambiguous.
# Step 4 of 10: Calculating expression estimates from count data...
# Skipped as user supplied expression via the "isoformRepExpression" argument...
# Step 5 of 10: Testing for unwanted effects...
# 
# SVA analysis failed. No unwanted effects were added.
# Step 6 of 10: Batch correcting expression estimates...
# Skipped as no batch effects were found or annoated...
# Step 7 of 10: Extracting data from each condition...
# Step 8 of 10: Making comparisons...
# Step 9 of 10: Making switchAnalyzeRlist object...
# Step 10 of 10: Guestimating differential usage...
# The GUESSTIMATED number of genes with differential isoform usage are:
#   comparison estimated_genes_with_dtu
# 1  RGC vs RO                112 - 186
# Done
# 
# Warning messages:
#   1: In importRdata(isoformCountMatrix = salmonQuant$counts, isoformRepExpression = salmonQuant$abundance,  :
#                       
#  There were estimated unwanted effects in your dataset but the automatic sva run failed.
#  We highly reccomend you run sva yourself, add the nessesary surrogate variables
# as extra columns in the "designMatrix" and re-run this function
#                     
# 2: In createSwitchAnalyzeRlist(isoformFeatures = isoAnnot, exons = isoformExonStructure,  :
# The gene_ids or isoform_ids were not unique - we identified multiple instances 
#of the same gene_id/isoform_id on different chromosomes.
#To solve this we removed 23 gene_id. 
#Please note there might still be duplicated gene_id located on the same chromosome. 
#Some of these could be due to fusion transcripts which IsoformSwitchAnalyzeR cannot handle.

#The switchAnalyzeRlist now contains all the information imported in separate entries o
#Of the switchAnalyzeRlist object:


######
head(aSwitchList$isoformFeatures,2)
head(aSwitchList$exons,2)


# Removal of single isoform genes is the default setting 
# in preFilter() since these genes, per definition, cannot
#have changes in isoform usage.


exampleSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

NF1_SwitchList_post_filter <- exampleSwitchListFiltered[[1]] |> dplyr::filter(gene_id == "NF1")


write.csv(NF1_SwitchList_post_filter, "/users/sparthib/retina_lrs/processed_data/dtu/NF1_pre_DEXSeq_Jan13.csv")
###### why DEXSeq #######
# DEXSeq : Using DEXSeq (state-of-the art) to test for differential isoform usage. 
# DEXSeq was originally designed for testing for differential exon usage but have 
# recently been shown to perform exceptionally well for differential isoform usage. 
# The default and recommend in IsoformSwitchAnalyzeR

# Two major challenges in testing differential isoform usage have been controlling 
# false discovery rates (FDR) and applying effect size cutoffs in experimental setups
# with confounding effects. 
# Recent studies such as Love at al highlights DEXSeq (developed by Anders et al., see 
# What To Cite — please remember to cite it) as being a good solution as it controls 
# FDR quite well. We have therefore implemented a DEXSeq based test as the default in IsoformSwitchAnalyzeR. 
# This test furthermore utilizes limma to produce effect sizes corrected for
# confounding effects.

##### use DEXSeq #######

exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = exampleSwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  alpha = 0.05 #FDR cutoff = 0.05
  
)

as.data.frame(table(exampleSwitchListAnalyzed[[1]]$iso_biotype))
# > nrow(exampleSwitchListAnalyzed[[1]])
# 
# [1] 5466



exampleSwitchListAnalyzed[[1]]

#write.csv(exampleSwitchListAnalyzed[[1]], "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/SwitchListAnalyzed_fdr_0.05_first_4.csv")




####### output #####
# A 'switchAnalyzeRlist' where the following have been modified:
#   
#   • '1': Two columns, 'isoform_switch_q_value' and
# 'gene_switch_q_value' in the 'isoformFeatures' entry have
# overwritten with the result of the test.
# 
# • '2': A 'data.frame' containing the details of the analysis
# have been added (called 'isoformSwitchAnalysis').
# 
# The data.frame added have one row per isoform per comparison of
# condition and contains the following columns:
#   
#   • 'iso_ref' : A unique reference to a specific isoform in a
# specific comparison of conditions. Enables easy handles to
# integrate data from all the parts of a 'switchAnalyzeRlist'.
# 
# • 'gene_ref' : A unique reference to a specific gene in a
# specific comparison of conditions. Enables easy handles to
# integrate data from all the parts of a 'switchAnalyzeRlist'.
# 
# • 'isoform_id': The name of the isoform analyzed. Matches the
# 'isoform_id' entry in the 'isoformFeatures' entry of the
# switchAnalyzeRlist
# 
# • 'condition_1': Condition 1 - the condition used as baseline.
# 
# • 'condition_2': Condition 2.
# 
# • 'dIF': The difference in IF values (IF2-IF1) - potentially
# corrected for confounding effects.
# 
# • 'pvalue': Isoform level P-values.
# 
# • 'padj': Isoform level False Discovery Rte (FDR) corrected
# P-values (q-values).
# 
# • 'IF1': Mean isoform fraction in condition 1 - potentially
# corrected for confounding effects.
# 
# • 'IF2': Mean isoform fraction in condition 2 - potentially
# corrected for confounding effects.



######## Plot TOP 50 DTU genes ########
extractSwitchSummary(exampleSwitchListAnalyzed)

top_genes <- extractTopSwitches(exampleSwitchListAnalyzed, n=50)$gene_name

pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/switch_Jan13.pdf")

# /users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR

for(gene_name in top_genes){
  plot <- switchPlot(
    exampleSwitchListAnalyzed,
    gene=gene_name,
    condition1 = 'RGC',
    condition2 = 'RO',
    plotTopology=FALSE
  ) 
  
}
print(plot)
dev.off()


###NF1 plot 
pdf("/users/sparthib/retina_lrs/plots/de/switch_analyzer/NF1_switch_Jan13.pdf")
NF1_plot <- switchPlot(
  exampleSwitchListFiltered,
  gene="NF1",
  condition1 = 'RGC',
  condition2 = 'RO',
  plotTopology=FALSE
)   
print(NF1_plot)
dev.off()


###### DGE #####

SwitchList_gene_expression <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = NULL,
  isoformExpressionCutoff = NULL,
  removeSingleIsoformGenes = FALSE
)

SwitchList_gene_expression[[1]]$gene_name <- SwitchList_gene_expression[[1]]$gene_id

geneCountMatrix <- extractGeneExpression(
  SwitchList_gene_expression,
  extractCounts = TRUE # set to FALSE for abundances
)

write.csv(geneCountMatrix,"/users/sparthib/retina_lrs/processed_data/de/gene_counts_matrix_Jan13.csv")

