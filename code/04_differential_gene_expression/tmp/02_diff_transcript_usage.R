library("ggplot2")
library("DRIMSeq")
# library(DEXSeq)
library("edgeR")
library("here")
library("readr")
library("dplyr")
library(data.table)
library(purrr)
library(tidyr)
# REFERENCE_GTF="/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf"
# CONFIG="/users/sparthib/retina_lrs/raw_data/data_paths.config"

tpm_counts_samples <- list.files("/dcs04/hicks/data/sparthib/casey/diff_expression_data/transcript_lengths",
                      full.names=TRUE)

tbl_fread <- 
  list.files("/dcs04/hicks/data/sparthib/casey/diff_expression_data/transcript_lengths",
             full.names=TRUE) |> map_df(~fread(.))

# library(GenomicFeatures)
# gtf <- "/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf"
# txdb <- makeTxDbFromGFF(gtf)
# 
# #get transcript features 
# transcript_lengths <- transcriptLengths(txdb, with.cds_len=FALSE,
#                   with.utr5_len=FALSE, with.utr3_len=FALSE)
# 
# transcript_lengths <- as.data.frame(transcript_lengths[2:5])
# 
# write.table(transcript_lengths,
#             file = "/dcs04/hicks/data/sparthib/transcript_lengths.tsv",
#             row.names=FALSE, sep="\t", quote = FALSE)



tbl_fread_tpm <- tbl_fread |> select(`#feature_id`, gene_id, TPM, sample_name)
tbl_fread_count <- tbl_fread |> select(`#feature_id`, gene_id, count, sample_name)
split_data_count <- split(tbl_fread_count, f = tbl_fread$sample_name)

countdata_wide <- split_data_count |> reduce(inner_join, by="#feature_id")

countdata_wide <- countdata_wide |> select(`#feature_id`, gene_id.x, count.x,
                                           count.y, count.x.x, count.y.y)
colnames(countdata_wide) <- c("feature_id", "gene_id", "DG-WT-hRGC", 
                              "EP1-BRN3B-RO",
                              "H9-BRN3B-RO",
                              "hRGC")

countdata_wide <- as.data.frame(countdata_wide)
sample_metadata <- data.frame(
  sample_id = c("DG-WT-hRGC", 
                "EP1-BRN3B-RO",
                "H9-BRN3B-RO",
                "hRGC"),
  group = c("RGC", "RO", "RO", "RGC"))
# levels(sample_metadata$group)
# NULL

d <- dmDSdata(counts = countdata_wide, samples = sample_metadata)

head(counts(d), 3)

plotData(d)


###filtering gene expression 

d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_expr = 2,
              min_gene_expr = 2, min_feature_expr = 2)

## Calculate precision
d <- dmPrecision(d, design = design_full)
head(mean_expression(d), 3)
head(genewise_precision(d))

## To make the analysis reproducible
set.seed(123)

## Create the design matrix
d <- dmFit(d, design = design_full, verbose = 1)

head(proportions(d))
head(coefficients(d))
head(coefficients(d), level = "feature")

d <- dmTest(d, coef = "groupRO", verbose = 1)
design(d)

head(results(d), 3)



design_null <- model.matrix(~ 1, data = samples(d))
design_null

d <- dmTest(d, design = design_null)
head(results(d), 3)


contrast <- c(0, 1)
d <- dmTest(d, contrast = contrast)
design(d)

head(results(d), 3)

head(results(d, level = "feature"), 3)
