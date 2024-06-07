library(edgeR)


bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation"

cpm <- read.table(file.path(bambu_dir, "CPM_transcript.txt"),
                  header = TRUE)
colnames(cpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(cpm))
colnames(cpm)[1] <- "isoform_id"
cpm <- cpm[, -2]
head(cpm)

targets  <- data.frame(Sample = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2") ,
                        Group = c( "FT", "FT", "RGC", "RGC"),
                       Replicate = c(1, 2, 1, 2),
                        stringsAsFactors = FALSE)


y <- DGEList(counts = cpm[2:5],
             samples = targets$Sample,
             group = targets$Group,
             genes = cpm[1])

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)

pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotMDS.pdf")
p <- plotMDS(y,col = c(1:2)[y$samples$group],labels = y$samples$Sample,xlim = c(-4,4))
print(p)
dev.off()

design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotBCV.pdf")
p <- plotBCV(y)
print(p)
dev.off()

fit <- glmQLFit(y, design, robust=TRUE)
pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotQLDisp.pdf")
p <- plotQLDisp(fit)
print(p)
dev.off()

contr <- makeContrasts(FT - RGC, levels=design)
qlf <- glmQLFTest(fit, contrast=contr)

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)

# The top set of most significant differentially expressed transcripts can be examined with top
# Tags, with a positive log-fold change representing up-regulation in expression levels in FT
# over RGC.

tt <- topTags(qlf,n = Inf)
nrow(tt) #10261
head(tt)

tt$table$isoform_id <- gsub("\\..*", "", tt$table$isoform_id)


###add gene names and ids

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
listAttributes(mart)[20:30,]

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_gene_id",
                "ensembl_transcript_id",
                "external_gene_name",
                "transcript_biotype"),
  filter="ensembl_transcript_id",
  values=tt$table$isoform_id,
  uniqueRows=TRUE)

colnames(annotLookup) <- c("gene_id", "isoform_id", "gene_name", "transcript_biotype")
tt$table <- merge(tt$table, annotLookup,
                                    by="isoform_id", all.x=TRUE)

unique(tt$table$transcript_biotype)

#order by FDR
tt$table <- tt$table[order(tt$table$FDR),]

write.table(tt$table, file = "/users/sparthib/retina_lrs/processed_data/dtu/edgeR/bambu/FT_vs_RGC/DTEs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# sessionInfo()
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices datasets  utils     methods   base
# 
# other attached packages:
#   [1] biomaRt_2.58.2 edgeR_4.0.14   limma_3.58.1
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.3          utf8_1.2.4              generics_0.1.3
# [4] xml2_1.3.6              bitops_1.0-7            RSQLite_2.3.5
# [7] stringi_1.8.3           lattice_0.22-5          hms_1.1.3
# [10] digest_0.6.34           magrittr_2.0.3          grid_4.3.2
# [13] fastmap_1.1.1           blob_1.2.4              progress_1.2.3
# [16] AnnotationDbi_1.64.1    GenomeInfoDb_1.38.5     DBI_1.2.1
# [19] httr_1.4.7              purrr_1.0.2             fansi_1.0.6
# [22] XML_3.99-0.16.1         Biostrings_2.70.2       cli_3.6.2
# [25] rlang_1.1.3             crayon_1.5.2            dbplyr_2.4.0
# [28] XVector_0.42.0          Biobase_2.62.0          bit64_4.0.5
# [31] splines_4.3.2           withr_3.0.0             cachem_1.0.8
# [34] tools_4.3.2             memoise_2.0.1           dplyr_1.1.4
# [37] filelock_1.0.3          locfit_1.5-9.8          GenomeInfoDbData_1.2.11
# [40] BiocGenerics_0.48.1     curl_5.2.0              vctrs_0.6.5
# [43] R6_2.5.1                png_0.1-8               stats4_4.3.2
# [46] lifecycle_1.0.4         BiocFileCache_2.10.1    zlibbioc_1.48.0
# [49] KEGGREST_1.42.0         stringr_1.5.1           S4Vectors_0.40.2
# [52] IRanges_2.36.0          bit_4.0.5               pkgconfig_2.0.3
# [55] pillar_1.9.0            glue_1.7.0              Rcpp_1.0.12
# [58] statmod_1.5.0           tidyselect_1.2.0        tibble_3.2.1
# [61] compiler_4.3.2          prettyunits_1.2.0       RCurl_1.98-1.14
