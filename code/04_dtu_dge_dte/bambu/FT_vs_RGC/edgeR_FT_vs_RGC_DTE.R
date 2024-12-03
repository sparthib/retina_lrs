library(edgeR)


bambu_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/"

counts <- read.table(file.path(bambu_dir, "counts_transcript.txt"),
                  header = TRUE)
colnames(counts) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(counts))
colnames(counts)[1] <- "isoform_id"
counts<- counts[, -2]
head(counts)

myDesign  <- data.frame(Sample = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2") ,
                        Group = c( "FT", "FT", "RGC", "RGC"),
                       Replicate = c(1, 2, 1, 2),
                        stringsAsFactors = FALSE)

y <- DGEList(counts = counts[2:5],
             samples = targets$Sample,
             group = targets$Group,
             genes = counts[1])

keep <- filterByExpr(y, )
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
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
summary(is.de)

# The top set of most significant differentially expressed transcripts can be examined with top
# Tags, with a positive log-fold change representing up-regulation in expression levels in FT
# over RGC.

# 1*FT -1*RGC

tt <- topTags(qlf,n = Inf)
nrow(tt) #10261
head(tt)

tt$table$isoform_id <- gsub("\\..*", "", tt$table$isoform_id)
tt$table <- tt$table[order(tt$table$FDR),]


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

write.table(tt$table, file = "/users/sparthib/retina_lrs/processed_data/dtu/edgeR/bambu/FT_vs_RGC/DTEs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotMDS.pdf")
# p <- plotMDS(y,col = c(1:2)[y$samples$group],labels = y$samples$Sample,xlim = c(-4,4))
# print(p)
# dev.off()
# 
# pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotBCV.pdf")
# p <- plotBCV(y)
# print(p)
# dev.off()
# 
# pdf("/users/sparthib/retina_lrs/plots/de/edgeR/bambu/FT_vs_RGC/plotQLDisp.pdf")
# p <- plotQLDisp(fit)
# print(p)
# dev.off()


