library(edgeR)

input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT"
output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/edgeR/isoquant/FT_vs_RGC"

tpm <- read.table(file.path(input_dir, "OUT.transcript_model_grouped_tpm.tsv"),
                  header = FALSE, sep = "\t")
head(tpm)
colnames(tpm) <- c( "isoform_id" , "H9.FT_1_primary_over_30_chr_only_sorted", "H9.FT_2_primary_over_30_chr_only_sorted", "H9.hRGC_1_primary_over_30_chr_only_sorted","H9.hRGC_2_primary_over_30_chr_only_sorted")
#remove "_primary_over_30_chr_only_sorted" in column names 
colnames(tpm) <- gsub("_primary_over_30_chr_only_sorted", "", colnames(tpm))
colnames(tpm)[1] <- "isoform_id"

head(tpm)


targets  <- data.frame(Sample = c("H9.FT_1","H9.FT_2", "H9.hRGC_1", "H9.hRGC_2") ,
                       Group = c( "FT", "FT", "RGC", "RGC"),
                       Replicate = c(1, 2, 1, 2),
                       stringsAsFactors = FALSE)


y <- DGEList(counts = tpm[2:5],
             samples = targets$Sample,
             group = targets$Group,
             genes = tpm[1])

keep <- filterByExpr(y, )
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)

output_plot_dir <- "/users/sparthib/retina_lrs/plots/de/edgeR/isoquant/FT_vs_RGC/"
pdf(paste0(output_plot_dir, "plotMDS.pdf"))
p <- plotMDS(y,col = c(1:2)[y$samples$group],labels = y$samples$Sample,xlim = c(-4,4))
print(p)
dev.off()

design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

pdf(paste0(output_plot_dir, "plotBCV.pdf"))
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

1*FT -1*RGC

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

output_data_dir <- "/users/sparthib/retina_lrs/processed_data/dtu/edgeR/isoquant/FT_vs_RGC/"
write.table(tt$table, file = paste0(output_data_dir, "DTEs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

