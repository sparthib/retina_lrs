library(edgeR)


input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT"

cpm <- read.table(file.path(input_dir, "OUT.transcript_model_grouped_tpm.tsv"),
                  header = FALSE)


colnames(cpm) <- c("isoform_id", "EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2", "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2")


targets  <- data.frame(Sample = c("EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2", "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2") ,
                       Group = c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"),
                       Replicate = c(1, 1, 1, 2, 2, 2 ,3),
                       stringsAsFactors = FALSE)

y <- DGEList(counts = cpm[2:8],
             samples = targets$Sample,
             group = targets$Group,
             genes = cpm[1])

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)


output_plot_dir <- "/users/sparthib/retina_lrs/plots/de/edgeR/isoquant/ROs"
pdf(paste0(output_plot_dir, "/plotMDS.pdf"))
p <- plotMDS(y,col = c(1:2)[y$samples$group],labels = y$samples$Sample,xlim = c(-4,4))
print(p)
dev.off()

design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion


pdf(paste0(output_plot_dir, "/plotBCV.pdf"))
p <- plotBCV(y)
print(p)
dev.off()

fit <- glmQLFit(y, design, robust=TRUE)
pdf(paste0(output_plot_dir, "/plotQLDisp.pdf"))
p <- plotQLDisp(fit)
print(p)
dev.off()

contr <- makeContrasts(D200_vs_D100 =  RO_D200 - RO_D100, 
                       D200_vs_D45 = RO_D200 - RO_D45,
                       D100_vs_D45 = RO_D100 - RO_D45,
                       levels=design)

# require("biomaRt")
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
# listAttributes(mart)[20:30,]

# > i = 1
# >   qlf <- glmQLFTest(fit, contrast = contr[,i])
# >   is.de <- decideTests(qlf, p.value=0.05)
# >   summary(is.de)
# -1*RO_D100 1*RO_D200
# Down                    560
# NotSig                11364
# Up                      700
# > i = 2
# >   qlf <- glmQLFTest(fit, contrast = contr[,i])
# >   is.de <- decideTests(qlf, p.value=0.05)
# >   summary(is.de)
# 1*RO_D200 -1*RO_D45
# Down                  1345
# NotSig                9912
# Up                    1367
# > i = 3
# >   qlf <- glmQLFTest(fit, contrast = contr[,i])
# >   is.de <- decideTests(qlf, p.value=0.05)
# >   summary(is.de)
# 1*RO_D100 -1*RO_D45
# Down                   754
# NotSig               11107
# Up                     763

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

for (i in 1:ncol(contr)){
  
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  # summary(qlf)
  is.de <- decideTests(qlf, p.value=0.05)
  summary(is.de)
  tt <- topTags(qlf,n = Inf)
  nrow(tt) #10261
  head(tt)
  
  tt$table$isoform_id <- gsub("\\..*", "", tt$table$isoform_id)
  
  
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
  
  tt$table <- tt$table[order(tt$table$FDR),]
  
  file <- paste0("/users/sparthib/retina_lrs/processed_data/dtu/edgeR/isoquant/ROs/", colnames(contr)[i], "_DTEs.tsv")
  write.table(tt$table, file = file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  
}

