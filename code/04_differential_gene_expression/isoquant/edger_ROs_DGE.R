library(edgeR)

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/ROs/OUT/"
counts <- read.table(file.path(isoquant_dir,  "OUT.gene_grouped_tpm.tsv"),
                     header = TRUE)


colnames(counts) <- c("gene_id", "EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
                      "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2")



group <- factor(c("RO_D200", "RO_D45", "RO_D100", "RO_D200", "RO_D100", "RO_D45", "RO_D100"))
y <- DGEList(counts=counts, group=group, genes=counts$gene_id,
             samples = c("EP_1_BRN3B_RO",  "EP1_WT_ROs_D45", "EP1_WT_hRO_2",
                         "H9_BRN3B_RO", "H9_BRN3B_hRO_2", "H9_CRX_ROs_D45", "H9_CRX_hRO_2"))

keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)


design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(D200_vs_D100 =  RO_D200 - RO_D100, 
                       D200_vs_D45 = RO_D200 - RO_D45,
                       D100_vs_D45 = RO_D100 - RO_D45,
                       levels=design)

for (i in 1:ncol(contr)){
  
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  # summary(qlf)
  is.de <- decideTests(qlf, p.value=0.05)
  summary(is.de)
  tt <- topTags(qlf,n = Inf)
  nrow(tt) #10261
  head(tt)
  
  tt$table$gene_id <- gsub("\\..*", "", tt$table$gene_id)
  
  
  annotLookup <- getBM(
    mart=mart,
    attributes=c( "ensembl_gene_id",
                  "external_gene_name",
                  "gene_biotype"),
    filter="ensembl_gene_id",
    values=tt$table$gene_id,
    uniqueRows=TRUE)
  
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype")
  tt$table <- merge(tt$table, annotLookup,
                    by="gene_id", all.x=TRUE)
  
  tt$table <- tt$table[order(tt$table$FDR),]
  
  file <- paste0("/users/sparthib/retina_lrs/processed_data/dge/edgeR/isoquant/ROs/", colnames(contr)[i], "_DGEs.tsv")
  write.table(tt$table, file = file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
}
