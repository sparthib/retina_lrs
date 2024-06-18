library(edgeR)

isoquant_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT"
counts <- read.table(file.path(isoquant_dir, "OUT.gene_grouped_tpm.tsv"),
                     header = FALSE)

#remove "_primary_over_30_chr_only_sorted" in column names 
colnames(counts) <- c("gene_id","H9-FT_1", "H9-FT_2", "H9-hRGC_1","H9-hRGC_2")

group <- factor(c("FT", "FT", "RGC", "RGC"))
y <- DGEList(counts=counts, group=group, genes=counts$gene_id,
             samples = c("H9-FT_1", "H9-FT_2", "H9-hRGC_1","H9-hRGC_2"))
keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)

y$samples
y <- estimateDisp(y)

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


tt <- topTags(qlf,n = Inf)

tt$table$gene_id <- gsub("\\..*", "", tt$table$gene_id)


###add gene names and ids
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

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

unique(tt$table$gene_biotype)

#order by FDR
tt$table <- tt$table[order(tt$table$FDR),]

write.table(tt$table, file = "/users/sparthib/retina_lrs/processed_data/dge/edgeR/isoquant/FT_vs_RGC/DGEs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
