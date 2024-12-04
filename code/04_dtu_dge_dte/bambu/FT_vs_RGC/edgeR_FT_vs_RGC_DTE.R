library(edgeR)

method <- "bambu"
comparison <- "FT_vs_RGC"
groups <- c( "FT", "FT", "RGC", "RGC")
##### DTE ######
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison)
counts <- file.path(matrix_dir, "isoform_counts.RDS") 
counts <- readRDS(counts)

y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

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

# > summary(is.de)
# 1*FT -1*RGC
# Down           378
# NotSig       29240
# Up             399


tt <- topTags(qlf,n = Inf)
nrow(tt) #10261
head(tt)
colnames(tt$table) <- c("isoform_id", "logFC", "logCPM", "F" , "PValue", "FDR")

tt$table$isoform_id <- gsub("\\..*", "", tt$table$isoform_id)
tt$table <- tt$table[order(tt$table$FDR),]


output_path <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison, "DTE_table.tsv" )

write.table(tt$table, file = output_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
##### DGE ######
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison)
counts <- file.path(matrix_dir, "gene_counts.RDS") 
counts <- readRDS(counts)

y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

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



tt <- topTags(qlf,n = Inf)
nrow(tt) #10261
head(tt)
colnames(tt$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")

tt$table$isoform_id <- gsub("\\..*", "", tt$table$gene_id)
tt$table <- tt$table[order(tt$table$FDR),]


output_path <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison, "DGE_table.tsv" )

write.table(tt$table, file = output_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
