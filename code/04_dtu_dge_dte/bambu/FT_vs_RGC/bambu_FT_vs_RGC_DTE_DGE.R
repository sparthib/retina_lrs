library(edgeR)
library(readr)

method <- "bambu"
comparison <- "FT_vs_RGC"
groups <- c( "FT", "FT", "RGC", "RGC")
##### DTE ######
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered_by_counts_and_biotype")
counts <- file.path(matrix_dir, "filtered_isoform_counts.RDS") 
counts <- readRDS(counts)
nrow(counts)

isoformFeatures <- read_tsv(file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                                      method, comparison, "protein_coding", "isoformFeatures.tsv"))

rownames(counts) <- gsub("\\..*", "", rownames(counts))
counts <- counts[rownames(counts) %in% isoformFeatures$isoform_id,]
nrow(counts)


y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))

y <- normLibSizes(y)

design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(RGC - FT, levels=design)
qlf <- glmQLFTest(fit, contrast=contr)

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)

# > summary(is.de)
# -1*FT 1*RGC
# Down           394
# NotSig       35662
# Up             370

tt <- topTags(qlf,n = Inf)
nrow(tt) 
head(tt)
colnames(tt$table) <- c("isoform_id", "logFC", "logCPM", "F" , "PValue", "FDR")

tt$table$isoform_id <- gsub("\\..*", "", tt$table$isoform_id)
tt$table <- tt$table[order(tt$table$FDR),]
tt$table$condition_1 <- names(contr[,1])[contr[,1] == -1]
tt$table$condition_2 <- names(contr[,1])[contr[,1] == 1]

output_path <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison, "protein_coding","DTE_table.tsv" )

write.table(tt$table, file = output_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

##### DGE ######
matrix_dir <- file.path("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/counts_matrices/",
                        method, comparison, "filtered_by_counts_and_biotype")
counts <- file.path(matrix_dir, "filtered_gene_counts.RDS") 
counts <- readRDS(counts)

y <- DGEList(counts = counts,
             samples = colnames(counts),
             group = groups,
             genes = rownames(counts))
y <- normLibSizes(y)


design <- model.matrix(~ 0 + group,data = y$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
fit <- glmQLFit(y, design, robust=TRUE)


contr <- makeContrasts(RGC - FT, levels=design)
qlf <- glmQLFTest(fit, contrast=contr)

is.de <- decideTests(qlf, p.value=0.05)
summary(is.de)

# > summary(is.de)
# -1*FT 1*RGC
# Down          1462
# NotSig       12443
# Up            1196

tt <- topTags(qlf,n = Inf)
nrow(tt) #10261
head(tt)
colnames(tt$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")

tt$table$gene_id <- gsub("\\..*", "", tt$table$gene_id)
tt$table <- tt$table[order(tt$table$FDR),]
tt$table$condition_1 <- names(contr[,1])[contr[,1] == -1]
tt$table$condition_2 <- names(contr[,1])[contr[,1] == 1]

output_path <- file.path("/users/sparthib/retina_lrs/processed_data/dtu/",
                         method, comparison,"protein_coding", "DGE_table.tsv" )

write.table(tt$table, file = output_path,
            sep = "\t", quote = FALSE, row.names = FALSE)



