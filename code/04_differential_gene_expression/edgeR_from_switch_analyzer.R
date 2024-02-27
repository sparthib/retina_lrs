library(edgeR)

geneCountMatrix <- read.csv("/users/sparthib/retina_lrs/processed_data/de/gene_counts_matrix_Jan13.csv")

group <- c("RGC", "RO","RO", "RO", 
           "RO", "RGC", "RGC","RGC")
y <- DGEList(counts=geneCountMatrix[,3:11],
             group= factor(group) )



y$samples
####filtering y
head(y$counts)
head(cpm(y))
apply(y$counts, 2, sum). #total gene count per sample

# DG.WT.hRGC   EP1.BRN3B.RO EP1.WT_ROs_D45    H9.BRN3B.RO H9.CRX_ROs_D45 
# 17125095       18629344       15864143       16291403       14863666 
# hRGC    YZ.15T_hRGC     YZ.3T_hRGC 
# 24776030       14822000       11639939 


dim(y)
#  55910     8
keep <- rowSums(cpm(y)>10) >= 2
#only keeping a gene if it has a cpm of 10 or greater for at least two samples.
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html


y_new <- y[keep,]
dim(y_new)

y_new$samples$lib.size <- colSums(y_new$counts)
y_new$samples

y_new <- calcNormFactors(y_new)

pdf("/users/sparthib/retina_lrs/plots/de/edgeR_plotmds_8samples_Jan13.pdf")
plotMDS(y_new, method="bcv", col=as.numeric(y_new$samples$group))
legend("bottomleft", as.character(unique(y_new$samples$group)), col=1:3, pch=20)
dev.off()


####fit GLM estimates of dispersion 

design.mat <- model.matrix(~ 0 + y_new$samples$group)
colnames(design.mat) <- levels(y_new$samples$group)
y2 <- estimateGLMCommonDisp(y_new,design.mat)
y2 <- estimateGLMTrendedDisp(y2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
y2 <- estimateGLMTagwiseDisp(y2,design.mat)

pdf("/users/sparthib/retina_lrs/plots/de/edgeR_plotBCV_8samples_Jan13.pdf")
plotBCV(y2)
dev.off()


et <- exactTest(y2)
topTags(et, n = 10 )

test <- as.data.frame(et)

write.table(test, file='/users/sparthib/retina_lrs/processed_data/de/DE_edgeR_cpm_filter_10_glm_dispersion_8samples_Jan13.tsv', quote=FALSE, sep='\t')


test[test$gene_name == "BSG",]


DIR=/dcs04/hicks/data/sparthib/casey/bams/DG-WT-hRGC
# 
# coolbox add XAxis - \
# add BAMCov $DIR/DG-WT-hRGC_sorted.bam - \
# add Title "bam" - \
# goto "chr19:562360-572228" - \
# plot /users/sparthib/retina_lrs/plots/coverage/tmp/test_coolbox_BSG.jpg
# 
