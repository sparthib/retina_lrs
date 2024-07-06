library(clusterProfiler)
library(readr)
library(dplyr)
library(org.Hs.eg.db)

sample_compare <- "Bambu FT_vs_RGC"
output_plot_dir <- "/users/sparthib/retina_lrs/plots/go_analysis/bambu/FT_vs_RGC/"
input_data_dir <- "/users/sparthib/retina_lrs/processed_data/dge/edgeR/bambu/FT_vs_RGC/"
output_data_dir <- "/users/sparthib/retina_lrs/processed_data/go_analysis/bambu/FT_vs_RGC/"
dge <- read_tsv(paste0(input_data_dir, "DGEs.tsv"),
                 col_names = TRUE)
colnames(dge)
# [1] "gene_id"      "genes"        "logFC"        "logCPM"       "F"
# [6] "PValue"       "FDR"          "gene_name"    "gene_biotype"


dge <- dge |> filter(FDR < 0.05) 
gene_list <- dge$logFC
names(gene_list) <- dge$gene_id
gene_list=sort(gene_list, decreasing = TRUE)

ggo_cc <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               keyType  = "ENSEMBL",
               level    = 3,
               readable = TRUE)

write.table(ggo_cc, file = paste0(output_data_dir,
                               "go_cc_classification.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

ggo_bp <- groupGO(gene     = gene,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  keyType  = "ENSEMBL",
                  level    = 3,
                  readable = TRUE)

write.table(ggo_bp, file = paste0(output_data_dir,
                                  "go_bp_classification.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


ggo_mf <- groupGO(gene     = gene,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "MF",
                  keyType  = "ENSEMBL",
                  level    = 3,
                  readable = TRUE)

write.table(ggo_mf, file = paste0(output_data_dir,
                                  "go_mf_classification.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


##### Gene Set Enrichment Analysis #####
#code source fhttps://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

ego <- gseGO(geneList     = gene_list,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              keyType      = "ENSEMBL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
             eps = 0
             )

write.table(ego, file = paste0(output_data_dir,
                               "go_gene_set_enrichment.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)


#dotplot 
require(DOSE)
output_plot_dir <- "/users/sparthib/retina_lrs/plots/go_analysis/bambu/FT_vs_RGC/"
pdf(paste0(output_plot_dir, "dotplot_go_analysis.pdf"))
p <- dotplot(ego, showCategory = 10, title = paste0(sample_compare, "GO Enrichment Analysis"), 
        split=".sign") + facet_grid(.~.sign)
print(p)
dev.off()

#emapplot 
x2 <- enrichplot::pairwise_termsim(ego)
pdf(paste0(output_plot_dir, "enrichment_map_go_analysis.pdf"))
p <- emapplot(x2,showCategory = 10, title=paste0(sample_compare, "GO Enrichment Map"))
print(p)
dev.off()

#cnet plot
pdf(paste0(output_plot_dir, "cnet_go_analysis.pdf"))
p <- cnetplot(ego, categorySize="pvalue", foldChange=gene_list, 
              title=paste0(sample_compare, "CNET plot"), showCategory = 3)
print(p)
dev.off()

#ridgeplot
pdf(paste0(output_plot_dir, "ridgeplot_go_analysis.pdf"))
p <- ridgeplot(ego) + ggplot2::labs(title = paste0(sample_compare, "Enrichment Distribution"))
print(p)
dev.off()

#### Over Representation Analysis ####
ego2_BP <- enrichGO(gene          = names(gene_list),
                # universe      = names(gene_list),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
write.table(ego2_BP, file = paste0(output_data_dir,
                               "go_gene_over_representation_bp.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

ego2_CC <- enrichGO(gene          = names(gene_list),
                    # universe      = names(gene_list),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = "CC",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
write.table(ego2_CC, file = paste0(output_data_dir,
                                   "go_gene_over_representation_cc.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

ego2_MF <- enrichGO(gene          = names(gene_list),
                    # universe      = names(gene_list),
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",
                    ont           = "MF",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
write.table(ego2_MF, file = paste0(output_data_dir,
                                   "go_gene_over_representation_mf.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

#goplot
pdf(paste0(output_plot_dir, "goplot_bp.pdf"))
p <- goplot(ego2_BP, showCategory = 10, title = paste0(sample_compare, "GO Overrepresentation Analysis BP"))
print(p)
dev.off()

pdf(paste0(output_plot_dir, "goplot_cc.pdf"))
p <- goplot(ego2_CC, showCategory = 10, title = paste0(sample_compare, "GO Overrepresentation Analysis CC"))
print(p)
dev.off()

pdf(paste0(output_plot_dir, "goplot_mf.pdf"))
p <- goplot(ego2_MF, showCategory = 10, title = paste0(sample_compare, "GO Overrepresentation Analysis MF"))
print(p)
dev.off()
