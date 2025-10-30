library(readr)
library(dplyr)
library(biomaRt)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


library(readr)
library(dplyr)
library(biomaRt)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


filtered_cpm_matrix <- read_tsv("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_cpm.tsv")
filtered_counts_matrix <- read_tsv("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_counts.tsv")

EP1_gene_cpm <- filtered_cpm_matrix[, grepl("EP1", colnames(filtered_cpm_matrix))]
EP1_gene_counts <- filtered_counts_matrix[, grepl("EP1", colnames(filtered_counts_matrix))]
# perform DGE 
alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2")
design <- data.frame(
  sample = colnames(EP1_gene_counts),
  allele = factor(alleles, levels = c("H1", "H2"))
)
design <- model.matrix(~allele , data = design)
#coef 4 (alleleH2)

y <- DGEList(counts = EP1_gene_counts,
             samples = colnames(EP1_gene_counts),
             allele = alleles,
             genes = rownames(EP1_gene_counts))

keep <- filterByExpr(y, design = design, min.count = 3)  # or other threshold
y <- y[keep , keep.lib.sizes=FALSE]
y <- normLibSizes(y)

y <- estimateDisp(y, design, robust=TRUE)

fit <- glmQLFit(y, design, robust=TRUE)
design
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
colnames(design)
colnames(fit$coefficients)

results <- topTags(qlf, n = Inf)$table
head(results)

sig_results <- results |> #genes highly expressed in H2 allele
  filter(FDR < 0.05 & abs(logFC) >= 1)

nrow(results)
nrow(sig_results)

saveRDS(sig_results, file = "/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/EP1_ASE_DGE_H2_sig_results.rds")
sig_results <- readRDS("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/EP1_ASE_DGE_H2_sig_results.rds")


dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices"

go_plot_dir <- file.path(dge_output_dir, "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


ora_plot <- function(genelist, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = genelist,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  readable      = TRUE) 
  
  ego_df <- as.data.frame(ego)
  
  if (nrow(ego_df) > 0) {
    write_tsv(
      ego_df,
      file.path(output_plot_dir,
                paste0("simplified_ORA_ASE_DGE_genes_", analysis_type, "_BP.tsv"))
    )
    
    pdf(file.path(output_plot_dir,
                  paste0("simplified_ORA_ASE_DGE_genes_", analysis_type, "_BP.pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  } else {
    message("No significant GO terms found for ", analysis_type, " in BP")
  }
}


names <- sig_results |> pull(genes) |> unique()
ora_plot(names, go_plot_dir, "EP1_ASE_DGE_H2_sig")

neg_results <- sig_results |> 
  filter(FDR < 0.05 & logFC <= -1)
names_neg <- neg_results |> pull(genes) |> unique()
ora_plot(names_neg, go_plot_dir, "EP1_ASE_DGE_H1_downregulated_in_H2")

pos_results <- sig_results |> 
  filter(FDR < 0.05 & logFC >= 1)
names_pos <- pos_results |> pull(genes) |> unique()
ora_plot(names_pos, go_plot_dir, "EP1_ASE_DGE_H2_upregulated_in_H2")

write_tsv(sig_results, 
          file = file.path(dge_output_dir, "EP1_ASE_DGE_H2_all_results.tsv"))
