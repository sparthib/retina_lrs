suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(rtracklayer)
  library(biomaRt)
  library(ggrepel)
  library(biomaRt)
  library(readxl)
  library(cowplot)
  library(here)
  library(paletteer)
})



brainRBPs = read.csv("https://github.com/gandallab/Dev_Brain_IsoSeq/raw/main/data/RBP_Data/CSVs/RBP_targets_v5.csv", header=TRUE);
brainRBPs = dplyr::select(brainRBPs, -c(MGI.symbol, ENSMUSG)) # |>
  #rename("hgnc_symbol"="HGNC.symbol", "ensembl_gene_id"="ENSG")


encodeRBPs = read.csv("https://github.com/gandallab/Dev_Brain_IsoSeq/raw/main/data/RBP_Data/CSVs/RBP_targets_ENCODE.csv", header=TRUE);
encodeRBPs = encodeRBPs |> filter(cell.type=="HepG2") |>
  rename("hgnc_symbol"="HGNC.symbol", "ensembl_gene_id"="ENSG")

rbp_targets = rbind(encodeRBPs, brainRBPs)


mart = useMart("ENSEMBL_MART_ENSEMBL","mmusculus_gene_ensembl")
f = listFilters(mart); a = listAttributes(mart)
featuresToGet = c("ensembl_gene_id", "external_gene_name", 
                  "hsapiens_homolog_ensembl_gene",
                  "hsapiens_homolog_associated_gene_name",
                  "hsapiens_homolog_orthology_type")
mouseHumanHomologs = getBM(attributes = featuresToGet,mart = mart)

human_mouse_bg = mouseHumanHomologs %>% as_tibble() %>% filter(hsapiens_homolog_orthology_type == "ortholog_one2one") %>% dplyr::select("hsapiens_homolog_ensembl_gene") %>% pull()


# tableS3.gene <- read_tsv("https://github.com/gandallab/Dev_Brain_IsoSeq/raw/main/output/tables/TableS3_v3.tsv.gz")

tableS3.gene <- read_tsv("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/DGE_DTU_DTE.tsv")

genesets = list("DTU"= tableS3.gene %>% filter(DTU) %>% mutate(gene_id = substr(gene_id,1,15)) %>% dplyr::select(gene_id) %>% pull(),
                "DTE" = tableS3.gene %>% filter(DTE) %>% mutate(gene_id = substr(gene_id,1,15)) %>% dplyr::select(gene_id) %>% pull(),
                "DGE" = tableS3.gene %>% filter(DGE) %>% mutate(gene_id = substr(gene_id,1,15)) %>% dplyr::select(gene_id) %>% pull(),
                "DTUnotDGE" = tableS3.gene %>% filter(DTU,DGE_pval>.05) %>% mutate(gene_id = substr(gene_id,1,15)) %>%
                  dplyr::select(gene_id) %>% pull())



#### Over-representation analysis functions ####
## Odds-ratio estimator
OR <- function(q,k,m,t) {
  ## 2 x 2 table:
  ##             inTest path    !inTest path 
  ## inRef  path     q            k 
  ## !inRef path     m            t
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

## count overlaps and run the analysis
ORA <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}
##### #####

## background genes, can be all genes in your sample, or all in reference GTF 
DTU.bg = tableS3.gene |> mutate(gene_id = substr(gene_id,1,15)) |> dplyr::select(gene_id) %>% pull() %>% unique()

df_fisher = data.frame()
for(i in 1:length(genesets)) {
  for(this_dataset in unique(na.omit(rbp_targets$dataset.id))) {
    this_rbp = rbp_targets %>% filter(dataset.id == this_dataset) %>% mutate(target = paste0(RBP, "_", data.type, "_", cell.type)) %>% dplyr::select(target) %>% unique()  %>% pull()
    target_genes = rbp_targets %>% filter(dataset.id == this_dataset) %>% dplyr::select(ENSG) %>% pull()
    
    if(grepl("Human",this_rbp)) {
      this_or = ORA(genesets[[i]], target_genes, DTU.bg, DTU.bg)
    } else {
      this_or = ORA(genesets[[i]], target_genes, DTU.bg, human_mouse_bg)
    }
    df_fisher = rbind(df_fisher, data.frame(set = names(genesets)[[i]], dataset = this_dataset, target = this_rbp, t(this_or)))
  }
}

df_fisher$OR = as.numeric(df_fisher$OR)
df_fisher$Fisher.p[df_fisher$OR<1] = 1
df_fisher$Fisher.p = p.adjust(as.numeric(df_fisher$Fisher.p),'fdr')


order.brainRBPs = read_excel(here("/users/sparthib/retina_lrs/raw_data/curatedRBPs_order.xlsx")) |> as_tibble()
order.Encode.TarReg = read_excel(here("/users/sparthib/retina_lrs/raw_data/ENCODE_vanNostrand_NatMeth2016_Fig2a_order.xlsx"), sheet=2);
order.Encode.TarReg = order.Encode.TarReg |> filter(name %in% df_fisher$target)
order = rbind(order.brainRBPs, order.Encode.TarReg)


df_fisher$target = factor(df_fisher$target, levels=order$name)
df_fisher$org = "Mouse"; df_fisher$org[grep("Human",df_fisher$target)] = "Human"
df_fisher$data.type = rbp_targets$data.type[match(df_fisher$dataset, rbp_targets$dataset.id)]
df_fisher$cell.type = rbp_targets$cell.type[match(df_fisher$dataset, rbp_targets$dataset.id)]
df_fisher$target.region = order$target.region[match(df_fisher$target, order$name)]
df_fisher$label = signif(df_fisher$OR,1)
df_fisher$label[df_fisher$Fisher.p>.05] = ''


Fig3H.1 <- ggplot(df_fisher |> filter(set == "DTU"), aes(x = target, y = -log10(Fisher.p), fill = OR)) +
  geom_bar(stat = 'identity', position = position_dodge2()) + 
  theme_bw() +
  geom_hline(yintercept = 1, lty = 'dashed', size = 0.5, color = 'red') + 
  labs(y = 'Enrichment\n(-log10 q-value)', x = '') + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_text(size = 5), 
    axis.title.y = element_text(size = 7), 
    legend.key.size = unit(0.3, 'cm'), 
    legend.text = element_text(size = 4), 
    legend.title = element_text(size = 6), 
    plot.margin = unit(c(0, 0, 0, 0), "pt"), 
    legend.position = c(0.95, 0.75), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line(color = "black", size = 0.2)
  ) + 
  coord_fixed(ratio = 1/3)  # Adjust aspect ratio

Fig3H.2 <- ggplot(df_fisher %>% filter(set == "DTU"), aes(x = target, label = target.region)) + 
  geom_tile(aes(y = factor(1), fill = target.region)) + 
  geom_point(aes(y = factor(1), shape = data.type), position = position_dodge2(width = 1), size = 0.5) + 
  scale_shape_manual(values = c(1:9)) + 
  scale_x_discrete(labels = sapply(strsplit(levels(df_fisher$target), "_"), "[[", 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 3), 
    legend.key.size = unit(0.25, 'cm'), 
    legend.text = element_text(size = 3), 
    legend.title = element_text(size = 5), 
    plot.margin = unit(c(0, 0, 0, 0), "pt"), 
    legend.position = "bottom", 
    legend.box = "horizontal", 
    legend.direction = "horizontal", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line(color = "black", size = 0.2)
  ) + 
  labs(x = '', y = '') + 
  paletteer::scale_fill_paletteer_d("rcartocolor::Vivid") + 
  guides(fill = guide_legend(order = 1)) + 
  coord_fixed(ratio = 1/1.5)
library(patchwork)
pdf("/users/sparthib/retina_lrs/processed_data/dtu/DTU_gandall/bambu/FT_vs_RGC/RBP_plot_1.pdf")
patchwork <- Fig3H.1 / Fig3H.2
print(patchwork)
dev.off()



