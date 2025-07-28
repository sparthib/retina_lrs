library(biomaRt)
library(readr)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(edgeR)

ora_plot <- function(genelist, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = names,
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  readable      = TRUE) 
  if(nrow(as.data.frame(ego)) != 0){
    ego <- enrichplot::pairwise_termsim(ego)
    # ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    write_tsv(as.data.frame(ego), file.path(output_plot_dir,
                                            paste0("simplified_ORA_ASE_DGE_genes_",analysis_type, "BP", ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("simplified_ORA_all_ASE_DGE_genes_",analysis_type, "BP", ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  } else {
    message("No significant GO terms found for ", analysis_type, " in ", "BP")}
}

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/"

go_plot_dir <- file.path(dge_output_dir, "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


samples <- c("H9-BRN3B_hRO_2", "H9-BRN3B-RO", "H9-CRX_hRO_2", "H9-CRX_ROs_D45",
             "H9-FT_1" , "H9-FT_2", "H9-hRGC_1", "H9-hRGC_2",
             "EP1-WT_ROs_D45", "EP1-BRN3B-RO", "EP1-WT_hRO_2") 

gene_counts_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/H9_EP1_gene_counts_all_samples"

# load all counts matrix in the directory
files <- list.files(gene_counts_dir,
                    pattern = "_counts.txt", 
                    full.names = TRUE)

# create a list of data frames

# Load all the count data into a list of data frames
counts_list <- lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1)  # Adjust if needed (e.g., use sep = "," for CSV)
})

# Combine all data frames into a single matrix by column binding them
counts_matrix <- do.call(cbind, counts_list)

## remove rows with 0 values

counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]

colnames(counts_matrix) <- gsub(".bam", "", colnames(counts_matrix))


#remove version number from row names
rownames(counts_matrix) <- gsub("\\..*", "", rownames(counts_matrix))

# only keep PTC genes

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(
  mart=mart,
  attributes=c( 
    "external_gene_name",
    "gene_biotype", "ensembl_gene_id",  "chromosome_name"),
  filter="ensembl_gene_id",
  values=rownames(counts_matrix),
  uniqueRows=TRUE)

head(annotLookup)
colnames(annotLookup) <- c("gene_name","gene_biotype", "gene_id", "chromosome_name" )
annotLookup <- annotLookup |> dplyr::distinct()
annotLookup |> nrow()
annotLookup <- annotLookup |> 
  dplyr::filter(gene_biotype == "protein_coding") 
annotLookup |> nrow()

# filter for protein coding genes
counts_matrix <- counts_matrix[rownames(counts_matrix) 
                        %in% annotLookup$gene_id, ]
nrow(counts_matrix)

culture_type <- c("EP1", "EP1", "EP1", "EP1", "EP1", "EP1",
                  "H9", "H9", "H9", "H9", "H9", "H9",
                  "H9", "H9", "H9", "H9", "H9", "H9",
                  "H9", "H9", "H9", "H9")

alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")

groups <- c("Stage3", "Stage3", 
            "Stage2", "Stage2",
            "Stage1", "Stage1",
            "Stage2", "Stage2",
            "Stage3", "Stage3",
            "Stage2", "Stage2",
            "Stage1", "Stage1",
            "FT", "FT",  
            "FT", "FT",
            "RGC", "RGC",
            "RGC", "RGC")

column_data <- data.frame(
  row.names = colnames(counts_matrix),
  culture_type = culture_type,
  allele = alleles,
  group = groups
)

design <- data.frame(
  sample = colnames(counts_matrix),
  group = factor(groups, levels = c("Stage1", "Stage2", "Stage3", "FT", "RGC")),
  allele = factor(alleles, levels = c("H1", "H2")),
  culture_type  = factor(culture_type , levels = c("H9", "EP1"))
)

design <- model.matrix(~ 0 + group*allele + allele*culture_type + group*allele,
                       data = design)
# > colnames(design)
# [1] "groupStage1"              "groupStage2"             
# [3] "groupStage3"              "groupFT"                 
# [5] "groupRGC"                 "alleleH2"                
# [7] "culture_typeEP1"          "groupStage2:alleleH2"    
# [9] "groupStage3:alleleH2"     "groupFT:alleleH2"        
# [11] "groupRGC:alleleH2"        "alleleH2:culture_typeEP1"

### H1 vs. H2 is coef 6
### H9 vs EP1 is coef 7

y <- DGEList(counts = counts_matrix,
             samples = colnames(counts_matrix),
             group = groups,
             allele = alleles,
             genes = rownames(counts_matrix))

y <- normLibSizes(y)


y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)

qlf_allele <- glmQLFTest(fit, coef=6)
qlf_cell_line <- glmQLFTest(fit, coef = 7)

is.de_allele <- decideTests(qlf_allele, p.value=0.05)
is.de_cell_line <- decideTests(qlf_cell_line, p.value=0.05)

# summary(is.de_allele)
# summary(is.de_cell_line)
# > summary(is.de_allele)
# alleleH2
# Down        255
# NotSig    17087
# Up          175
# > summary(is.de_cell_line)
# culture_typeEP1
# Down              2952
# NotSig           12273
# Up                2292


# Extract the results
tt_allele <- topTags(qlf_allele,n = Inf)
nrow(tt_allele) 
head(tt_allele)
colnames(tt_allele$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")


tt_allele$condition_1 <- "H1"
tt_allele$condition_2 <- "H2"

tt_cell_line <- topTags(qlf_cell_line, n = Inf)
nrow(tt_cell_line)
colnames(tt_cell_line$table) <- c("gene_id", "logFC", "logCPM", "F" , "PValue", "FDR")
tt_cell_line$condition_1 <- "H9"
tt_cell_line$condition_2 <- "EP1"


merge_table <- function(tt, annotLookup) {
  # Merge the annotation with the results
  merged_results <- merge(tt$table, annotLookup,
                          by = "gene_id", all.x = TRUE)
  
  # Filter significant results
  merged_results_sig <- merged_results |> 
    dplyr::filter(FDR < 0.05 & abs(logFC) >= 1)
  
  # Add -log10 FDR for plotting
  merged_results <- merged_results |> 
    mutate(neg_log10_FDR = -log10(FDR))
  
  # Separate out the first 20 genes for labeling
  merged_results$label <- NA
  merged_results <- merged_results |>
    arrange(desc(abs(logFC)), FDR)
  merged_results$label[1:20] <- merged_results$gene_name[1:20]
  
  merged_results$significant <- merged_results$FDR < 0.05 & 
    abs(merged_results$logFC) > 1
  
  return(merged_results)
}
  
merged_results_allele <- merge_table(tt_allele, annotLookup)
merged_results_cell_line <- merge_table(tt_cell_line, annotLookup)

write_tsv(merged_results_allele,
          file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_DE_results.tsv"))
write_tsv(merged_results_cell_line,
          file = file.path(dge_output_dir, "ASE_H9_vs_EP1_PTC_DE_results.tsv"))

write_tsv(counts_matrix,
          file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_gene_counts.tsv"))

# Volcano plot
pdf(file = file.path(dge_output_dir, "ASE_H1_vs_H2_PTC_DE_volcano_plot.pdf"),
    width = 8, height = 6)

ggplot(merged_results_allele, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "H1 vs H2, all samples",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()
dev.off()

pdf(file = file.path(dge_output_dir, "ASE_H9_vs_EP1_PTC_DE_volcano_plot.pdf"),
    width = 8, height = 6)
ggplot(merged_results_cell_line, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(aes(color = significant), size = 1.5) +  # TRUE/FALSE coloring
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "H9 vs EP1, all samples",
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme_minimal()
dev.off()



# GO plot 
tt <- merged_results_allele
tt <- tt |> filter(significant == TRUE)
names <- tt |> pull(gene_id)
values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
names(values) <- names
head(values)

ora_plot(genelist = names, 
         output_plot_dir = go_plot_dir,
         analysis_type = "H1_vs_H2")

tt <- merged_results_cell_line
tt <- tt |> filter(significant == TRUE)
names <- tt |> pull(gene_id)
values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
names(values) <- names
head(values)

ora_plot(genelist = names, 
         output_plot_dir = go_plot_dir,
         analysis_type = "H9_vs_EP1")


##### H9 vs EP1 and H1 vs H2 ####
dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/DGE/H9_vs_EP1/"
dir.create(dge_output_dir, showWarnings = FALSE)


contr <- makeContrasts( H9_H1_vs_H2 = H9_H2 - H9_H1, 
                        EP1_H1_vs_H2 = EP1_H2 - EP1_H1, 
                        H2_H9_vs_EP1 =  EP1_H2- H9_H2 ,
                        H1_H9_vs_EP1 = EP1_H1- H9_H1 ,
                        levels=design)

for (i in seq_len(ncol(contr))) {
  qlf <- glmQLFTest(fit, contrast = contr[,i])
  is.de <- decideTests(qlf, p.value=0.05)
  
  tt <- topTags(qlf, n = Inf)$table
  tt$gene_id <- gsub("\\..*", "", tt$genes)
  
  annotLookup <- getBM(mart=mart, 
                       attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"), 
                       filter="ensembl_gene_id", 
                       values=tt$gene_id, uniqueRows=TRUE)
  colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype", "chromosome_name")
  
  tt <- merge(tt, annotLookup, by="gene_id", all.x=TRUE)
  tt <- tt[order(tt$FDR), ]
  
  tt$condition_1 <- names(contr[,i])[contr[,i] == -1]
  tt$condition_2 <- names(contr[,i])[contr[,i] == 1]
  tt$neg_log10_FDR <- -log10(tt$FDR)
  tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
  tt <- tt |>
    arrange(desc(abs(logFC)), FDR)
  tt$label <- NA
  tt$label[1:20] <- tt$gene_name[1:20]
  file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  write_tsv(tt, file)
}

                       
for (i in seq_len(ncol(contr))) {
  read_file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  tt <- read_tsv(read_file)
  pdf(file = paste0(dge_output_dir, colnames(contr)[i], "_volcano_plot.pdf"),
      width = 8, height = 6)
  
  print(
    ggplot(tt, aes(x = logFC, y = neg_log10_FDR)) +
      geom_point(aes(color = significant), size = 1.5) +
      geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      scale_color_manual(values = c("gray", "red")) +
      labs(title = colnames(contr)[i],
           x = "log2 Fold Change",
           y = "-log10 FDR") +
      theme_minimal()
  )
  
  dev.off()
}





for (i in c(1,2,3,4)) {
  read_file <- paste0(dge_output_dir, colnames(contr)[i], "_DGEs.tsv")
  tt <- read_tsv(read_file)
  tt <- tt |> filter(significant == TRUE)
  # check for NAs in gene_id
  if (any(is.na(tt$gene_id))) {
    message("NAs found in gene_id column of the DGE results.")
  }
  
  # GO plot 
  names <- tt |> pull(gene_id)
  values <- tt |> pull(logFC) |> sort(decreasing = TRUE)
  names(values) <- names
  head(values)
  
  ora_plot(genelist = names, 
           output_plot_dir = go_plot_dir,
           analysis_type = colnames(contr)[i])
}







