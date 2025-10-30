library(dplyr)
library(biomaRt)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

source("/users/sparthib/retina_lrs/code/08_ASE/short_reads/02_feature_counts_analysis/volcano_helper.R")


gene_counts <- read_tsv("/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_counts.tsv")

samples <- colnames(gene_counts)

ncol(gene_counts)
nrow(gene_counts)

rownames(gene_counts) <- gene_counts$gene_id
counts_matrix <- as.matrix(gene_counts[,
                                       -which(names(gene_counts) == "gene_id")])

rownames(counts_matrix) <- gene_counts$gene_id
counts_matrix <- counts_matrix[,1:14]

alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")

stages <- c("Stage_3", "Stage_3", 
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "Stage_2", "Stage_2",
            "Stage_3", "Stage_3",
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1")


cell_lines <-  c("EP1", "EP1", "EP1", "EP1", "EP1", "EP1",
                 "H9", "H9", "H9", "H9", "H9", "H9",
                 "H9", "H9")

groups <- paste0(stages,"_", alleles)

samples <- colnames(counts_matrix)


# Create a proper samples data frame
samples_df <- data.frame(
  sample = colnames(counts_matrix),
  cell_line = factor(cell_lines, levels = c("H9", "EP1")),
  group = factor(groups),
  allele = factor(alleles, levels = c("H1", "H2")),
  stage = factor(stages, levels = c("Stage_1", "Stage_2", "Stage_3"))
)

# Create DGEList with proper structure
y <- DGEList(
  counts = counts_matrix,
  samples = samples_df
)

# Now create the design matrix
design <- model.matrix(~ allele*cell_line, 
                       data = y$samples)

y <- normLibSizes(y)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

dge_output_dir <- "/users/sparthib/retina_lrs/processed_data/ASE/bambu_counts_matrices/DGE/cell_line_vs_allele"
dir.create(dge_output_dir, showWarnings = FALSE)

# Fix the column names to make them syntactically valid
colnames(design) <- make.names(colnames(design))

# Check the new names
colnames(design)

contrasts_H1vH2 <- makeContrasts(
  H1_H9_vs_H2_H9 = alleleH2,
  H1_EP1_vs_H2_EP1 = alleleH2 + alleleH2.cell_lineEP1,
  levels = design
)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

run_dge_contrasts <- function(fit, contrasts, mart, dge_output_dir) {
  # Ensure output directory exists
  if (!dir.exists(dge_output_dir)) {
    dir.create(dge_output_dir, recursive = TRUE)
  }
  
  for (i in seq_len(ncol(contrasts))) {
    contrast_name <- colnames(contrasts)[i]
    message("Running contrast: ", contrast_name)
    
    
    # Run QL F-test
    qlf <- glmQLFTest(fit, contrast = contrasts[, i])
    
    # Get full result table
    tt <- topTags(qlf, n = Inf)$table
    colnames(tt) <- c("logFC", "logCPM", "F", "PValue", "FDR")
    tt$gene_id <- rownames(tt)
    
    # Annotate using biomaRt
    annotLookup <- tryCatch({
      getBM(
        mart = mart,
        attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"),
        filter = "ensembl_gene_id",
        values = tt$gene_id,
        uniqueRows = TRUE
      )
    }, error = function(e) {
      warning("Annotation failed for ", contrast_name, ": ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(annotLookup)) {
      colnames(annotLookup) <- c("gene_id", "gene_name", "gene_biotype", "chromosome_name")
      tt <- merge(tt, annotLookup, by = "gene_id", all.x = TRUE)
    } else {
      tt$gene_name <- NA
      tt$gene_biotype <- NA
      tt$chromosome_name <- NA
    }
    
    # Order and compute extra stats
    tt <- tt[order(tt$FDR), ]
    tt$neg_log10_FDR <- -log10(tt$FDR)
    tt$significant <- tt$FDR < 0.05 & abs(tt$logFC) > 1
    
    # Parse condition names if formatted as X_vs_Y
    tt$condition_1 <- gsub("_vs_.*", "", contrast_name)
    tt$condition_2 <- gsub(".*_vs_", "", contrast_name)
    
    # Add top labels
    tt <- tt |> dplyr::arrange(desc(abs(logFC)), FDR)
    tt$label <- NA
    tt$label[1:min(20, nrow(tt))] <- tt$gene_name[1:min(20, nrow(tt))]
    
    # Write output file
    outfile <- file.path(dge_output_dir, paste0(contrast_name, "_DGEs.tsv"))
    readr::write_tsv(tt, outfile)
  }
  
  message("✅ All contrasts processed and results saved in: ", dge_output_dir)
}


run_dge_contrasts(fit, contrasts_H1vH2, mart, dge_output_dir)


contrast_names_H1vH2 <- c("H1_H9_vs_H2_H9", "H1_EP1_vs_H2_EP1")

# Generate volcano plots
generate_allele_volcano_plots(
  input_dir = dge_output_dir,
  contrast_names = contrast_names_H1vH2,
  table_type = "DGE"
)

# Run GO plots for H1 vs H2 contrasts
H1_v_H2_go_plot_dir <- file.path( dge_output_dir, "GO_plots")
dir.create(H1_v_H2_go_plot_dir, showWarnings = FALSE)
run_go_plots(contrasts_H1vH2, dge_output_dir,H1_v_H2_go_plot_dir)


