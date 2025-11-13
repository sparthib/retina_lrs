library(biomaRt)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)

#helper function for volcano plots
code_dir <- Sys.getenv("retina_lrs_code")
data_dir <- Sys.getenv("retina_lrs_dir")
ref_dir <- Sys.getenv("references_dir")

source(file.path(code_dir,"code/08_ASE/short_reads/02_bambu_analysis/volcano_helper.R"))


gene_counts <- read_tsv(file.path(code_dir,
                                  "processed_data/ASE/bambu_counts_matrices/bambu_ptc_gene_counts.tsv"))

samples <- colnames(gene_counts)

ncol(gene_counts)
nrow(gene_counts)

rownames(gene_counts) <- gene_counts$gene_id
counts_matrix <- as.matrix(gene_counts[,
                                       -which(names(gene_counts) == "gene_id")])

rownames(counts_matrix) <- gene_counts$gene_id

alleles <- c("H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2", "H1", "H2",
             "H1", "H2")

stages <- c("Stage_3", "Stage_3", 
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "Stage_2", "Stage_2",
            "Stage_3", "Stage_3",
            "Stage_2", "Stage_2",
            "Stage_1", "Stage_1",
            "FT", "FT",  
            "FT", "FT",
            "RGC", "RGC",
            "RGC", "RGC")


cell_lines <-  c("EP1", "EP1", "EP1", "EP1", "EP1", "EP1",
                 "H9", "H9", "H9", "H9", "H9", "H9",
                 "H9", "H9", "H9", "H9", "H9", "H9",
                 "H9", "H9", "H9", "H9")

groups <- paste0(stages,"_", alleles)

samples <- colnames(counts_matrix)


# Create a proper samples data frame
samples_df <- data.frame(
  sample = colnames(counts_matrix),
  cell_line = factor(cell_lines, levels = c("H9", "EP1")),
  group = factor(groups),
  allele = factor(alleles, levels = c("H1", "H2")),
  stage = factor(stages, levels = c("Stage_1", "Stage_2", "Stage_3", "FT", "RGC"))
)

# Create DGEList with proper structure
y <- DGEList(
  counts = counts_matrix,
  samples = samples_df
)

# Now create the design matrix
design <- model.matrix(~ allele*stage, data = y$samples)

y <- normLibSizes(y)

y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)


dge_output_dir <- file.path(code_dir,"processed_data/ASE/bambu_counts_matrices/DGE")
dir.create(dge_output_dir, showWarnings = FALSE)


# Fix the column names to make them syntactically valid
colnames(design) <- make.names(colnames(design))

# Check the new names
colnames(design)

contrasts_H1vH2 <- makeContrasts(
  H1_Stage1_vs_H2_Stage1 = alleleH2,
  H1_Stage2_vs_H2_Stage2 = alleleH2 + alleleH2.stageStage_2,
  H1_Stage3_vs_H2_Stage3 = alleleH2 + alleleH2.stageStage_3,
  H1_FT_vs_H2_FT   = alleleH2 + alleleH2.stageFT,
  H1_RGC_vs_H2_RGC    = alleleH2 + alleleH2.stageRGC,
  levels = design
)

contrasts_H1 <- makeContrasts(
  H1_Stage1_vs_H1_Stage2 = stageStage_2,
  H1_Stage2_vs_H1_Stage3 = stageStage_3 - stageStage_2,
  H1_Stage1_vs_H1_Stage3 = stageStage_3,
  H1_FT_vs_H1_RGC         = stageRGC - stageFT,
  H1_Stage1_vs_H1_RGC     = stageRGC,      
  H1_Stage2_vs_H1_RGC     = stageRGC - stageStage_2,
  H1_Stage3_vs_H1_RGC     = stageRGC - stageStage_3,
  levels = design
)

contrasts_H2 <- makeContrasts(
  H2_Stage1_vs_H2_Stage2 = (stageStage_2 + alleleH2.stageStage_2),
  H2_Stage2_vs_H2_Stage3 = (stageStage_3 + alleleH2.stageStage_3) - (stageStage_2 + alleleH2.stageStage_2),
  H2_Stage1_vs_H2_Stage3 = (stageStage_3 + alleleH2.stageStage_3),
  H2_FT_vs_H2_RGC         = (stageRGC + alleleH2.stageRGC) - (stageFT + alleleH2.stageFT),
  H2_Stage1_vs_H2_RGC     = (stageRGC + alleleH2.stageRGC),
  H2_Stage2_vs_H2_RGC     = (stageRGC + alleleH2.stageRGC) - (stageStage_2 + alleleH2.stageStage_2),
  H2_Stage3_vs_H2_RGC     = (stageRGC + alleleH2.stageRGC) - (stageStage_3 + alleleH2.stageStage_3),
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



# Run DGE contrasts for H1 vs H2 at each stage
H1_v_H2_output_dir <- file.path(dge_output_dir, "H1_vs_H2_contrasts")
run_dge_contrasts(fit, contrasts_H1vH2, mart, H1_v_H2_output_dir)


contrast_names_H1vH2 <- c("H1_Stage1_vs_H2_Stage1", "H1_Stage2_vs_H2_Stage2", "H1_Stage3_vs_H2_Stage3",
 "H1_FT_vs_H2_FT" ,"H1_RGC_vs_H2_RGC")

# Generate volcano plots
generate_allele_volcano_plots(
  input_dir = H1_v_H2_output_dir,
  contrast_names = contrast_names_H1vH2,
  table_type = "DGE"
)

#Run DGE for H1 across stages
H1_output_dir <- file.path(dge_output_dir, "between_stage_in_H1_allele")
dir.create(H1_output_dir, showWarnings = FALSE)
run_dge_contrasts(fit, contrasts_H1, mart, H1_output_dir)
contrast_names_H1 <- c("H1_Stage1_vs_H1_Stage2", "H1_Stage2_vs_H1_Stage3", "H1_Stage1_vs_H1_Stage3",
                        "H1_FT_vs_H1_RGC", "H1_Stage1_vs_H1_RGC", "H1_Stage2_vs_H1_RGC", "H1_Stage3_vs_H1_RGC")
# Generate volcano plots
generate_allele_volcano_plots(
  input_dir = H1_output_dir,
  contrast_names = contrast_names_H1,
  table_type = "DGE"
)

#Run DGE for H2 across stages
H2_output_dir <- file.path(dge_output_dir, "between_stage_in_H2_allele")
dir.create(H2_output_dir, showWarnings = FALSE)
run_dge_contrasts(fit, contrasts_H2, mart, H2_output_dir)
contrast_names_H2 <- c("H2_Stage1_vs_H2_Stage2", "H2_Stage2_vs_H2_Stage3", "H2_Stage1_vs_H2_Stage3",
                        "H2_FT_vs_H2_RGC", "H2_Stage1_vs_H2_RGC", "H2_Stage2_vs_H2_RGC", "H2_Stage3_vs_H2_RGC")

# Generate volcano plots
generate_allele_volcano_plots(
  input_dir = H2_output_dir,
  contrast_names = contrast_names_H2,
  table_type = "DGE"
)

ora_plot <- function(genelist, output_plot_dir, analysis_type){
  
  ego <- enrichGO(gene          = genelist,
                  OrgDb         = org.Hs.eg.db,
                  keyType  = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  readable      = TRUE) 
  if(nrow(as.data.frame(ego)) != 0){
    #ego <- enrichplot::pairwise_termsim(ego)
    # ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    write_tsv(as.data.frame(ego), file.path(output_plot_dir,
                                            paste0("simplified_ORA_ASE_DGE_genes_",analysis_type, "BP", ".tsv")))
    
    
    pdf(file.path(output_plot_dir, paste0("simplified_ORA_all_ASE_DGE_genes_",analysis_type, "BP", ".pdf")))
    print(dotplot(ego, showCategory = 15))
    dev.off()
    
  } else {
    message("No significant GO terms found for ", analysis_type, " in ", "BP")}
}


go_plot_dir <- file.path(dge_output_dir, "plots", "GO_plots")
if (!dir.exists(go_plot_dir)) {
  dir.create(go_plot_dir, recursive = TRUE)
}


run_go_plots <- function(contrasts, dge_output_dir, go_plot_dir) {
  for (i in seq_len(ncol(contrasts))) {
    contrast_name <- colnames(contrasts)[i]
    read_file <- file.path(dge_output_dir, paste0(contrast_name, "_DGEs.tsv"))
    
    if (!file.exists(read_file)) {
      message("⚠️ File not found: ", read_file)
      next
    }
    
    tt <- readr::read_tsv(read_file, show_col_types = FALSE)
    
    if (!("significant" %in% colnames(tt))) {
      message("⚠️ No 'significant' column in: ", contrast_name)
      next
    }
    
    tt <- tt |> dplyr::filter(significant == TRUE)
    
    if (nrow(tt) == 0) {
      message("ℹ️ No significant genes for: ", contrast_name)
      next
    }
    
    # Prepare ranked gene list
    values <- tt |> dplyr::pull(logFC)
    names(values) <- tt |> dplyr::pull(gene_id)
    values <- sort(values, decreasing = TRUE)
    
    # Run GO enrichment plot
    ora_plot(
      genelist = as.vector(names(values)),
      output_plot_dir = go_plot_dir,
      analysis_type = contrast_name
    )
    
    message("✅ GO plot done for: ", contrast_name)
  }
}

# Run GO plots for H1 vs H2 contrasts
H1_v_H2_go_plot_dir <- file.path( H1_v_H2_output_dir, "GO_plots")
dir.create(H1_v_H2_go_plot_dir, showWarnings = FALSE)
run_go_plots(contrasts_H1vH2, H1_v_H2_output_dir,H1_v_H2_go_plot_dir)

# Run GO plots for H1 across stages
H1_go_plot_dir <- file.path(H1_output_dir, "GO_plots")
dir.create(H1_go_plot_dir, showWarnings = FALSE)
run_go_plots(contrasts_H1, H1_output_dir,H1_go_plot_dir)

# Run GO plots for H2 across stages
H2_go_plot_dir <- file.path(H2_output_dir, "GO_plots")
dir.create(H2_go_plot_dir, showWarnings = FALSE)
run_go_plots(contrasts_H2, H2_output_dir,H2_go_plot_dir)

