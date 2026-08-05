#!/usr/bin/env Rscript
# 10_filter_genes_isoforms_by_biotype.R
# Short-read (salmon) analog of 04_dtu_dge_dte/01c_filter_by_gene_biotypes.R.
# Produces protein-coding + expression-filtered gene and isoform count matrices
# (and TMM-CPM) from the raw salmon/tximport matrices. NO differential testing.
#
# Mirrors 01c exactly:
#   * biotype key = gene_biotype: keep genes with gene_biotype=="protein_coding";
#     keep ALL isoforms belonging to those protein-coding genes.
#   * gene matrix : filterByExpr(group, min.count=10) + TMM
#   * isoform matrix: filterByExpr(group, min.count=2)  + TMM
# Differences from LR: salmon has known ENST only (no novel/Bambu isoforms), so
# there is no common-isoform (01b) step. gene_biotype is read from the local
# GENCODE release-46 primary-assembly GTF (offline; = same annotation salmon used).
# group = the 7 timepoints (filterByExpr grouping factor only, not a DE contrast).
#
# Run with: module load conda_R/4.5.x && Rscript 10_filter_genes_isoforms_by_biotype.R

.libPaths("/users/sparthib/R/4.5.x")
suppressMessages(library(edgeR))
strip <- function(x) sub("\\.[0-9]+$", "", x)

dd  <- "/dcs04/hicks/data/sparthib/retina_lrs"
ref <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
sal <- file.path(dd, "10_short_reads/salmon")
gtf <- file.path(ref, "release_46_primary_assembly.gtf")
out_dir <- file.path(sal, "filtered_by_counts_and_biotype")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- gene_biotype from GTF gene lines ----
gl <- system(paste0("awk -F'\\t' '$3==\"gene\"' ", gtf,
                    " | sed -E 's/.*gene_id \"([^\"]+)\".*gene_type \"([^\"]+)\".*/\\1\\t\\2/'"),
             intern = TRUE)
gm <- do.call(rbind, strsplit(gl, "\t"))
gene_biotype <- setNames(gm[,2], strip(gm[,1]))
pc_genes <- names(gene_biotype)[gene_biotype == "protein_coding"]
cat("GTF genes:", length(gene_biotype), " protein_coding:", length(pc_genes), "\n")

## ---- isoform -> gene map (salmon tx2gene) ----
t2g <- read.table(file.path(sal, "tx2gene_release46.tsv"), header = FALSE,
                  sep = "\t", stringsAsFactors = FALSE)
iso2gene <- setNames(strip(t2g$V2), strip(t2g$V1))

## ---- raw matrices ----
gmat <- as.matrix(read.csv(file.path(sal, "gene_counts_matrix.csv"), row.names = 1, check.names = FALSE))
imat <- as.matrix(read.csv(file.path(sal, "transcript_counts_matrix.csv"), row.names = 1, check.names = FALSE))
rownames(gmat) <- strip(rownames(gmat)); rownames(imat) <- strip(rownames(imat))
group <- factor(sub("_.*", "", colnames(gmat)))            # 7 timepoints
cat("group levels:", paste(levels(group), collapse=","), " | sizes:",
    paste(table(group), collapse=","), "\n")
stopifnot(identical(colnames(gmat), colnames(imat)))

## ---- filter helper: biotype subset -> filterByExpr -> TMM ----
filter_counts <- function(mat, keep_ids, min.count, group) {
  m <- mat[rownames(mat) %in% keep_ids, , drop = FALSE]
  pre <- nrow(m)
  dge <- DGEList(counts = m)
  keep <- filterByExpr(dge, group = group, min.count = min.count)
  dge <- dge[keep, ]
  dge <- calcNormFactors(dge)                                # TMM
  list(counts = m[keep, , drop = FALSE], cpm = cpm(dge, normalized.lib.sizes = TRUE),
       pre = pre, post = sum(keep))
}

## ---- genes: keep protein-coding genes, min.count=10 ----
g <- filter_counts(gmat, pc_genes, 10, group)
cat("GENES  : biotype-subset", g$pre, "-> expr-filtered", g$post, "\n")
saveRDS(g$counts, file.path(out_dir, "filtered_gene_counts.RDS"))
saveRDS(g$cpm,    file.path(out_dir, "filtered_gene_counts_cpm.RDS"))

## ---- isoforms: keep isoforms of protein-coding genes, min.count=2 ----
pc_isoforms <- names(iso2gene)[iso2gene %in% pc_genes]
i <- filter_counts(imat, pc_isoforms, 2, group)
cat("ISOFORM: biotype-subset", i$pre, "-> expr-filtered", i$post, "\n")
saveRDS(i$counts, file.path(out_dir, "filtered_isoform_counts.RDS"))
saveRDS(i$cpm,    file.path(out_dir, "filtered_isoform_cpm.RDS"))

## ---- report: genes represented in each ----
g_in_iso <- length(unique(iso2gene[rownames(i$counts)]))
summ <- data.frame(
  matrix = c("gene", "isoform"),
  biotype_subset = c(g$pre, i$pre),
  expr_filtered  = c(g$post, i$post),
  unique_genes   = c(nrow(g$counts), g_in_iso))
write.csv(summ, file.path(out_dir, "filter_summary.csv"), row.names = FALSE)
cat("\nSUMMARY\n"); print(summ)
cat("\ngenes in gene matrix:", nrow(g$counts),
    " | genes represented in isoform matrix:", g_in_iso, "\n")
cat("\nDONE\n")
