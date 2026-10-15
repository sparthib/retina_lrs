#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=ph_cov_merge
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=ph_cov_merge_%j.log
#
# 09_phased_hets_coverage_merge.R (sbatch wrapper)
# ---------------------------------------------------------------------------
# Merge the 11 per-BAM .depth files (from 08) onto the full phased-het BED
# (positions absent from a BAM -> depth 0), annotate each position with
# genomic location (exonic / intronic / intergenic) and overlapping gene via
# the release-46 GTF, and write the per-BAM coverage table.
#
# INPUT : /users/sparthib/.depth_arr/*.depth , /users/sparthib/phased_hets.bed
#         /dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
# OUTPUT: /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz
#         cols: chrom pos loc gene_id gene_name + 11 per-sample depth columns
#
# Run after 08, before 10.
# ---------------------------------------------------------------------------
set -eo pipefail
module load conda_R/4.5.x 2>/dev/null
Rscript - <<'RS'
.libPaths(c("/users/sparthib/R/4.5.x","/users/sparthib/R/4.4.x",.libPaths()))
suppressMessages({library(rtracklayer);library(GenomicRanges);library(data.table)})
arrdir  <- "/users/sparthib/.depth_arr"
gtf     <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
outfile <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz"
samples <- c("H9-BRN3B_hRO_2","H9-BRN3B-RO","H9-CRX_hRO_2","H9-CRX_ROs_D45","H9-FT_1","H9-FT_2","H9-hRGC_1","H9-hRGC_2","EP1-BRN3B-RO","EP1-WT_hRO_2","EP1-WT_ROs_D45")
colnm   <- gsub("-","_",samples)
bed <- fread("/users/sparthib/phased_hets.bed", header=FALSE, col.names=c("chrom","start","pos"))
bed[, start:=NULL]; setkey(bed, chrom, pos)
dt <- copy(bed)
for(i in seq_along(samples)){
  d <- fread(file.path(arrdir, paste0(samples[i], ".depth")), header=FALSE, col.names=c("chrom","pos","depth"))
  setkey(d, chrom, pos); dt <- d[dt]; setnames(dt, "depth", colnm[i]); dt[is.na(get(colnm[i])), (colnm[i]):=0]
}
g <- import(gtf); genes <- g[g$type=="gene"]; exons <- g[g$type=="exon"]
vg <- GRanges(dt$chrom, IRanges(dt$pos, width=1))
dt[, loc := ifelse(overlapsAny(vg,exons),"exonic", ifelse(overlapsAny(vg,genes),"intronic","intergenic"))]
h <- findOverlaps(vg, genes, select="first")
dt[, gene_id   := ifelse(is.na(h), NA, genes$gene_id[h])]
dt[, gene_name := ifelse(is.na(h), NA, genes$gene_name[h])]
setcolorder(dt, c("chrom","pos","loc","gene_id","gene_name", colnm))
fwrite(dt, outfile, sep="\t")
cat("rows:", nrow(dt), "-> ", outfile, "\n")
RS
