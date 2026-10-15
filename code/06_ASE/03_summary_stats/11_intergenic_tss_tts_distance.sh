#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=ig_tss_tts
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=ig_tss_tts_%j.log
#
# 11_intergenic_tss_tts_distance.sh
# ---------------------------------------------------------------------------
# For each INTERGENIC phased het, compute genomic distance to the nearest
# canonical TSS and nearest TTS (strand-aware gene ends from the release-46
# GTF), take the closer one, and plot the distance distribution colored by
# which boundary (TSS vs TTS) is nearest.
#
# INPUT : /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz  (loc column; from 09)
#         /dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
# OUTPUT: /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/intergenic_hets_tss_tts_distance.tsv.gz
#           cols: chrom pos dTSS dTTS nearest_type nearest_dist
#         /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/intergenic_hets_tss_tts_distance.png
# ---------------------------------------------------------------------------
set -eo pipefail
module load conda_R/4.5.x 2>/dev/null
Rscript - <<'RS'
.libPaths(c("/users/sparthib/R/4.5.x","/users/sparthib/R/4.4.x",.libPaths()))
suppressMessages({library(rtracklayer);library(GenomicRanges);library(data.table);library(ggplot2)})
tab <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz"; gtf <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf"
d <- fread(tab, select=c("chrom","pos","loc")); ig <- d[loc=="intergenic"]
cat("intergenic hets:", nrow(ig), "\n")
g <- import(gtf); genes <- g[g$type=="gene"]; st <- as.character(strand(genes))
tss <- GRanges(seqnames(genes), IRanges(ifelse(st=="-", end(genes), start(genes)), width=1))
tts <- GRanges(seqnames(genes), IRanges(ifelse(st=="-", start(genes), end(genes)), width=1))
vg <- GRanges(ig$chrom, IRanges(ig$pos, width=1))
nt <- distanceToNearest(vg, tss, ignore.strand=TRUE); nq <- distanceToNearest(vg, tts, ignore.strand=TRUE)
ig[, dTSS := NA_integer_]; ig[queryHits(nt), dTSS := mcols(nt)$distance]
ig[, dTTS := NA_integer_]; ig[queryHits(nq), dTTS := mcols(nq)$distance]
ig[, nearest_type := ifelse(dTSS<=dTTS, "TSS", "TTS")]
ig[, nearest_dist := pmin(dTSS, dTTS, na.rm=TRUE)]
fwrite(ig[,.(chrom,pos,dTSS,dTTS,nearest_type,nearest_dist)], "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/intergenic_hets_tss_tts_distance.tsv.gz", sep="\t")
pl <- ig[nearest_dist>=1]
p <- ggplot(pl, aes(x=nearest_dist, colour=nearest_type)) + geom_density(linewidth=0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6), labels=c("1","10","100","1k","10k","100k","1M")) +
  scale_colour_manual(values=c(TSS="#1f78b4", TTS="#e31a1c"), name="nearest boundary") +
  labs(x="distance to nearest gene boundary (bp, log scale)", y="density",
       title=sprintf("Intergenic phased hets: distance to nearest canonical TSS/TTS (n=%s)", format(nrow(ig),big.mark=","))) +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), plot.title=element_text(size=10,hjust=0))
ggsave("/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/intergenic_hets_tss_tts_distance.png", plot=p, width=8, height=4.5, dpi=150)
cat("WROTE png + tsv\n")
RS
