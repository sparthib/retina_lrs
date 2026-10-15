#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=ph_cov_bars
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=ph_cov_bars_%j.log
#
# 10_phased_hets_coverage_bars.sh
# ---------------------------------------------------------------------------
# Per-sample read-depth histograms (true bars, so depths with no sites stay
# empty) at phased hets, faceted sample (11 rows) x location (3 cols).
# Intronic shown 1-36, intergenic 1-17 (their ~99th pcts), exonic full range
# on a log10 x (spans 1 to ~53k). Counts on y, free y per sample.
#
# INPUT : /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz   (from 09)
# OUTPUT: /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_coverage_bars.png
# ---------------------------------------------------------------------------
set -eo pipefail
module load conda_R/4.5.x 2>/dev/null
Rscript - <<'RS'
.libPaths(c("/users/sparthib/R/4.5.x","/users/sparthib/R/4.4.x",.libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
tab <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_perBAM_coverage.tsv.gz"
png <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/phased_hets_coverage_bars.png"
d <- fread(tab)
samples <- c("H9_BRN3B_hRO_2","H9_BRN3B_RO","H9_CRX_hRO_2","H9_CRX_ROs_D45","H9_FT_1","H9_FT_2","H9_hRGC_1","H9_hRGC_2","EP1_BRN3B_RO","EP1_WT_hRO_2","EP1_WT_ROs_D45")
m <- melt(d, id.vars=c("chrom","pos","loc"), measure.vars=samples, variable.name="sample", value.name="depth")
m <- m[depth>=1]; m[, sample := factor(sample, levels=samples)]
pal <- c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#ff7f00","#fdbf6f","#6a3d9a","#cab2d6","#b15928"); names(pal) <- samples
mk <- function(dd, brk, logx, ttl){
  p <- ggplot(dd, aes(x=depth, fill=sample))
  if(logx) p <- p + geom_histogram(bins=35, colour=NA) + scale_x_log10(breaks=c(1,10,100,1000,10000), labels=c("1","10","100","1k","10k"))
  else     p <- p + geom_histogram(binwidth=1, colour="grey30", linewidth=0.1) + scale_x_continuous(breaks=brk)
  p + facet_grid(sample ~ ., scales="free_y", switch="y") + scale_fill_manual(values=pal, guide="none") +
    labs(x=ttl, y=NULL) + theme_bw(base_size=9) +
    theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          strip.text.y.left=element_text(angle=0,size=6.5), strip.background=element_rect(fill="grey94"),
          axis.text.y=element_text(size=5.5))
}
p_in <- mk(m[loc=="intronic" & depth<=36], seq(0,36,6), FALSE, "Intronic: read depth (1-36)")
p_ig <- mk(m[loc=="intergenic" & depth<=17], seq(0,17,3), FALSE, "Intergenic: read depth (1-17)")
p_ex <- mk(m[loc=="exonic"], NULL, TRUE, "Exonic: read depth (log, full range to ~53k)")
P <- (p_in | p_ig | p_ex) + plot_annotation(title="Per-sample read-depth histograms at WhatsHap-phased hets (bars = site counts)", theme=theme(plot.title=element_text(size=11)))
ggsave(png, plot=P, width=13, height=11, dpi=150); cat("WROTE", png, "\n")
RS
