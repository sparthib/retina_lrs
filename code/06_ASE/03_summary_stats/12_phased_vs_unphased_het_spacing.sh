#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=jhpce
#SBATCH --job-name=het_spacing
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=het_spacing_%j.log
#
# 12_phased_vs_unphased_het_spacing.sh
# ---------------------------------------------------------------------------
# Compare genomic spacing of phased vs unphased heterozygous sites, to show
# that unphased hets are more isolated (farther from a neighbouring het, so
# less likely to share a read that WhatsHap could use to link them).
#
#   Phased hets   -> distance to nearest OTHER PHASED het (what WhatsHap could
#                    actually link against)
#   Unphased hets -> distance to nearest het of ANY kind (best-case partner)
#
# INPUT : all_samples_H9_and_EP1_phased.vcf.gz   (phased VCF)
# OUTPUT: <OUT>/phased_vs_unphased_het_spacing.png
#         <OUT>/phased_vs_unphased_het_spacing.tsv.gz  (dist, set)
# ---------------------------------------------------------------------------
set -eo pipefail
VCF=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output/all_samples_H9_and_EP1_phased.vcf.gz
OUT=/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output
module load htslib bcftools 2>/dev/null
bcftools query -f '%CHROM\t%POS\t[%GT]\n' $VCF \
  | awk -F'\t' '$3=="0|1"||$3=="1|0"||$3=="0/1"||$3=="1/2"' > /users/sparthib/.hets_gt.tsv
module load conda_R/4.5.x 2>/dev/null
Rscript - <<'RS'
.libPaths(c("/users/sparthib/R/4.5.x","/users/sparthib/R/4.4.x",.libPaths()))
suppressMessages({library(data.table); library(ggplot2)})
OUT <- "/dcs04/hicks/data/sparthib/retina_lrs/09_ASE/H9_DNA_Seq_data/whatshap_output"
h <- fread("/users/sparthib/.hets_gt.tsv", header=FALSE, col.names=c("chrom","pos","gt"))
h[, phased := gt %in% c("0|1","1|0")]; setorder(h, chrom, pos)
# phased -> nearest other phased
ph <- h[phased==TRUE, .(chrom,pos)]; setorder(ph,chrom,pos)
ph[, d_prev := pos - shift(pos,1,type="lag"),  by=chrom]
ph[, d_next := shift(pos,1,type="lead") - pos, by=chrom]
ph[, nn := pmin(d_prev,d_next,na.rm=TRUE)]; ph <- ph[is.finite(nn)]
# unphased -> nearest any het
allpos <- h[, .(chrom,pos)]; setorder(allpos,chrom,pos)
allpos[, ap_prev := pos - shift(pos,1,type="lag"),  by=chrom]
allpos[, ap_next := shift(pos,1,type="lead") - pos, by=chrom]
allpos[, nn_any := pmin(ap_prev,ap_next,na.rm=TRUE)]
un <- merge(h[phased==FALSE,.(chrom,pos)], allpos[,.(chrom,pos,nn_any)], by=c("chrom","pos"))
un <- un[is.finite(nn_any)]
cat("PHASED n:", nrow(ph), " UNPHASED n:", nrow(un), "\n")
cat("PHASED quantiles:\n");   print(round(quantile(ph$nn, c(.1,.25,.5,.75,.9,.99))))
cat("UNPHASED quantiles:\n"); print(round(quantile(un$nn_any, c(.1,.25,.5,.75,.9,.99))))
pl <- rbind(data.table(dist=ph$nn, set="Phased het \u2192 nearest phased het"),
            data.table(dist=un$nn_any, set="Unphased het \u2192 nearest het (any)"))[dist>=1]
p <- ggplot(pl, aes(x=dist, colour=set)) + geom_density(linewidth=0.9) +
  scale_x_log10(breaks=c(1,10,100,1000,1e4,1e5,1e6), labels=c("1","10","100","1k","10k","100k","1M")) +
  scale_colour_manual(values=c("Phased het \u2192 nearest phased het"="#1f78b4",
                               "Unphased het \u2192 nearest het (any)"="#e31a1c"), name=NULL) +
  labs(x="distance to nearest neighbouring het (bp, log scale)", y="density",
       title="Spacing of phased vs unphased heterozygous sites") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), legend.position=c(0.72,0.85),
        legend.background=element_rect(fill=alpha("white",0.7), colour=NA),
        plot.title=element_text(size=11,hjust=0))
ggsave(file.path(OUT,"phased_vs_unphased_het_spacing.png"), plot=p, width=8, height=4.6, dpi=150)
fwrite(pl, file.path(OUT,"phased_vs_unphased_het_spacing.tsv.gz"), sep="\t")
cat("WROTE png + tsv\n")
RS
echo "=== DONE ==="
