#!/usr/bin/env Rscript
# 08d_length_gc_distributions.R
# Run with: module load conda_R/4.5.x
#
# Distributions of gene GC content and transcript length for genes
# UP- vs DOWN-regulated in long reads, at each matched platform-DE stage.
#
# Significance / direction (user ruling): |logFC| >= 1 AND FDR < 0.05.
#   up_in_LR   = logFC >= 1  & FDR < 0.05   (higher in long read)
#   down_in_LR = logFC <= -1 & FDR < 0.05   (higher in short read)
# logFC sign: platform factor levels c("SR","LR"), coef platformLR ->
#   positive logFC = higher in long read.
#
# Stage binning:
#   Stage 1 (W_S1): LR D45  vs SR D25 + D65
#   Stage 2 (W_S2): LR D100 vs SR D65 + D100
#   Stage 3 (W_S3): LR D200 vs SR D180 + D280

.libPaths("/users/sparthib/R/4.5.x")
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})

outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/10_short_reads/08b_platform_de"
stages <- c(W_S1="Stage 1 (LR D45 vs SR D25+D65)",
            W_S2="Stage 2 (LR D100 vs SR D65+D100)",
            W_S3="Stage 3 (LR D200 vs SR D180+D280)")

LFC <- 1; FDR_CUT <- 0.05
all_df <- list(); summ <- data.frame()
for (w in names(stages)) {
  de <- read.csv(file.path(outdir, paste0("08b_platform_DE_stage_", w, ".csv")), stringsAsFactors=FALSE)
  de$direction <- with(de, ifelse(FDR < FDR_CUT & logFC >=  LFC, "Up in long read",
                           ifelse(FDR < FDR_CUT & logFC <= -LFC, "Down in long read", "Not sig")))
  de$stage <- stages[[w]]
  keep <- de[de$direction != "Not sig" & !is.na(de$median_gc_pct), ]
  all_df[[w]] <- keep
  for (d in c("Up in long read","Down in long read")) {
    s <- keep[keep$direction==d, ]
    summ <- rbind(summ, data.frame(stage=stages[[w]], direction=d, n=nrow(s),
                  median_gc=round(median(s$median_gc_pct),2),
                  median_len=round(median(s$median_tx_length),1)))
  }
}
df <- bind_rows(all_df)
df$direction <- factor(df$direction, levels=c("Up in long read","Down in long read"))
df$stage <- factor(df$stage, levels=unname(stages))
write.csv(summ, file.path(outdir, "08d_length_gc_distribution_summary.csv"), row.names=FALSE)
cat("Counts / medians:\n"); print(summ)

cols <- c("Up in long read"="#D7301F", "Down in long read"="#2166AC")

# --- GC content: density by direction, faceted by stage ---
p_gc <- ggplot(df, aes(median_gc_pct, fill=direction, colour=direction)) +
  geom_density(alpha=0.35, linewidth=0.5) +
  facet_wrap(~stage, nrow=1) +
  scale_fill_manual(values=cols) + scale_colour_manual(values=cols) +
  labs(x="Median gene GC content (%)", y="Density", fill=NULL, colour=NULL,
       title="A. GC content of platform-DE genes (|logFC|>=1, FDR<0.05)") +
  theme_bw(base_size=10) + theme(legend.position="top", panel.grid.minor=element_blank())

# --- Transcript length: density on log10 x, faceted by stage ---
p_len <- ggplot(df, aes(median_tx_length, fill=direction, colour=direction)) +
  geom_density(alpha=0.35, linewidth=0.5) +
  facet_wrap(~stage, nrow=1) +
  scale_x_log10() +
  scale_fill_manual(values=cols) + scale_colour_manual(values=cols) +
  labs(x="Median transcript length (bp, log10)", y="Density", fill=NULL, colour=NULL,
       title="B. Transcript length of platform-DE genes (|logFC|>=1, FDR<0.05)") +
  theme_bw(base_size=10) + theme(legend.position="top", panel.grid.minor=element_blank())

fig <- p_gc / p_len
ggsave(file.path(outdir, "08d_length_gc_distributions.png"), plot=fig,
       width=11, height=7, dpi=200)
cat("SAVED 08d_length_gc_distributions.png\nDONE\n")
