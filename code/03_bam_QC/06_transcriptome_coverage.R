library(dplyr)

sample <- commandArgs(trailingOnly = TRUE)


coverage_file <- "/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/transcriptome/"
coverage_file <- paste0(coverage_file, sample, ".coverage.txt" )

coverage_table <- read.table(coverage_file, 
                             header = FALSE, sep = "\t")

# Column names 
colnames(coverage_table) <- c("transcript", "position", "coverage")
head(coverage_table)

#check dimensions
dim(coverage_table)

#Group by transcript and remove groups where all coverage is 0, 
#and number of rows less than 50
coverage_table <- coverage_table |>
  group_by(transcript) |>
  filter(sum(coverage) > 0) |>
  filter(n() > 50)


#For each group, plot the coverage using ggplot2 
#where x-axis is position and y-axis is coverage
#label each plot by transcript name
#save all plots as one pdf file
#don't use facet_wrap or facet_grid

pdf(paste0( "/users/sparthib/retina_lrs/", sample, "_coverage_plots.pdf"))

for (i in 1:nrow(coverage_table)) {
  plot <- ggplot(coverage_table[i,], aes(x = position, y = coverage)) +
    geom_line() +
    ggtitle(coverage_table[i,]$transcript)
  print(plot)
}

dev.off()


