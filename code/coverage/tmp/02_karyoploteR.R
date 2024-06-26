#####load libraries ######

library(karyoploteR)
library(GenomicAlignments)
library(GenomicRanges)
library(readr)
library(purrr)

##### load bam dir and gene of interest files #####

bam_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE/YZ-15T_hRGC/sorted"

genes <- read_tsv("/users/sparthib/retina_lrs/code/coverage/coverage_genes_of_interest.tsv")

head(genes)


###### custom karyploteR function for plotting coverage from bam ######

# toGRanges("chr4:325000-400000"))

create_karyoplote <- function(gene, chr_num, chr_region){
  
  pdf(file = paste0("/users/sparthib/retina_lrs/plots/plotCoverage/karyoploteR_renamed/", gene, "_chr_", chr_region, ".pdf"))
  zoom.region <- toGRanges(paste0("chr",chr_region) )
  kp <- plotKaryotype(genome = "hg38", chromosomes = chr_num, zoom=zoom.region)
  kp <- kpAddBaseNumbers(kp, tick.dist = 1e4, add.units = TRUE)
  kp <- kpPlotBAMCoverage(kp, data = paste0(bam_dir, "/chr", chr_num , ".bam_sorted.bam"), max.valid.region.size = 2e6) 
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
  kpAddMainTitle(kp, main=paste0(gene, " chr", chr_region), cex=0.5)
  print(kp)
  dev.off()
}

##### use apply to plot cov for gene of interest #####

selected_columns <- genes |> dplyr::select("Gene Symbol", "Chr", "Chr:region")
selected_columns$Chr <- gsub(' ', '', selected_columns$Chr)
selected_columns$`Chr:region` <- gsub(' ', '', selected_columns$`Chr:region`)
selected_columns$`Chr:region` <- gsub(',', '', selected_columns$`Chr:region`)
selected_columns$`Chr:region` <- gsub('\n', '', selected_columns$`Chr:region`)
selected_columns <- as.matrix(selected_columns)
apply(selected_columns, 1, function(row)  create_karyoplote(row[1], row[2], row[3]))


##### TODO ######
#1. include gene name in argument for title
#2. fix y axis
#3. clean dataframe (remove spaces and commas)
#4. colors? 

#### TEMP CODE ######
# size <- 
#   
# What to insert in shell script: 
# chr_num 
# bam file
# size 
# BAM Coverage 
#   
  
# kp <- plotKaryotype(genome = "hg38", chromosomes = "$chr_#", zoom=toGRanges("$InsertedRange"))
# kp <- 
# 
# 
# kp <- plotKaryotype(genome = "dm6", chromosomes = "chr4", zoom=toGRanges("chr4:325000-400000"))
# kpAddBaseNumbers(kp, tick.dist = 25000, add.units = TRUE)
# kp <- kpPlotBAMCoverage(kp, data=bam1, col="gold", border="red", density=20)
# kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)