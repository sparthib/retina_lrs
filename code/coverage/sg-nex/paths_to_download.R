library(readr)
library(dplyr)

sample_data <- read_tsv("/dcs04/hicks/data/sparthib/sg_nex_data/samples.tsv")
colnames(sample_data) 
direct_cDNA_data  <- sample_data |> filter(protocol == "direct-cDNA")


direct_cDNA_data$transcriptome_bam_path
#save the paths as a file
writeLines(direct_cDNA_data$transcriptome_bam_path, "/dcs04/hicks/data/sparthib/sg_nex_data/transcriptome_bam_paths_to_download.txt")

#download the files using wget
# wget -i transcriptome_bam_paths_to_download.txt -P /dcs04/hicks/data/sparthib/sg_nex_data/transcriptome_bam_files

