library(dplyr)
library(ggplot2)

#load all samples in input folder that have suffix count_distribution.txt

input_dir <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/CIGAR_sequences"

files <- list.files(input_dir, pattern = "count_distribution.txt", full.names = TRUE)

#for each file, add a column of the sample name 
files <- lapply(files, function(x) { 
  sample_name <- gsub(".*/", "", x)
  sample_name <- gsub("_N_count_distribution.txt", "", sample_name)
  data <- read.table(x, header = TRUE, sep = ",")
  data$sample <- sample_name
  return(data)
})


# combine all files into one data frame
data <- do.call(rbind, files)

#group by sample and convert frequency to proportion 
unique(data$sample)
RO_data <- data |> filter(sample %in% c("EP1-BRN3B-RO","EP1-WT_hRO_2","EP1-WT_ROs_D45","H9-BRN3B_hRO_2","H9-BRN3B-RO" ,"H9-CRX_hRO_2","H9-CRX_ROs_D45"))
#remove Unique_Count = 0
RO_data <- RO_data |> filter(Unique_Count != 0)
#summarize proportion
RO_data <- RO_data |> group_by(Unique_Count) |> summarise(Frequency = sum(Frequency))  |> mutate(proportion = Frequency / sum(Frequency))

#plot the proportion for unique count = 1 to 10

pdf("/users/sparthib/retina_lrs/plots/exon_exon/RO_num_exon_junctions.pdf")
p <- RO_data |> filter(Unique_Count <= 10) |> ggplot(aes(x = Unique_Count, y = proportion)) + geom_bar(stat = "identity")
print(p)
dev.off()


FT_vs_RGC_data <- data |> filter(sample %in% c("H9-FT_1","H9-FT_2","H9-hRGC_1" ,"H9-hRGC_2" ))
#remove Unique_Count = 0
FT_vs_RGC_data <- FT_vs_RGC_data |> filter(Unique_Count != 0)
#summarize proportion
FT_vs_RGC_data <- FT_vs_RGC_data |> group_by(Unique_Count) |> summarise(Frequency = sum(Frequency))  |> mutate(proportion = Frequency / sum(Frequency))

pdf("/users/sparthib/retina_lrs/plots/exon_exon/FT_vs_RGC_num_exon_junctions.pdf")
p <- FT_vs_RGC_data |> filter(Unique_Count <= 10) |> ggplot(aes(x = Unique_Count, y = proportion)) + geom_bar(stat = "identity") 
print(p)
dev.off()


#all samples 
#remove Unique_Count = 0
non_zero_data <- data |> filter(Unique_Count != 0)
#summarize proportion
non_zero_data <- non_zero_data |> group_by(Unique_Count) |> summarise(Frequency = sum(Frequency))  |> mutate(proportion = Frequency / sum(Frequency))

pdf("/users/sparthib/retina_lrs/plots/exon_exon/all_samples_num_exon_junctions.pdf")
p <- non_zero_data |> filter(Unique_Count <= 10) |> ggplot(aes(x = Unique_Count, y = proportion)) + geom_bar(stat = "identity") 
print(p)
dev.off()

