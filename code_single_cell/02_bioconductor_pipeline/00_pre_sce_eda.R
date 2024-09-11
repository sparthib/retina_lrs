library(stringr)



#get all genome and transcriptome flagstat files 

genome_files <- list.files("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/",
                           pattern = "genome_bam_flagstat.txt",
                           recursive = T,
                           full.names = T)

transcriptome_files <- list.files("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/",
                                  pattern = "transcriptome_bam_flagstat.txt",
                                  recursive = T,
                                  full.names = T)

#create function that reads file, grep only line with "primary mapped" and return the number of reads aligned

get_reads <- function(file) {
  line <- readLines(file)[grep("primary mapped", readLines(file))]
  num_reads <- str_extract(line, "\\d+")
  as.numeric(num_reads)
  
}

#apply function to genome files 
genome_reads <- sapply(genome_files, get_reads)

names(genome_reads) <- gsub(names(genome_reads), 
                            pattern = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output//", replacement = "")
names(genome_reads) <- gsub(names(genome_reads), pattern = "/genome_bam_flagstat.txt", replacement = "")


#repeat for transcriptome files 

transcriptome_reads <- sapply(transcriptome_files, get_reads)

names(transcriptome_reads) <- gsub(names(transcriptome_reads), 
                            pattern = "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output//", replacement = "")
names(transcriptome_reads) <- gsub(names(transcriptome_reads), pattern = "/transcriptome_bam_flagstat.txt", replacement = "")


transcriptome_reads


#combine into dataframe
primary_mapped_reads_df <- data.frame(genome = genome_reads, transcriptome = transcriptome_reads)

#Primary Unmapped Reads
#                   genome transcriptome
# 10x_D100-EP1_A1   18117032       9797100
# 10x_D100-EP1_A2     699749        407313
# 10x_D100-EP1_B1   32457291      16332639
# 10x_D100-EP1_B2    3874547       2247286
# 10x_D200-EP1-1_A1 31997316      10797588
# 10x_D200-EP1-1_A2  2519136       1328475
# 10X_D200-EP1-1_B1 35374512      22593735
# 10x_D200-EP1-1_B2  2136278       1106423
# 10x_D200-EP1-2_A1 19896090       6927127
# 10x_D200-EP1-2_A2   539154        289590
# 10x_D200-EP1-2_B1 20918783       7094870
# 10x_D200-EP1-2_B2   654685        349144




