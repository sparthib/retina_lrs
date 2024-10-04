



day100_b1_genome <- read.table("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10x_D100-EP1_B1/genome_read_lengths.txt",header = F)
day100_b1_genome$sample <- "day100_B1"
day100_b1_genome$alignment <- "genome"


day100_b2_genome <- read.table("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10x_D100-EP1_B2/genome_read_lengths.txt",header = F)

day100_b2_genome$sample <- "day100_B2"
day100_b2_genome$alignment <- "genome"

day100_b1_transcriptome <- read.table("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10x_D100-EP1_B1/transcript_read_lengths.txt",header = F)

day100_b1_transcriptome$sample <- "day100_B1"
day100_b1_transcriptome$alignment <- "transcriptome"

day100_b2_transcriptome <-read.table("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10x_D100-EP1_B2/transcript_read_lengths.txt",header = F)


day100_b2_transcriptome$sample <- "day100_B2"
day100_b2_transcriptome$alignment <- "transcriptome"


data <- rbind(day100_b1_genome,day100_b2_genome,day100_b1_transcriptome,day100_b2_transcriptome)

colnames(data) <- c("read_length","sample","alignment")


#make a histogram of the read lengths
pdf("/users/sparthib/retina_lrs/single_cell_plots/temp_read_length_bams.pdf")
p<-ggplot(data, aes(x = read_length)) +
  geom_histogram(binwidth = 100, color = "black", fill = "skyblue") +
  labs(title = "Histogram of Read Lengths by Sample and Alignment",
       x = "Read Length",
       y = "Frequency") +
  theme_minimal() +
  facet_grid(sample ~ alignment) + 
  xlim(0,25000) +
  ylim(0,500000)
print(p)
dev.off()



##########



pdf("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/transcriptome/01_quality_controlled/size_factors.pdf")

# 
# new_sce_list <- lapply(sce_list, function(x) {
#   x <- x[discard=FALSE]
#   cluster_x <- quickCluster(x) 
#   x <- computeSumFactors(x, cluster=cluster_x, min.mean=0.1)
#   print("computed sum factors")
#   x <- logNormCounts(x)
#   print(x$sample_replicate[1])
#   x
# })

sce_list <- lapply(sce_list, function(x) {
  x <- x[discard=FALSE]
  x
})

#10x_D100_EP1_B1
# Remove anything after the first dot
d100_b1_transcripts <- gsub("\\..*$", "", rownames(sce_list[[3]]))
d100_b2_transcripts <- gsub("\\..*$", "", rownames(sce_list[[4]]))



library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transcript_id_to_length <- function(transcript_id) {
  # Fetch gene symbols using getBM
  transcript_length <- getBM(attributes = c("ensembl_transcript_id", "transcript_length"),
                             filters = "ensembl_transcript_id",
                             values = transcript_id,
                             mart = mart)
  transcript_length
}

transcript_id_to_gc_content <- function(transcript_id){
  gene_gc_content <- getBM(attributes = c("ensembl_transcript_id", "percentage_gene_gc_content"),
                           filters = "ensembl_transcript_id",
                           values = transcript_id,
                           mart = mart)
  
}

transcript_length_d100_b1 <- transcript_id_to_length(d100_b1_transcripts)
transcript_length_d100_b1$sample_replicate <-  "10x_D100_EP1_B1"
transcript_length_d100_b2$sample_replicate <- "10x_D100_EP1_B2"

transcript_gc_content_d100_b1 <- transcript_id_to_gc_content(d100_b1_transcripts)
transcript_gc_content_d100_b2 <- transcript_id_to_gc_content(d100_b2_transcripts)

transcript_gc_content_d100_b1$sample_replicate <-  "10x_D100_EP1_B1"
transcript_gc_content_d100_b2$sample_replicate <- "10x_D100_EP1_B2"


transcript_lengths <- rbind(transcript_length_d100_b1, transcript_length_d100_b2)
transcript_gc_content <- rbind(transcript_gc_content_d100_b1, transcript_gc_content_d100_b2)

head(transcript_lengths)
#create histogram of transcript lengths bin by 100

pdf("/users/sparthib/retina_lrs/single_cell_plots/temp_transcript_lengths.pdf")
p <- ggplot(transcript_lengths, aes(x = transcript_length, y = ..count../sum(..count..)*100)) +  # Convert counts to percentages
  geom_histogram(binwidth = 100, color = "black", fill = "skyblue") +
  labs(title = "Histogram of Transcript Lengths by Sample Replicate",
       x = "Transcript Length",
       y = "Percentage") +
  theme_minimal() +
  facet_wrap(~ sample_replicate) 
print(p)
dev.off()


pdf("/users/sparthib/retina_lrs/single_cell_plots/temp_gc_content.pdf")
p <- ggplot(transcript_gc_content , aes(x = percentage_gene_gc_content)) +
  geom_histogram(binwidth = 10, color = "black", fill = "skyblue") +
  labs(title = "Histogram of Gene GC content of trancsript by Sample Replicate",
       x = "GC content",
       y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ sample_replicate)
print(p)
dev.off()

transcript_length_d100_b2 <- transcript_id_to_length(d100_b2_transcripts)

listAttributes(mart)[10:20,]
#get transcript length from biomart for transcript ensembl IDs


