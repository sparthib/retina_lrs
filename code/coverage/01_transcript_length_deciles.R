library(readr)
library(dplyr)
library(biomaRt)

sample <- commandArgs(trailingOnly = TRUE)

read_asgts <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/read_asgts_for_gene_body.tsv",
                         sep="\t", header=FALSE, skip=3)
colnames(read_asgts) <- c("read_id", "transcript_id", "gene_id")

sample_reads <- read.table(paste0("/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/", 
                                  sample, "_read_ids.txt"))

# Filter out reads in read_asgts if they are in sample_reads
read_asgts <- read_asgts |> filter(read_id %in% sample_reads$V1)

# Remove version number in read_asgts$transcript_id
read_asgts$transcript_id <- gsub("\\..*", "", read_asgts$transcript_id)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "transcript_length"),
  filter="ensembl_transcript_id",
  values=read_asgts$transcript_id,
  uniqueRows=TRUE)

colnames(annotLookup) <- c("transcript_id", "transcript_length")

read_asgts <- read_asgts |> left_join(annotLookup, by="transcript_id")
read_asgts <- read_asgts |> distinct()

# Remove reads where transcript length is NA 
read_asgts_new <- read_asgts |> filter(!is.na(transcript_length))

write.table(read_asgts_new, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/reads_and_transcript_lengths_sample.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

read_asgts_new <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/reads_and_transcript_lengths_sample.tsv",
                             sep="\t", header=TRUE)

read_asgts_new <- read_asgts_new |> filter(transcript_length <= 10000)
hist(read_asgts_new$transcript_length, breaks=100, 
     main="Transcript Length Distribution", 
     xlab="Transcript Length", ylab="Frequency")

# Group reads by transcript_length deciles
read_asgts_new <- read_asgts_new |>
  mutate(decile = ntile(transcript_length, 10)) |>
  group_by(decile) 

dir.create(paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/", sample, "_deciles"), showWarnings = FALSE)

# Separate df by decile and save as a separate tsv 
for (i in 1:10) {
  read_asgts_new_decile <- read_asgts_new |> filter(decile == i) 
  write.table(read_asgts_new_decile, paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/",
                                            sample, "_deciles/", "decile_", i, ".tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}



