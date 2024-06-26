library(readr)
library(dplyr)
library(biomaRt)

read_asgts <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/read_asgts_for_gene_body.tsv",
                       sep="\t", header=FALSE,
                       skip=3)
colnames(read_asgts) <- c("read_id", "transcript_id", "gene_id")

H9_FT_2_reads <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_2_read_ids.txt")
#filter out reads in read_asgts if they are in H9_FT_2_reads
read_asgts <- read_asgts |> filter(read_id %in% H9_FT_2_reads$V1)
#remove version number in read_asgts$transcript_id
read_asgts$transcript_id <- gsub("\\..*", "", read_asgts$transcript_id)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c( "ensembl_transcript_id",
                 "transcript_length"),
  filter="ensembl_transcript_id",
  values=read_asgts$transcript_id,
  uniqueRows=TRUE)

colnames(annotLookup) <- c("transcript_id", "transcript_length")

read_asgts <- read_asgts |> left_join (annotLookup, by="transcript_id")
read_asgts <- read_asgts |> distinct()

#remove reads where transcript length is NA 
read_asgts_new <- read_asgts |> filter(!is.na(transcript_length))
nrow(read_asgts_new)
#9736802

#get list of transcripts that have missing transcript length
read_asgts_missing <- read_asgts |> filter(is.na(transcript_length)) |> distinct()
nrow(read_asgts_missing)

range(read_asgts_new$transcript_length)
# range(read_asgts_new$transcript_length)
# [1]     48 347561

write.table(read_asgts_new, "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/reads_and_transcript_lengths_H9_FT_2.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

#######
library(dplyr)
read_asgts_new <- read.table("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/reads_and_transcript_lengths_H9_FT_2.tsv",
                             sep="\t", header=TRUE)

read_asgts_new <- read_asgts_new |> dplyr::filter(transcript_length <= 10000)
hist(read_asgts_new$transcript_length, breaks=100, 
     main="Transcript Length Distribution", 
     xlab="Transcript Length", ylab="Frequency")

#group reads by transcript_length deciles
read_asgts_new <- read_asgts_new |>
  mutate(decile = ntile(transcript_length, 10)) |>
  group_by(decile) 

dir.create("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/H9_FT2_deciles", showWarnings = FALSE)
# View the result
# separate df by decile and save as a separate tsv 
for (i in 1:10) {
  read_asgts_new_decile <- read_asgts_new |> filter(decile == i) 
  write.table(read_asgts_new_decile, paste0("/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/H9_FT2_deciles/decile_", i, ".tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}

#check number of reads in each decile
read_asgts_new |> count(decile)



