library(readr)
library(dplyr)


cpm <- read_tsv('/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/bambu/RGC_FT_extended_annotation/CPM_transcript.txt',
                col_names = TRUE)

#filter transcripts that start with ENST 
ENST_transcripts <- cpm[grepl("^ENST", cpm$TXNAME),]

ENST_transcripts <- ENST_transcripts |>
  dplyr::select(TXNAME, `H9-FT_1_primary_over_30_chr_only_sorted` )

#remove version number from transcript
ENST_transcripts$ensembl_transcript_id <- gsub("\\..*", "", ENST_transcripts$TXNAME)


require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")

mart <- useDataset("hsapiens_gene_ensembl", mart)
listAttributes(mart)[1:20,]


annotLookup <- getBM(
  mart=mart,
  attributes=c( "transcript_length",
                "ensembl_transcript_id"),
  filter="ensembl_transcript_id",
  values=ENST_transcripts$ensembl_transcript_id,
  uniqueRows=TRUE)

colnames(ENST_transcripts)[2] <- "H9-FT_1"
#merge the two dataframes
ENST_transcripts <- merge(ENST_transcripts, annotLookup, by = "ensembl_transcript_id")

nrow(ENST_transcripts)



## remove rows where count is 0
ENST_transcripts <- ENST_transcripts[ENST_transcripts$`H9-FT_1` > 0,]
nrow(ENST_transcripts)
#transcripts less than 1000bp

short_transcripts <- ENST_transcripts[ENST_transcripts$transcript_length < 1000,]
nrow(short_transcripts)

write.table(short_transcripts, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/FT_vs_RGC/short_transcripts.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
#transcripts greater than 1000bp < 5000bp
medium_transcripts <- ENST_transcripts[ENST_transcripts$transcript_length >= 1000 & ENST_transcripts$transcript_length < 5000,]

nrow(medium_transcripts)
write.table(medium_transcripts, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/FT_vs_RGC/medium_transcripts.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#long transcripts greater than 5000bp
long_transcripts <- ENST_transcripts[ENST_transcripts$transcript_length >= 5000,]
nrow(long_transcripts)
write.table(long_transcripts, file = "/users/sparthib/retina_lrs/processed_data/dtu/IsoformSwitchAnalyzeR/bambu/FT_vs_RGC/long_transcripts.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)






