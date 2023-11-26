
library(readr)
library(dplyr)
library(tidyr)
library("sessioninfo")


sample_num <- as.numeric(Sys.getenv(SLURM_ARRAY_TASK_ID))

# colnames(transcript_lengths) = c("tx_id", "gene_id", "nexons", "length")
# nrow(transcript_lengths)
# write_tsv(transcript_lengths, "/dcs04/hicks/data/sparthib/transcript_lengths_sorted.tsv",
#           append=FALSE)

config <- read_tsv("/users/sparthib/retina_lrs/config.tsv", col_names = TRUE)

sample <- config[sample_num , 2]
IsoQuant_dir <- paste0("/dcs04/hicks/data/sparthib/casey/IsoQuant_output/", sample,"/OUT")
IsoQuant_tpm <- readr::read_tsv(paste0(IsoQuant_dir, "/OUT.transcript_tpm.tsv"),
                                col_names = T)
IsoQuant_counts <-  readr::read_tsv( paste0(IsoQuant_dir, "/OUT.transcript_counts.tsv"),
                                     col_names = T)
transcript_lengths <- readr::read_tsv("/dcs04/hicks/data/sparthib/transcript_lengths_sorted.tsv",
                                      col_names = TRUE)
# nrow(transcript_lengths)
# [1] 117748

tpm_counts <- merge(IsoQuant_tpm, IsoQuant_counts, by = "#feature_id")


# nrow(IsoQuant_tpm)
nrow(IsoQuant_tpm)
#  nrow(IsoQuant_counts)
nrow(IsoQuant_tpm)
# nrow(tpm_counts)
nrow(tpm_counts)

tx_length_counts <- merge(tpm_counts, transcript_lengths, by.x="#feature_id",
                          by.y="tx_id", all.x = TRUE)

tx_length_counts$sample_name <- sample

# nrow(tx_length_counts)
nrow(tx_length_counts)

write_tsv(tx_length_counts, paste0("/dcs04/hicks/data/sparthib/casey/diff_expression_data/transcript_lengths/", 
                 sample,"_tx_length_counts.tsv" ), append=FALSE)

#head(tx_length_counts)
head(tx_length_counts)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

