library(tximport)
library(dplyr)
library(readr)
library('biomaRt')
library(EnsDb.Hsapiens.v79)
library(stringr)

# sample_num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

sample_num <- 4
config <- read_tsv("~/retina_lrs/config.tsv", col_names = TRUE)
sample <- config$sample_name[sample_num]

tx2gene <-  read_tsv("/dcs04/hicks/data/sparthib/transcript_lengths_sorted.tsv")

tx2gene <- tx2gene |> select(tx_id, gene_id)

quant_path <- paste0("/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level/", sample,
                     "/quant.sf")
txi <- tximport(quant_path,
                type = "salmon",
                tx2gene)

class(txi$counts)

txi_counts <- as.data.frame(txi$counts)
colnames(txi_counts) <- "count"
txi_counts_filtered <- txi_counts |> dplyr::filter(count != 0)
nrow(txi_counts_filtered)
#EP1-BRN3B-RO
#14935
# head(txi_counts)
# count
# ENSG00000000003.16 489.957
# ENSG00000000005.6    0.000
# ENSG00000000419.14 116.547


genes <- tx2gene$gene_id
get_part_before_dot <- function(vec) {
  parts <- sapply(strsplit(vec, "\\."), function(x) x[[1]])
  return(parts)
}

genes_truncated <- get_part_before_dot(genes)

G_list <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes_truncated , keytype = "GENEID", columns = c("SYMBOL","GENEID"))


txi_counts$genes_truncated <- get_part_before_dot(rownames(txi_counts))

txi_counts_symbols <- merge(txi_counts, G_list, by.x = "genes_truncated",
                            by.y = "GENEID", all.x = TRUE)

print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "NF1"))
print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "BSG"))

# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "NF1"))
# genes_truncated   count SYMBOL
# 1 ENSG00000196712 236.567    NF1
# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "BSG"))
# genes_truncated    count SYMBOL
# 1 ENSG00000172270 1511.756    BSG


##sample 2 "H9-BRN3B-RO"
nrow(txi_counts_filtered) 14825

# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "NF1"))
# genes_truncated   count SYMBOL
# 1 ENSG00000196712 123.375    NF1
# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "BSG"))
# genes_truncated    count SYMBOL
# 1 ENSG00000172270 1200.811    BSG


##sample 3 14572  "hRGC"

# genes_truncated   count SYMBOL
# 1 ENSG00000196712 457.086    NF1
# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "BSG"))
# genes_truncated    count SYMBOL
# 1 ENSG00000172270 1697.994    BSG

##sample 4 
# 14119
# genes_truncated   count SYMBOL
# 1 ENSG00000196712 361.314    NF1
# > print(txi_counts_symbols |>  dplyr::filter(SYMBOL == "BSG"))
# genes_truncated   count SYMBOL
# 1 ENSG00000172270 940.891    BSG


####Check which transcripts 

tx2gene <-  read_tsv("/dcs04/hicks/data/sparthib/transcript_lengths_sorted.tsv")

tx2gene <- tx2gene |> dplyr::select(tx_id, gene_id)


sample_num <- 4
config <- read_tsv("~/retina_lrs/config.tsv", col_names = TRUE)
sample <- config$sample_name[sample_num]


quant_path <- paste0("/dcs04/hicks/data/sparthib/casey/salmon_outputs_transcript_level/", sample,
                     "/quant.sf")
txi <- tximport(quant_path,
                type = "salmon",
                txOut = TRUE)

txi_counts <- as.data.frame(txi$counts)
txi_counts$tx <- rownames(txi$counts)
colnames(txi_counts) <- c("count", "tx_id")
nrow(txi_counts |> dplyr::filter(count != 0))

#number of isoforms expressed 
#sample 1  51218 EP1-BRN3B-RO
#sample 2  47318 H9-BRN3B-RO
#sample 3 52772 hRGC 
#sample 4 50761 DG-WT-hRGC


txi_counts <- merge(txi_counts, tx2gene, by.x = "tx_id", 
                    by.y = "tx_id", all.x = TRUE)

txi_counts_filtered <- na.omit(txi_counts)
txi_counts_filtered$genes_truncated <- get_part_before_dot(txi_counts_filtered$gene_id)
G_list <- ensembldb::select(EnsDb.Hsapiens.v79, keys= txi_counts_filtered$genes_truncated , 
                            keytype = "GENEID", columns = c("SYMBOL","GENEID"))

txi_counts_filtered <- merge(txi_counts_filtered, G_list, by.x = "genes_truncated",
                            by.y = "GENEID", all.x = TRUE)

txi_counts_filtered |> dplyr::filter(SYMBOL == "NF1")
txi_counts_filtered |> dplyr::filter(SYMBOL =="BSG")

##sample 1 EP1-BRN3B-RO
# > txi_counts_filtered |> dplyr::filter(SYMBOL == "NF1")
# genes_truncated             tx_id  count            gene_id SYMBOL
# 1  ENSG00000196712 ENST00000687027.1 42.521 ENSG00000196712.20    NF1
# 2  ENSG00000196712 ENST00000684826.1  2.522 ENSG00000196712.20    NF1
# 3  ENSG00000196712 ENST00000693617.1  0.000 ENSG00000196712.20    NF1
# 4  ENSG00000196712 ENST00000696138.1 21.462 ENSG00000196712.20    NF1
# 5  ENSG00000196712 ENST00000691014.1  2.988 ENSG00000196712.20    NF1
# 6  ENSG00000196712 ENST00000490416.2  0.000 ENSG00000196712.20    NF1
# 7  ENSG00000196712 ENST00000356175.7 39.039 ENSG00000196712.20    NF1
# 8  ENSG00000196712 ENST00000358273.9 62.137 ENSG00000196712.20    NF1
# 9  ENSG00000196712 ENST00000431387.8 11.321 ENSG00000196712.20    NF1
# 10 ENSG00000196712 ENST00000487476.5 54.577 ENSG00000196712.20    NF1
# > txi_counts_filtered |> dplyr::filter(SYMBOL =="BSG")
# genes_truncated             tx_id   count            gene_id SYMBOL
# 1  ENSG00000172270 ENST00000679472.1  16.568 ENSG00000172270.22    BSG
# 2  ENSG00000172270 ENST00000333511.9 439.238 ENSG00000172270.22    BSG
# 3  ENSG00000172270 ENST00000353555.9 787.662 ENSG00000172270.22    BSG
# 4  ENSG00000172270 ENST00000545507.6   3.192 ENSG00000172270.22    BSG
# 5  ENSG00000172270 ENST00000573784.6   0.000 ENSG00000172270.22    BSG
# 6  ENSG00000172270 ENST00000680552.1   0.000 ENSG00000172270.22    BSG
# 7  ENSG00000172270 ENST00000618006.4   2.977 ENSG00000172270.22    BSG
# 8  ENSG00000172270 ENST00000680326.1   2.764 ENSG00000172270.22    BSG
# 9  ENSG00000172270 ENST00000346916.9   0.000 ENSG00000172270.22    BSG
# 10 ENSG00000172270 ENST00000576984.3   2.657 ENSG00000172270.22    BSG
# 11 ENSG00000172270 ENST00000614867.2   0.000 ENSG00000172270.22    BSG
# 12 ENSG00000172270 ENST00000680065.1 256.698 ENSG00000172270.22    BSG


##sample 2 "H9-BRN3B-RO"

# genes_truncated             tx_id  count            gene_id SYMBOL
# 1  ENSG00000196712 ENST00000687027.1  9.547 ENSG00000196712.20    NF1
# 2  ENSG00000196712 ENST00000684826.1  0.000 ENSG00000196712.20    NF1
# 3  ENSG00000196712 ENST00000693617.1  0.000 ENSG00000196712.20    NF1
# 4  ENSG00000196712 ENST00000696138.1 13.707 ENSG00000196712.20    NF1
# 5  ENSG00000196712 ENST00000691014.1  8.162 ENSG00000196712.20    NF1
# 6  ENSG00000196712 ENST00000490416.2  0.000 ENSG00000196712.20    NF1
# 7  ENSG00000196712 ENST00000356175.7 23.809 ENSG00000196712.20    NF1
# 8  ENSG00000196712 ENST00000358273.9 19.878 ENSG00000196712.20    NF1
# 9  ENSG00000196712 ENST00000431387.8  5.177 ENSG00000196712.20    NF1
# 10 ENSG00000196712 ENST00000487476.5 43.095 ENSG00000196712.20    NF1
# > txi_counts_filtered |> dplyr::filter(SYMBOL =="BSG")
# genes_truncated             tx_id   count            gene_id SYMBOL
# 1  ENSG00000172270 ENST00000679472.1   2.852 ENSG00000172270.22    BSG
# 2  ENSG00000172270 ENST00000333511.9 412.097 ENSG00000172270.22    BSG
# 3  ENSG00000172270 ENST00000353555.9 644.236 ENSG00000172270.22    BSG
# 4  ENSG00000172270 ENST00000545507.6   0.000 ENSG00000172270.22    BSG
# 5  ENSG00000172270 ENST00000573784.6   7.822 ENSG00000172270.22    BSG
# 6  ENSG00000172270 ENST00000680552.1   0.000 ENSG00000172270.22    BSG
# 7  ENSG00000172270 ENST00000618006.4   5.545 ENSG00000172270.22    BSG
# 8  ENSG00000172270 ENST00000680326.1   0.000 ENSG00000172270.22    BSG
# 9  ENSG00000172270 ENST00000346916.9   1.040 ENSG00000172270.22    BSG
# 10 ENSG00000172270 ENST00000576984.3   2.157 ENSG00000172270.22    BSG
# 11 ENSG00000172270 ENST00000614867.2   1.473 ENSG00000172270.22    BSG
# 12 ENSG00000172270 ENST00000680065.1 123.589 ENSG00000172270.22    BSG


## sample 3 "hRGC"

# genes_truncated             tx_id   count            gene_id SYMBOL
# 1  ENSG00000196712 ENST00000687027.1  38.006 ENSG00000196712.20    NF1
# 2  ENSG00000196712 ENST00000684826.1   2.422 ENSG00000196712.20    NF1
# 3  ENSG00000196712 ENST00000693617.1   0.000 ENSG00000196712.20    NF1
# 4  ENSG00000196712 ENST00000696138.1  47.718 ENSG00000196712.20    NF1
# 5  ENSG00000196712 ENST00000691014.1  28.739 ENSG00000196712.20    NF1
# 6  ENSG00000196712 ENST00000490416.2   0.000 ENSG00000196712.20    NF1
# 7  ENSG00000196712 ENST00000356175.7 134.507 ENSG00000196712.20    NF1
# 8  ENSG00000196712 ENST00000358273.9  15.287 ENSG00000196712.20    NF1
# 9  ENSG00000196712 ENST00000431387.8  28.001 ENSG00000196712.20    NF1
# 10 ENSG00000196712 ENST00000487476.5 162.406 ENSG00000196712.20    NF1
# > txi_counts_filtered |> dplyr::filter(SYMBOL =="BSG")
# genes_truncated             tx_id    count            gene_id SYMBOL
# 1  ENSG00000172270 ENST00000679472.1   13.969 ENSG00000172270.22    BSG
# 2  ENSG00000172270 ENST00000333511.9    7.528 ENSG00000172270.22    BSG
# 3  ENSG00000172270 ENST00000353555.9 1277.376 ENSG00000172270.22    BSG
# 4  ENSG00000172270 ENST00000545507.6    2.082 ENSG00000172270.22    BSG
# 5  ENSG00000172270 ENST00000573784.6    2.246 ENSG00000172270.22    BSG
# 6  ENSG00000172270 ENST00000680552.1    1.015 ENSG00000172270.22    BSG
# 7  ENSG00000172270 ENST00000618006.4    4.081 ENSG00000172270.22    BSG
# 8  ENSG00000172270 ENST00000680326.1    0.000 ENSG00000172270.22    BSG
# 9  ENSG00000172270 ENST00000346916.9    0.000 ENSG00000172270.22    BSG
# 10 ENSG00000172270 ENST00000576984.3    0.000 ENSG00000172270.22    BSG
# 11 ENSG00000172270 ENST00000614867.2    1.369 ENSG00000172270.22    BSG
# 12 ENSG00000172270 ENST00000680065.1  388.328 ENSG00000172270.22    BSG

## sample 4 "DG-WT-hRGC"
# > txi_counts_filtered |> dplyr::filter(SYMBOL == "NF1")
# genes_truncated             tx_id   count            gene_id SYMBOL
# 1  ENSG00000196712 ENST00000687027.1   0.000 ENSG00000196712.20    NF1
# 2  ENSG00000196712 ENST00000684826.1   0.000 ENSG00000196712.20    NF1
# 3  ENSG00000196712 ENST00000693617.1  41.206 ENSG00000196712.20    NF1
# 4  ENSG00000196712 ENST00000696138.1  27.928 ENSG00000196712.20    NF1
# 5  ENSG00000196712 ENST00000691014.1  75.489 ENSG00000196712.20    NF1
# 6  ENSG00000196712 ENST00000490416.2   1.000 ENSG00000196712.20    NF1
# 7  ENSG00000196712 ENST00000356175.7 144.224 ENSG00000196712.20    NF1
# 8  ENSG00000196712 ENST00000358273.9   8.808 ENSG00000196712.20    NF1
# 9  ENSG00000196712 ENST00000431387.8   4.943 ENSG00000196712.20    NF1
# 10 ENSG00000196712 ENST00000487476.5  57.716 ENSG00000196712.20    NF1

# > txi_counts_filtered |> dplyr::filter(SYMBOL =="BSG")
# genes_truncated             tx_id   count            gene_id SYMBOL
# 1  ENSG00000172270 ENST00000679472.1   7.828 ENSG00000172270.22    BSG
# 2  ENSG00000172270 ENST00000333511.9   3.303 ENSG00000172270.22    BSG
# 3  ENSG00000172270 ENST00000353555.9 666.726 ENSG00000172270.22    BSG
# 4  ENSG00000172270 ENST00000545507.6   0.000 ENSG00000172270.22    BSG
# 5  ENSG00000172270 ENST00000573784.6   4.010 ENSG00000172270.22    BSG
# 6  ENSG00000172270 ENST00000680552.1   0.000 ENSG00000172270.22    BSG
# 7  ENSG00000172270 ENST00000618006.4   1.041 ENSG00000172270.22    BSG
# 8  ENSG00000172270 ENST00000680326.1   2.019 ENSG00000172270.22    BSG
# 9  ENSG00000172270 ENST00000346916.9   0.000 ENSG00000172270.22    BSG
# 10 ENSG00000172270 ENST00000576984.3   1.077 ENSG00000172270.22    BSG
# 11 ENSG00000172270 ENST00000614867.2   1.044 ENSG00000172270.22    BSG
# 12 ENSG00000172270 ENST00000680065.1 253.843 ENSG00000172270.22    BSG

