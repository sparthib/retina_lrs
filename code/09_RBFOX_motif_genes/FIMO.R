
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(universalmotif))
library(memes)
library(Biostrings)


hm_genome <-  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  
class(hm_genome)
# "BSgenome"

## convert to Biostrings::XStringSet
chroms <- seqnames(hm_genome)
chroms <- chroms[1:25]

dna_set <- Biostrings::getSeq(hm_genome, names = chroms)


# https://rnasysu.com/encori/motif_short/hg38/SBDH1753/homerResults/motif1.motif


# Create the 6x4 matrix
mat <- matrix(c(
  0.001, 0.050, 0.001, 0.948,
  0.001, 0.001, 0.997, 0.001,
  0.001, 0.997, 0.001, 0.001,
  0.997, 0.001, 0.001, 0.001,
  0.001, 0.001, 0.001, 0.997,
  0.001, 0.001, 0.997, 0.001
), nrow = 6, ncol = 4, byrow = TRUE)

# Transpose the matrix
t_mat <- t(mat)

# Print the transposed matrix
print(t_mat)

motif <- universalmotif::create_motif(t_mat, 
                                      alphabet = "DNA",
                                      type = "PPM")


fimo_results <- runFimo(dna_set, motif, thresh = 1e-3)
class(fimo_results)
fimo_df <- fimo_results %>% as.data.frame()
fimo_df <- fimo_df |> dplyr::filter(matched_sequence == "TGCATG")
nrow(fimo_df)
### grep for TGCATG in attribute column
fimo_df$matched_sequence <- grepl("TGCATG", fimo_df$attribute)
write.csv(fimo_df, "/users/sparthib/retina_lrs/processed_data/dtu/fimo_results.csv",
          row.names = FALSE)
fimo_df <- readr::read_csv("/users/sparthib/retina_lrs/processed_data/dtu/fimo_results.csv")

### save genomic ranges as gtf file
library(rtracklayer)

# Add required GTF fields
mcols(fimo_results)$source <- "motif_scan"            # or any descriptive source
mcols(fimo_results)$type <- "motif"                   # feature type
mcols(fimo_results)$score <- as.character(mcols(fimo_results)$score)  # GTF needs score as character
mcols(fimo_results)$frame <- "."                      # no reading frame for motifs

# Construct the attribute field
mcols(fimo_results)$attribute <- paste0(
  'motif_id "', mcols(fimo_results)$motif_id, '"; ',
  'motif_alt_id "', mcols(fimo_results)$motif_alt_id, '"; ',
  'pvalue "', mcols(fimo_results)$pvalue, '"; ',
  'qvalue "', mcols(fimo_results)$qvalue, '"; ',
  'matched_sequence "', mcols(fimo_results)$matched_sequence, '";'
)

# Keep only GTF-relevant metadata columns
mcols(fimo_results) <- mcols(fimo_results)[, c("source", "type", "score", "frame", "attribute")]

# Export to GTF
export(fimo_results, "/users/sparthib/retina_lrs/processed_data/dtu/fimo_motifs.gtf", format = "gtf")




##### add gene ID info ####

library(data.table)

gtf <- rtracklayer::import("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly_protein_coding_lncRNA.gtf")
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df |> dplyr::filter(gene_type == "protein_coding")
gtf_df |> colnames()


# Convert data.frames to data.tables
setDT(fimo_df)
setDT(gtf_df)

# Rename columns to match foverlaps convention
setnames(fimo_df, c("start", "end"), c("motif_start", "motif_end"))
setnames(gtf_df, c("start", "end"), c("gene_start", "gene_end"))

# Add keys for overlap
setkey(gtf_df, seqnames, gene_start, gene_end)
setkey(fimo_df, seqnames, motif_start, motif_end)

colnames(fimo_df) <- c("seqnames", "motif_start", "motif_end", "width",
                       "motif_strand", "motif_id", "motif_alt_id","score",
                       "pvalue", "qvalue", "matched_sequence")

gtf_df <- gtf_df |> dplyr::select(seqnames, gene_start, gene_end, 
                                   gene_id, gene_name)
# Perform the overlap join
result <- foverlaps(fimo_df, gtf_df, by.x = c("seqnames", "motif_start", "motif_end"),
                    by.y = c("seqnames", "gene_start", "gene_end"),
                    type = "within", nomatch = 0)

# If needed, restore original column names
setnames(result, c("query_start", "query_end", "subject_start", "subject_end"),
         c("motif_start", "motif_end", "gene_start", "gene_end"))


result_df <- as.data.frame(result)
result_df <- result_df |> dplyr::select(seqnames, motif_start, score,
                                        motif_end,  gene_id, 
                                        gene_start, gene_end,
                                        gene_name)
result_df <- result_df |> dplyr::filter(score == 11.87)
nrow(result_df)
# > nrow(result_df)
# [1] 816394
readr::write_csv(result_df, 
                 "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_gene_overlap.csv")




