
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

motif <- universalmotif::create_motif("GCATG", 
                      nsites = 50,
                      pseudocount = 1)


fimo_results <- runFimo(dna_set, motif, thresh = 1e-3)
class(fimo_results)
fimo_df <- fimo_results %>% as.data.frame()
write.csv(fimo_df, "/users/sparthib/retina_lrs/processed_data/dtu/fimo_results.csv", row.names = FALSE)

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



# tutorial
# data("example_chip_summits", package = "memes")
# 
# dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3
# 
# # Take 100bp windows around ChIP-seq summits
# summit_flank <- example_chip_summits %>% 
#   plyranges::anchor_center() %>% 
#   plyranges::mutate(width = 100) 
# 
# # Get sequences in peaks as Biostring::BStringSet
# sequences <- summit_flank %>% 
#   get_sequence(dm.genome)
# 
# names(sequences)[1:2]
# 
# 
# e93_motif <- MotifDb::MotifDb %>% 
#   # Query the database for the E93 motif using it's gene name
#   MotifDb::query("Eip93F") %>% 
#   # Convert from motifdb format to universalmotif format
#   universalmotif::convert_motifs() %>% 
#   # The result is a list, to simplify the object, return it as a single universalmotif
#   .[[1]]
# 
# e93_motif["name"] <- "E93_FlyFactor"
# 
# fimo_results <- runFimo(sequences, e93_motif)
# 
# fimo_results |> dplyr::select(pvalue, matched_sequence) |> unique()
# 
# fimo_results |> dplyr::select( matched_sequence) |> unique() 
# 




#### matrix ####

motif_probabilities <- as.matrix(rbind(
  c(0, 0, 1, 0, 0),
  c(0, 1, 0, 0, 0),
  c(1, 0, 0, 0, 1),
  c(0, 0, 0, 1, 0)
))

background_probabilities <- as.matrix(rbind(
  c(0.25, 0.25, 0.25, 0.25),
  c(0.25, 0.25, 0.25, 0.25),
  c(0.25, 0.25, 0.25, 0.25),
  c(0.25, 0.25, 0.25, 0.25)
))

log2(0.99^5 / 0.25^5)

# create a position weight matrix (PWM) for these sequences
# 
# ATTTCAGCGA
# ATATGGCGAA
# AGTTCAGCGA



