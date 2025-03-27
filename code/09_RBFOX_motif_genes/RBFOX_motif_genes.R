library(Biostrings)
library(GenomicRanges)
library(rtracklayer)

motifs <- c("TGCATGT", "TGCATG", "GCATGT", "GCATG")  # DNA versions of (U)GCAUG(U)


fasta <- readDNAStringSet("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa")  # replace with your path
names(fasta) <- sapply(strsplit(names(fasta), " "), `[`, 1)  # clean up chromosome names if needed
fasta <- fasta[1:25]

motif_hits <- lapply(names(fasta), function(seqname) {
  seq <- fasta[[seqname]]
  hits_list <- lapply(motifs, function(motif) {
    match <- matchPattern(motif, seq)
    if (length(match) > 0) {
      GRanges(seqnames = seqname,
              ranges = ranges(match),
              strand = "*",
              motif = motif)
    } else {
      NULL
    }
  })
  do.call(c, hits_list)
})

motif_hits <- Filter(Negate(is.null), motif_hits)
motif_hits_flat <- unlist(motif_hits, recursive = FALSE)
saveRDS(motif_hits_flat, "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_hits_flat.rds")
motif_hits_flat <- readRDS("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_hits_flat.rds")

df_list <- lapply(motif_hits_flat, 
                  as.data.frame)

df <- dplyr::bind_rows(df_list)

saveRDS(df, "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_hits_df.rds")


### get gtf file 

gtf <- rtracklayer::import("/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly_protein_coding_lncRNA.gtf")
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df |> dplyr::filter(gene_type == "protein_coding")
gtf_df |> colnames()


# head(df)
# seqnames start   end width strand   motif
# 1     chr1 11812 11818     7      * TGCATGT
# 2     chr1 12139 12145     7      * TGCATGT
# 3     chr1 23173 23179     7      * TGCATGT
# 4     chr1 33992 33998     7      * TGCATGT
# 5     chr1 37286 37292     7      * TGCATGT
# 6     chr1 40486 40492     7      * TGCATGT

# head(gtf_df)
# seqnames start    end width strand source type score phase           gene_id
# 1     chr1 11869  14409  2541      + HAVANA gene    NA    NA ENSG00000290825.1
# 2     chr1 29554  31109  1556      + HAVANA gene    NA    NA ENSG00000243485.5
# 3     chr1 34554  36081  1528      - HAVANA gene    NA    NA ENSG00000237613.2
# 4     chr1 57598  64116  6519      + HAVANA gene    NA    NA ENSG00000290826.1
# 5     chr1 65419  71585  6167      + HAVANA gene    NA    NA ENSG00000186092.7
# 6     chr1 89295 133723 44429      - HAVANA gene    NA    NA ENSG00000238009.6

library(data.table)

# Convert data.frames to data.tables
setDT(df)
setDT(gtf_df)

# Rename columns to match foverlaps convention
setnames(df, c("start", "end"), c("query_start", "query_end"))
setnames(gtf_df, c("start", "end"), c("subject_start", "subject_end"))

# Add keys for overlap
setkey(gtf_df, seqnames, subject_start, subject_end)
setkey(df, seqnames, query_start, query_end)

# Perform the overlap join
result <- foverlaps(df, gtf_df, by.x = c("seqnames", "query_start", "query_end"),
                    by.y = c("seqnames", "subject_start", "subject_end"),
                    type = "within", nomatch = 0)

# If needed, restore original column names
setnames(result, c("query_start", "query_end", "subject_start", "subject_end"),
         c("motif_start", "motif_end", "gene_start", "gene_end"))


result_df <- as.data.frame(result)
result_df <- result_df |> dplyr::select(seqnames, motif_start, 
                                        motif_end, motif, gene_id, 
                                        gene_start, gene_end,
                                        gene_name, strand)

nrow(result_df)

readr::write_csv(result_df, 
          "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/motif_gene_overlap.csv")





