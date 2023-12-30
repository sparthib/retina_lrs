BiocManager::install("wiggleplotr")

library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")

# https://biodatascience.github.io/compbio/bioc/GRL.html

plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE)

suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

columns(edb)
# [1] "ENTREZID"            "EXONID"             
# [3] "EXONIDX"             "EXONSEQEND"         
# [5] "EXONSEQSTART"        "GENEBIOTYPE"        
# [7] "GENEID"              "GENENAME"           
# [9] "GENESEQEND"          "GENESEQSTART"       
# [11] "INTERPROACCESSION"   "ISCIRCULAR"         
# [13] "PROTDOMEND"          "PROTDOMSTART"       
# [15] "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE"
# [17] "PROTEINID"           "PROTEINSEQUENCE"    
# [19] "SEQCOORDSYSTEM"      "SEQLENGTH"          
# [21] "SEQNAME"             "SEQSTRAND"          
# [23] "SYMBOL"              "TXBIOTYPE"          
# [25] "TXCDSSEQEND"         "TXCDSSEQSTART"      
# [27] "TXID"                "TXNAME"             
# [29] "TXSEQEND"            "TXSEQSTART"         
# [31] "UNIPROTDB"           "UNIPROTID"          
# [33] "UNIPROTMAPPINGTYPE" 


txdf <- select(edb,
               keys=keys(edb, "GENEID"),
               columns=c("GENEID","TXID"),
               keytype="GENEID")

ebt <- exonsBy(edb, by="tx")

class(ebt)

#TRY ebt[1] vs ebt[[1]]

# Let’s add the gene identifiers:

table(names(ebt) %in% txdf$TXID)

mcols(ebt)$GENEID <- mapIds(edb, names(ebt), "GENEID", "TXID")
mcols(ebt)

sapply(ebt[1:5], length)


mcols(ebt)$num_exons <- lengths(ebt)
mcols(ebt)