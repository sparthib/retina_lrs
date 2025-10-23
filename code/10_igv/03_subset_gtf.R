

library(GenomicFeatures)
library(rtracklayer)
library(stringr)

# Load the GTF file
gtf_dir <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly"
gtf_file <- file.path(gtf_dir, "release_46_primary_assembly.gtf")

# RO_gtf_file <- file.path(gtf_dir, "ROs_protein_coding_annotations.gtf")
# FT_vs_RGC_gtf_file <- file.path(gtf_dir, "FT_vs_RGC_protein_coding_annotations.gtf")
# 
# RO_gtf <- import(RO_gtf_file)
# FT_vs_RGC_gtf <- import(FT_vs_RGC_gtf_file)
# 
# # Subset to only include transcripts
transcripts_only <- gtf[gtf$type == "transcript" | 
                           gtf$type == "exon" | 
                           gtf$type == "gene"]
# transcripts_only$transcript_id <- sub("\\..*", "",
#                                       transcripts_only$transcript_id)
# 


# only keep these transcripts 
transcripts_of_interest <- c("ENST00000678899", "ENST00000616122", "ENST00000369622", 
                              "ENST00000355238", "ENST00000398125", "ENST00000264758", 
                              "ENST00000442998", "ENST00000374467", 
                              "ENST00000353555", "ENST00000333511",
                             "ENST00000395479", "ENST00000233596",
                             )

transcripts_of_interest_gtf <- gtf[
    str_remove(mcols(gtf)$transcript_id, "\\.\\d+$") %in% transcripts_of_interest
]

## check number of type == transcript


# only keep type = transcript if grep transcripts of interest
# SYNCRIP (ROs) has 4 major DTU isoforms: ENST00000678899, ENST00000616122, ENST00000369622, ENST00000355238
# For ADD1 (FT vs RGCs), I think there is 2 DTU transcripts: ENST00000398125, ENST00000264758
# For BAK1 (FT vs RGCs), there should only be 2 and both were DTU: ENST00000442998, ENST00000374467
# Last one for now I think may be BSG (ROs) also with 2 DTU isoforms:, ENST00000353555, ENST00000333511


#remove version numbers from transcript IDs in gtf 

# transcripts_of_interest_gtf <- transcripts_only[transcripts_only$transcript_id %in% transcripts_of_interest]

length(transcripts_of_interest_gtf)

# Save the filtered GTF file
output_gtf_file <- file.path(gtf_dir, "IGV_transcripts.gtf")
#SAVE to output_gtf_file
export(transcripts_of_interest_gtf, con = output_gtf_file, format = "gtf")



