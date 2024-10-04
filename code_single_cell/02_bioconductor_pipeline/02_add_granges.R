library(scuttle)
library(data.table)
library(Matrix)
library(biomaRt)
library(rtracklayer)


sample_names <- c( "10x_D100-EP1_A1",
                   "10x_D100-EP1_A2",
                   "10x_D100-EP1_B1",
                   "10x_D100-EP1_B2",
                   "10x_D200-EP1-1_A1",
                   "10x_D200-EP1-1_A2",
                   "10x_D200-EP1-1_B1",
                   "10x_D200-EP1-1_B2",
                   "10x_D200-EP1-2_A1",
                   "10x_D200-EP1-2_A2",
                   "10x_D200-EP1-2_B1",
                   "10x_D200-EP1-2_B2")

output_dir <- "/dcs04/hicks/data/sparthib/retina_single_cell_lrs/06_sce_rds_files/"


require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

create_granges <- function(sample){ 
  sce <- readRDS(file.path(output_dir, paste0(sample,'_sce.rds')))

  }



mart_attributes <- listAttributes(mart)

# chromosome_name, start_position, end_position, transcript_length,
# percentage_gene_gc_content, gene_biotype, transcript_biotype, hgnc_symbol

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id",
               "external_gene_name",
               "transcript_biotype"),
  filter="ensembl_transcript_id",
  values=top_switches_gtf$transcript_id,
  uniqueRows=TRUE)

colnames(annotLookup) <- c("transcript_id", "gene_name", "transcript_biotype")
top_switches_gtf <- merge(top_switches_gtf, annotLookup,
                          by="transcript_id", all.x=TRUE)

