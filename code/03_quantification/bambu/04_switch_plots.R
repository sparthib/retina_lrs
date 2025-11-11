library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtranscript)
library(biomaRt)

data_dir <- Sys.getenv("retina_lrs_dir")
code_dir <- Sys.getenv("retina_lrs_code")
ref_dir <- Sys.getenv("references_dir")

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

mart_attributes <- listAttributes(mart)
#write to file
write.table(mart_attributes, file = file.path(data_dir, "processed_data/mart_attributes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


gtf <- read.table(file.path(data_dir, "06_quantification/bambu/RGC_FT_extended_annotation/extended_annotations.gtf"),
                  header = FALSE, sep = "\t")
colnames(gtf) <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#split attribute column into multiple columns by semicolon
gtf <- gtf |> separate(attribute, into = c("gene_id", "transcript_id", "exon_number"), sep = ";") 

exon_gtf <- gtf |> filter(feature == "exon")
#remove "gene_id" and "transcript_id" from the attribute column
exon_gtf$gene_id <- gsub("gene_id ", "", exon_gtf$gene_id)
exon_gtf$transcript_id <- gsub("transcript_id ", "", exon_gtf$transcript_id)
exon_gtf$exon_number <- gsub("exon_number ", "", exon_gtf$exon_number)

#remove version number from gene and transcript id
exon_gtf$gene_id <- gsub("\\..*", "", exon_gtf$gene_id)
exon_gtf$transcript_id <- gsub("\\..*", "", exon_gtf$transcript_id)

#remove space from gene_id and transcript_id
exon_gtf$gene_id <- gsub(" ", "", exon_gtf$gene_id)
exon_gtf$transcript_id <- gsub(" ", "", exon_gtf$transcript_id)


# extract exons
switches_dir <- file.path(code_dir, "processed_data/dtu/IsoformSwitchAnalyzeR/bambu/FT_vs_RGC")
output_plot_dir <- file.path(code_dir,"plots/de/switch_analyzer/bambu/FT_vs_RGC")


#read switches.tsv file
switches <- read.table(paste0(switches_dir, "/FT_vs_RGCswitches.tsv"), header = TRUE, sep = "\t")
switches$transcript_id <- gsub("\\..*", "", switches$isoform_id)

top_switches_gtf <- exon_gtf |> filter(transcript_id %in% switches$transcript_id)
nrow(top_switches_gtf)



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
#for transcript id starting with "Bambu",replace transcript_biotype with "Bambu Novel"
top_switches_gtf$transcript_biotype <- ifelse(grepl("^Bambu", top_switches_gtf$transcript_id), "Bambu Novel", top_switches_gtf$transcript_biotype)
#replace NA with unknown biotype
top_switches_gtf$transcript_biotype <- ifelse(is.na(top_switches_gtf$transcript_biotype), "Unknown biotype", top_switches_gtf$transcript_biotype)
top_switches_gtf$gene_name <- ifelse(top_switches_gtf$gene_name == "", "Unknown name", top_switches_gtf$gene_name)


#generate pdf with plot for each gene in the top_genes_gtf

pdf(file.path(output_plot_dir, "switches.pdf"), width = 10, height = 10)

for (gene in unique(top_switches_gtf$gene_id)) {
  exons <- top_switches_gtf |> filter(gene_id == gene)
 p <-  exons |> ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
    geom_range(
      aes(fill = transcript_biotype)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id"),
      aes(strand = strand), arrow.min.intron.length = 200
    ) +
    ggtitle(paste(gene, exons$gene_name))
  print(p)
}
dev.off()
