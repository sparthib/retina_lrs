# To get started, the pipeline needs access to a gene annotation
# file in GFF3 or GTF format, a directory containing one or more FASTQ 
# files (which will be merged as pre-processing), a genome FASTA file,
# as well as the file path to minimap2, and the file path to the directory 
# to hold output files.
library("FLAMES")


gtf <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
genome_fasta <- "/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa"
outdir <- "/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flames/test"
genome_bam <- "/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/H9-hRGC_1_sorted.bam"

config_file <- FLAMES::create_config(outdir, do_barcode_demultiplex = FALSE)
config <- jsonlite::fromJSON(config_file)

# variants <- find_variants(
#   bam_path = genome_bam,
#   reference = genome_fasta,
#   annotation = gtf,
#   min_nucleotide_depth = 100,
#   homopolymer_window = 3,
#   annotated_region_only = FALSE,
#   names_from = "gene_name",
#   threads = 1
# )
# 


isoforms <- find_isoform(annotation = gtf, genome_fa = genome_fasta,
  genome_bam = genome_bam,
  outdir = outdir, config = config)
  