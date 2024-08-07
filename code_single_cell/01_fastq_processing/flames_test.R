library("FLAMES")
library("sessioninfo")

#have all input files uncompressed (Related to this issue: https://github.com/mritchielab/FLAMES/issues/31)
sample <- commandArgs(trailingOnly = TRUE)
fastq <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/", sample, "matched_reads.fastq")
genome_fa <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa")
annot <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf")
minimap_path <- file.path("/users/sparthib/minimap2")


outdir <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/", sample)
config_file <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/config/config_file_3766646.json")
# config_file <- FLAMES::create_config(config_dir, type = "sc_3end",
#                                      do_barcode_demultiplex = FALSE,
#                                      do_genome_alignment = FALSE,
#                                      do_gene_quantification = TRUE,
#                                      do_isoform_identification = TRUE,
#                                      do_read_realignment = TRUE,
#                                      do_transcript_quantification = TRUE,
#                                      min_sup_cnt = 5,
#                                      threads = 20)
genome_bam <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/", paste0(sample, "_sorted.bam"))

#run minimap2

if (!any(is.na(sys_which(c("minimap2", "k8"))))) {
  sce <- sc_long_pipeline(
    annotation = annot, fastq = fastq,
    genome_fa = genome_fa,
    outdir = outdir, 
    config_file = config_file,
    expect_cell_number = 2000)
  
  sce <- create_sce_from_dir(outdir = outdir, annotation = annot)
  
}

##TIPS:
#if you run blaze yourself, create a symlink to the unzipped blaze processed fastq file in the flames output directory called matched_reads.fastq
#set bambu_isoform_quantification = TRUE in the config file (something to do with the bam or gtf being unordered.)
#replace "." or "None" values in annot and isoform annot files with "*" before running quantify_transcript


sessionInfo()




