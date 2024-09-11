library("FLAMES")
library("sessioninfo")

#have all input files uncompressed (Related to this issue: https://github.com/mritchielab/FLAMES/issues/31)
sample <- commandArgs(trailingOnly = TRUE)
fastq <- paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/", sample, "/matched_reads.fastq")
genome_fa <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa")
annot <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf")
minimap_path <- file.path("/users/sparthib/minimap2")
# genome_bam <- paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/", 
#                      sample, "_primary_over_30_chr_only_sorted.bam")

outdir <- paste0("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/", sample)
config_dir <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/config")
config_file <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/config/full_pipeline.json")

if (!any(is.na(sys_which(c("minimap2", "k8"))))) {

  #UNCOMMENT IF YOU WANT TO RUN PIPELINE FROM SCRATCH. 
  #CHECK CONFIG FILE TO MAKE SURE YOU HAVE SET THE RIGHT PARAMETERS BEFORE STARTING THE PIPELINE.
  # sce <- sc_long_pipeline(
  #   annotation = annot, fastq = fastq,
  #   genome_fa = genome_fa,
  #   outdir = outdir,
  #   config_file = config_file,
  #   expect_cell_number = 2000)
  #saveRDS(sce, file = paste0(outdir, "/sce.rds"))
  
  #below is given all required files post isoform quantification
  #already exist in the output directory
  sce <- create_sce_from_dir(outdir, annotation = annot)
  saveRDS(sce, file = paste0(outdir, "/sce.rds"))
  
  
}
##TIPS:
#if you run blaze yourself, create a symlink to the unzipped blaze processed fastq file in the flames output directory called matched_reads.fastq
#set bambu_isoform_quantification = TRUE in the config file (something to do with the bam or gtf being unordered.)
#replace "." or "None" values in annot and isoform annot files with "*" before running quantify_transcript

sessionInfo()




