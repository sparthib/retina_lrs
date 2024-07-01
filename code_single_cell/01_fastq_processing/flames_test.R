library("FLAMES")
library("sessioninfo")

## blaze FASTQ files 


fastq <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_blaze_processed/raw/high_sensitivity/10x_D200-EP1-1_B2_matched_reads.fastq.gz")
genome_fa <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa")
annot <- file.path("/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation-short.gtf")
minimap_path <- file.path("/users/sparthib/minimap2")

outdir <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output")
config_dir <- file.path(outdir, "config")
config_file <- FLAMES::create_config(config_dir, type = "sc_3end",
                                     do_barcode_demultiplex = FALSE,
                                     min_sup_cnt = 5,
                                     threads = 20)

#run minimap2

if (!any(is.na(sys_which(c("minimap2", "k8"))))) {

  config <- jsonlite::fromJSON(config_file)
  genome_bam <- rownames(minimap2_align(
    config = config, fa_file = genome_fa, fq_in = fastq, annot = annot,
    outdir = outdir
  ))
  print("completed genome alignment")
  
  find_isoform(
    annotation = annot, genome_fa = genome_fa,
    genome_bam = genome_bam, outdir = outdir, config = config
  )
  print("completed isoform identification")
  minimap2_realign(
    config = config, fq_in = fastq,
    outdir = outdir
  )
  print("created minimap2 realigned bam")
  quantify_transcript(annotation = annot, outdir = outdir, config = config)
  print("completed transcript quantification")
  sce <- create_sce_from_dir(outdir = outdir, annotation = annot)
  print("created SingleCellExperiment object")
  
}


sessionInfo()



# sce <- sc_long_pipeline(
#   annotation = annot, fastq = fastq, genome_fa = genome_fa,
#   outdir = outdir, config_file = config_file, expect_cell_number = 63)
