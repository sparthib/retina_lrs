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
                                     do_genome_alignment = FALSE,
                                     do_gene_quantification = FALSE,
                                     min_sup_cnt = 5,
                                     threads = 20)
# genome_bam <- file.path("/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/align2genome.bam")

#run minimap2

if (!any(is.na(sys_which(c("minimap2", "k8"))))) {
  sce <- sc_long_pipeline(
    annotation = annot, fastq = fastq,
    genome_fa = genome_fa,
    outdir = outdir, 
    config_file = config_file,
    expect_cell_number = 63)

  
  #sce pipeline failed because it didn't find the matched_reads.fastq file
  #https://github.com/mritchielab/FLAMES/issues/30
  #work around is to create a symlink to the blaze processed fastq file 
  #ln -s /dcs04/hicks/data/sparthib/retina_single_cell_lrs/03_blaze_processed/raw/high_sensitivity/10x_D200-EP1-1_B2_matched_reads.fastq /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/matched_reads.fastq

}


sessionInfo()




