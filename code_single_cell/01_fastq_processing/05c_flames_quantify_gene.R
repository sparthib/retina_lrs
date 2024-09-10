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

#quantify gene is an internal function
FLAMES:::quantify_gene(
  annotation = annot,
  outdir,
  infq = fastq,
  n_process = 15, #reduced to less than allotted to avoid memory issues
  pipeline = "sc_single_sample",
  samples = NULL)

print("done quantifying genes")

sessionInfo()
