
module load singularity


singularity pull docker://lingminhao/bambusc:beta1.2


reference_fasta=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
reference_gtf=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf

nextflow run GoekeLab/bambu-singlecell-spatial \
  --reads $PWD/examples/reads_chr9_1_1000000.fastq.gz \
  --genome $reference_dir \
  --annotation $reference_gtf \
  --chemistry 10x3v2 \
  --ncore 16 --outdir output \
  -with-singularity lingminhao/bambusc:beta1.2
