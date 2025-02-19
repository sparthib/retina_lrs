#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=bambu_clump
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/bambu_clump.txt
#SBATCH -e logs/bambu_clump.txt
#SBATCH --time=7-00:00:00


ml load singularity
ml load nextflow

reference_fasta=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa
reference_gtf=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf
output_dir=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/05c_bambu_clump/

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

singularity pull docker://lingminhao/bambusc:beta1.2
nextflow run GoekeLab/bambu-singlecell-spatial \
  -r main \
  --bams /users/sparthib/retina_lrs/raw_data/single_cell_samples.csv \
  --genome $reference_dir \
  --annotation $reference_gtf \
  --chemistry 10x3v2 \
  --ncore 19 --outdir $output_dir \
  -with-singularity lingminhao/bambusc:beta1.2 
  
echo "****Job Ends****"
date +"%Y-%m-%d %T"




