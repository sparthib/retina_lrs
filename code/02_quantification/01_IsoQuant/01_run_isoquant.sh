#!/bin/bash

#SBATCH -p shared
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=isoquant
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/isoquant.%a.txt
#SBATCH -e logs/isoquant.%a.txt
#SBATCH --array=1-6,8,9,10,11,12,13,14,15
#SBATCH --time=7-00:00:00

#try running for all chromosomes

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


source activate isoquant 

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/primary_assembly/high_quality
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_assembly.gtf.gz
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/primary_assembly/release_46_primary_genome.fa.gz


OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/primary_assembly/${sample}
mkdir -p $OUTPUT_FOLDER

isoquant.py --reference $REFERENCE_FASTA --data_type ont --genedb $REFERENCE_GTF \
--bam ${BAM_FOLDER}/${sample}_primary_over_30_sorted.bam \
  --output $OUTPUT_FOLDER --clean_start --count_exons -t ${SLURM_CPUS_PER_TASK} --complete_genedb
  
  
conda deactivate

echo "**** Job ends ****"
date +"%Y-%m-%d %T"