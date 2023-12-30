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
#SBATCH --array=1-4

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
cd /users/sparthib/IsoQuant/

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo "${sample}"
BAM_FOLDER=/dcs04/hicks/data/sparthib/casey/bams
REFERENCE_GTF=/dcs04/hicks/data/sparthib/GENCODE_GTF.gtf.gz 
REFERENCE_FASTA=/dcs04/hicks/data/sparthib/GENCODE_FASTA.fa.gz 


OUTPUT_FOLDER=/dcs04/hicks/data/sparthib/casey/IsoQuant_output/${sample}
rm -r $OUTPUT_FOLDER
mkdir -p $OUTPUT_FOLDER

isoquant.py --reference $REFERENCE_FASTA --data_type ont --genedb $REFERENCE_GTF --bam ${BAM_FOLDER}/${sample}_sorted.bam \
  --output $OUTPUT_FOLDER --count_exons --clean_start -t ${SLURM_CPUS_PER_TASK} --complete_genedb
  
  
conda deactivate

echo "**** Job ends ****"
jobstats
date +"%Y-%m-%d %T"