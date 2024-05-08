#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --job-name=flair
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flair.%a.txt
#SBATCH -e logs/flair.%a.txt
#SBATCH --array=1-15


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"



CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

REFERENCE_FASTA=/dcs04/hicks/data/sparthib/references/genome/GENCODE/GRCh38.p14.genome.fa
REFERENCE_GTF=/dcs04/hicks/data/sparthib/references/genome/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf

source activate flair 
 
### take the bed files created from bedtools using minimap2 bams and input the bed12 into flair correct 
bed12_file=/dcs04/hicks/data/sparthib/retina_lrs/05b_beds/genome/GENCODE_splice/${sample}.bed12

correction_output=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/flair2/correction_output/${sample}
mkdir -p $correction_output


flair correct -q $bed12_file -f $REFERENCE_GTF -g $REFERENCE_FASTA --output $correction_output --print_check


echo "**** Job ends ****"
date +"%Y-%m-%d %T"