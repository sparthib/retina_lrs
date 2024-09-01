#!/bin/bash

#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH -c 10
#SBATCH --job-name=strand
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/strand_log.%a.txt
#SBATCH -e logs/strand_log.%a.txt
#SBATCH --array=1-15
#SBATCH -t 4-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"

ml load samtools 
input_dir=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only
input_bam=${input_dir}/${sample}_primary_over_30_chr_only_sorted.bam
#get strand info for each read
# Count positive strand reads
pos_count=$(samtools view $input_bam | awk '{if(and($2, 16) == 0) {count++}} END {print count}')

# Count negative strand reads
neg_count=$(samtools view $input_bam | awk '{if(and($2, 16) != 0) {count++}} END {print count}')

# Output the results
echo "Positive strand reads: $pos_count"
echo "Negative strand reads: $neg_count"

echo "**** Job ends ****"
date +"%Y-%m-%d %T"


