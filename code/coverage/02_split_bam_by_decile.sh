#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 5
#SBATCH --job-name=split_by_decile
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/split_bam_by_decile.txt
#SBATCH -e logs/split_bam_by_decile.txt
#SBATCH --array=1-15
#SBATCH --time=3-00:00:00

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
BAM_FOLDER=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/
sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo "$sample"


ml load samtools 


# Define the folder containing the decile files
decile_folder=/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/genome/isoquant/${sample}_deciles
mkdir -p $decile_folder

# Define the input BAM file
input_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/${sample}_primary_over_30_chr_only_sorted.bam

# Loop over each decile file and process it
for i in {1..10}
do
    # Define the current decile file and output BAM file
    decile_file=$decile_folder/decile_${i}.tsv
    output_bam=$decile_folder/$decile_${i}_$sample.bam

    # Subset the BAM file based on the current decile file
    samtools view -h $input_bam | \
    awk 'NR==FNR {reads[$1]; next} ($1 in reads) || ($1 ~ /^@/)' $decile_file - | \
    samtools view -Sb - > $output_bam
    samtools index $output_bam > ${output_bam}.bai
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
