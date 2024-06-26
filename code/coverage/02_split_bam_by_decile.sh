#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 5
#SBATCH --job-name=highqualbam
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/split_bam_by_decile.txt
#SBATCH -e logs/split_bam_by_decile.txt
#SBATCH -t 7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml load samtools 

# Define the folder containing the decile files
decile_folder=/dcs04/hicks/data/sparthib/retina_lrs/06_quantification/isoquant/FT_RGC/OUT/H9_FT2_deciles

# Define the input BAM file
input_bam=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9-FT_2_primary_over_30_chr_only_sorted.bam

# Loop over each decile file and process it
for i in {1..10}
do
    # Define the current decile file and output BAM file
    decile_file=$decile_folder/decile_${i}.tsv
    output_bam=$decile_folder/H9_FT2_genome_decile_${i}.bam

    # Subset the BAM file based on the current decile file
    samtools view -h $input_bam | \
    awk 'NR==FNR {reads[$1]; next} ($1 in reads) || ($1 ~ /^@/)' $decile_file - | \
    samtools view -Sb - > $output_bam
done

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
