#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=add_barcode
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/add_barcode.%a.txt
#SBATCH -e logs/add_barcode.%a.txt
#SBATCH --array=7
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
num_cells=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $5}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

input_bam=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/${sample}_primary_over_30_chr_only_sorted.bam
output_bam=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/tagged/${sample}_deduped.bam
rm $output_bam
ml load python/3.10.13
python3 01_add_barcode_tags.py $input_bam $output_bam


#sort tagged file 
input_bam=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/tagged/${sample}_deduped.bam
output_bam=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/tagged/${sample}_deduped_sorted.bam
output_bai=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/primary_over_30_chr_only/tagged/${sample}_deduped_sorted.bam.bai
rm $output_bam
rm $output_bai

ml load samtools

samtools sort $input_bam -o $output_bam
# rm ${BAM_FOLDER}/${sample}.bam
samtools index $output_bam $output_bai



echo "**** Job ends ****"
date +"%Y-%m-%d %T"
