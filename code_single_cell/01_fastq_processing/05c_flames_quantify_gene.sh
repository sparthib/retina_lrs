#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=flames
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/flames/flames.%a.txt
#SBATCH -e logs/flames/flames.%a.txt
#SBATCH --array=7
#SBATCH --time=7-00:00:00

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
num_cells=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $5}' $CONFIG)
path=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
seq_sum=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
echo $sample

#RUNNING JUST QUANTIFY GENE FUNCTION FOR SAMPLE 7

#create symlink for bam file 
# ln -s /dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/${sample}_sorted.bam /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10X_D200-EP1-1_B1/align2genome.bam
# ln -s /dcs04/hicks/data/sparthib/retina_single_cell_lrs/04_minimap2_output/genome/bams/${sample}_sorted.bam.bai /dcs04/hicks/data/sparthib/retina_single_cell_lrs/05_flames_output/10X_D200-EP1-1_B1/align2genome.bam.bai

ml load conda_R/4.4.x
Rscript 05c_flames_quantify_gene.R "${sample}" 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"