#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=plot_gc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/plot_gc/transcriptome/plot_gc.%a.txt
#SBATCH -e logs/plot_gc/transcriptome/plot_gc.%a.txt
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
ref_gene_model=/dcs04/hicks/data/sparthib/references/rseqc_gencode_44_comp.bed
BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/${sample}_chromosome_level/$sample.bam
output=/users/sparthib/retina_lrs/plots/gene_body_percentile/
mkdir $output

ml load rseqc/3.0.1

read_GC.py -i $BAM_FILE -r $ref_gene_model -o $output


conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"