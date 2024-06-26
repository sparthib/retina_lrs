#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=gene_body_cov
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/plot_gene_body_cov_transcript_length.txt
#SBATCH -e logs/plot_gene_body_cov_transcript_length.txt


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
# echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"

ref_gene_model=/dcs04/hicks/data/sparthib/references/rseqc_gencode_44_comp.bed
bam_paths=/users/sparthib/retina_lrs/code/coverage/H9-FT2_bam_paths.txt
path_to_output=/users/sparthib/retina_lrs/plots/coverage/H9_FT2_deciles
mkdir -p $path_to_output

ml load rseqc/3.0.1
geneBody_coverage.py -i $bam_paths -r $ref_gene_model -o $path_to_output

echo "**** Job ends ****"
date +"%Y-%m-%d %T"