#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=gene_body_cov
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/plot_gene_body_cov/length.txt
#SBATCH -e logs/plot_gene_body_cov/length.txt
#SBATCH --time=7-00:00:00



echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
# echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"


# CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
# sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
# echo $sample
ref_gene_model=/dcs04/hicks/data/sparthib/references/rseqc_gencode_44_comp.bed
SHORT_BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_short_transcripts.bam      
MEDIUM_BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_medium_transcripts.bam   
LONG_BAM_FILE=/dcs04/hicks/data/sparthib/retina_lrs/05_bams/genome/GENCODE_splice/primary_over_30_chr_only/H9_FT_1_long_transcripts.bam   
short_output=/users/sparthib/retina_lrs/plots/gene_body_percentile/short/
medium_output=/users/sparthib/retina_lrs/plots/gene_body_percentile/medium/
long_output=/users/sparthib/retina_lrs/plots/gene_body_percentile/long/
mkdir $short_output
mkdir $medium_output
mkdir $long_output

ml load rseqc/3.0.1

echo "**** Plotting gene body coverage for short length reads****"
geneBody_coverage.py -i $SHORT_BAM_FILE -r $ref_gene_model -o $short_output
echo "**** Plotting gene body coverage for medium length reads****"
geneBody_coverage.py -i $MEDIUM_BAM_FILE -r $ref_gene_model -o $medium_output
echo "**** Plotting gene body coverage for long length reads****"
geneBody_coverage.py -i $LONG_BAM_FILE -r $ref_gene_model -o $long_output


conda deactivate 

echo "**** Job ends ****"
date +"%Y-%m-%d %T"