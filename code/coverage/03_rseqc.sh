#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=30G
#SBATCH -c 10
#SBATCH --job-name=gene_body_cov
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/plot_gene_body_cov_transcript_length.%a.txt
#SBATCH -e logs/plot_gene_body_cov_transcript_length.%a.txt
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

ref_gene_model=/dcs04/hicks/data/sparthib/references/rseqc_gencode_44_comp.bed
bam_paths=/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/genome/isoquant/${sample}_deciles/${sample}_bam_paths.txt
touch $bam_paths
path_to_output=/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/genome/isoquant/${sample}_deciles/

cd /dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/genome/isoquant/${sample}_deciles/
rm $bam_paths
files_array=($(ls *.bam))
for file in "${files_array[@]}"; 
do 
    echo "/dcs04/hicks/data/sparthib/retina_lrs/05c_coverage/genome/isoquant/${sample}_deciles/$file" >> $bam_paths
done


ml load rseqc/3.0.1
geneBody_coverage.py -i $bam_paths -r $ref_gene_model -o $path_to_output

echo "**** Job ends ****"
date +"%Y-%m-%d %T"