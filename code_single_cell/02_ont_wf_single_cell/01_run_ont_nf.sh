#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=100G
#SBATCH --job-name=ont_sc_wf
#SBATCH -c 20
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o /dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_ont_wf_single_cell_output/logs/ont_sc_wf.%a.txt
#SBATCH -e /dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_ont_wf_single_cell_output/logs/ont_sc_wf.%a.txt
#SBATCH --array=1

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

ml load nextflow
ml load singularity

CONFIG=/users/sparthib/retina_lrs/raw_data/single_cell.config
INPUT_FOLDER=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $3}' $CONFIG)
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
echo $sample
FASTQ=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/01_input_fastqs/${sample}.fastq.gz
output_folder=/dcs04/hicks/data/sparthib/retina_single_cell_lrs/02_ont_wf_single_cell_output/

mkdir -p $output_folder/logs/$sample

nextflow -log $output_folder/logs/$sample run epi2me-labs/wf-single-cell \
    --fastq $FASTQ \
    --expected_cells 100 \
    --kit_name '3prime' \
    --kit_version 'v3' \
    --ref_genome_dir '/dcs04/hicks/data/sparthib/references/genome/GENCODE/' \
    -profile singularity \
    --out_dir $output_folder

echo "**** Job ends ****"
date +"%Y-%m-%d %T"