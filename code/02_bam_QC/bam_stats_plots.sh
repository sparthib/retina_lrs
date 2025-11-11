#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=pomoxis_bam_qc
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/pomoxis/plots.%a.txt
#SBATCH -e logs/pomoxis/plots.%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --array=1-15


echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"

ENV_FILE="../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

sample=$(awk -v Index=${SLURM_ARRAY_TASK_ID} '$1==Index {print $2}' $CONFIG)
echo $sample

input_dir="$retina_lrs_dir/05_bams/genome/primary_assembly/pomoxis_stats"
output_dir="$retina_lrs_code/plots/bam_qc/pomoxis_plots/"

ml load python/3.10.13
python3 bam_stats_plots.py $sample $input_dir $output_dir

input_dir="$retina_lrs_dir/05_bams/genome/primary_assembly/high_quality/pomoxis_stats"
output_dir="$retina_lrs_code/plots/bam_qc/pomoxis_plots/high_quality"

mkdir -p $output_dir
python3 bam_stats_plots.py $sample $input_dir $output_dir

echo "**** Job ends ****"
date +"%Y-%m-%d %T"
