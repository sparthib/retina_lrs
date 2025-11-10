#!/bin/bash

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=minIONQC
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/02_minIONQC/minIONQC.%a.txt
#SBATCH -e logs/02_minIONQC/minIONQC.%a.txt
#SBATCH --array=1-15

ENV_FILE="../../.env"
if [ -f $ENV_FILE ]; then
    set -a
    source $ENV_FILE
    set +a
fi

sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
guppy_summary_file=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
minIONQC_OUTPUT=$retina_lrs_dir/02_MinIONQC/${sample}

mkdir -p $minIONQC_OUTPUT


module load conda_R/4.3

Rscript 02_MinIONQC.R -i $guppy_summary_file  -o $minIONQC_OUTPUT \
  -p ${SLURM_CPUS_PER_TASK} -q 10


echo "**** Job ends ****"
date +"%Y-%m-%d %T"