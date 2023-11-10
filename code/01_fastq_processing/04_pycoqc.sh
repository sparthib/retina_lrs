#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=80G
#SBATCH --job-name=pycoqc
#SBATCH -c 10
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/pycoqc.%a.txt
#SBATCH -e logs/pycoqc.%a.txt
#SBATCH --array=1-4

echo "**** Job starts ****"
date +"%Y-%m-%d %T"
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"


##sequencing_summary_PAQ25549_580ccca1_aa394d34.txt

CONFIG=/users/sparthib/retina_lrs/raw_data/data_paths.config
sample=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $2}' $CONFIG)
guppy_summary_file=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print $4}' $CONFIG)
PYCOQC_OUTPUT=/dcs04/hicks/data/sparthib/casey/pyqoqc_outputs

mkdir -p $PYCOQC_OUTPUT

source activate pycoqc

cd /users/sparthib/pycoQC/pycoQC
./pycoQC.py –f $guppy_summary_file –o $PYCOQC_OUTPUT/${sample}.html


conda deactivate 


echo "**** Job ends ****"
date +"%Y-%m-%d %T"