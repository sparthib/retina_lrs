#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH -c 5
#SBATCH --job-name=chr_distribution
#SBATCH --mail-user=sparthi1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH -o logs/chr_distribution.%a.txt
#SBATCH -e logs/chr_distribution.%a.txt
#SBATCH --array=1-15
#SBATCH --time=7-00:00:00


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
# for stem_cell in DG-WT-hRGC EP1-BRN3B_iPSC EP1-WT_iPSC H9-BRN3B_ESC H9-BRN3B_iPSC H9-CRX_ESC H9-CRX_iPSC YZ_iPSC
# do
#     mkdir -p /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/stem_cell_vcfs/${stem_cell}
# done
stem_cell=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $1}' $CONFIG)
sample=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $2}' $CONFIG)

cp /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/01_longshot_vcfs/${sample}/${sample}.vcf /dcs04/hicks/data/sparthib/retina_lrs/09_ASE/02_g2gtools/stem_cell_vcfs/${sample}.vcf
  


ml load bcftools/1.18




